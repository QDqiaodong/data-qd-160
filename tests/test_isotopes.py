"""Exact checks for the reaction-isotopes joint mass/tag accounting.

Small cases are verified against an *independent per-position* oracle that
enumerates every atom's isotope choice (with multiplicity, so two atoms taking
one light and one heavy isotope contribute both orders) and accumulates the
exact joint probability of ``(mass, tag count)``.  The oracle never imports the
production DP, never uses a recurrence and never reuses the production channel
assignment; it classifies every enumerated outcome by hand against the input
channels.  Scale cases (80 atoms) cannot be enumerated and are checked against
the documented closed-form expectations and exact conservation identities.
"""

import hashlib
import itertools
import json
import os
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction
from math import comb
from pathlib import Path

from reaction_batch.isotopes import analyze

ROOT = Path(__file__).resolve().parents[1]


def element(token, isotopes):
    return {'element': token, 'isotopes': isotopes}


def iso(mass_number, abundance, tagged):
    return {'mass_number': mass_number,
            'abundance': abundance, 'tagged': tagged}


def channel(channel_id, mass, tags):
    return {'channel_id': channel_id,
            'mass': {'lower': mass[0], 'upper': mass[1]},
            'tag_count': {'lower': tags[0], 'upper': tags[1]}}


def query(query_id, formula_text, cutoff, channels):
    return {'query_id': query_id, 'formula': formula_text,
            'cutoff': cutoff, 'channels': channels}


def a2_table():
    return [element('A', [iso(10, '0.5', False), iso(11, '0.5', True)])]


def enumerate_joint(elements, formula_text, cutoff):
    """Enumerate all atom-position assignments into exact (mass, tag) weights.

    Returns ``(states, base)`` where ``states`` maps ``(mass, tags)`` to a
    ``Fraction`` probability summed over every distinct assignment, including
    the multiplicity of atom positions and repeated isotope entries.
    """
    from reaction_batch.chemistry import formula
    counts = formula(formula_text)
    table = {entry['element']: entry['isotopes'] for entry in elements}
    base = sum(number * min(item['mass_number']
                            for item in table[token])
               for token, number in counts.items())
    atoms = []
    for token, number in counts.items():
        for _ in range(number):
            atoms.append([(item['mass_number'],
                           item['tagged'],
                           Fraction(item['abundance']))
                          for item in table[token]])
    states = {}
    for choices in itertools.product(*atoms):
        mass = sum(c[0] for c in choices)
        tags = sum(1 for c in choices if c[1])
        probability = Fraction(1)
        for c in choices:
            probability *= c[2]
        if 0 <= mass - base <= cutoff:
            key = (mass, tags)
            states[key] = states.get(key, Fraction(0)) + probability
    return states, base


def moments(states):
    """``(p, u, v)`` for a ``{(mass, tags): probability}`` mapping."""
    probability = sum(states.values(), Fraction(0))
    first = sum((tags * p for (_m, tags), p in states.items()), Fraction(0))
    second = sum((tags * (tags - 1) * p
                  for (_m, tags), p in states.items()), Fraction(0))
    return probability, first, second


def stats(p, u, v):
    if p == 0:
        return {'probability': '0', 'first_moment': '0',
                'second_factorial_moment': '0',
                'mean_tags': None, 'variance_tags': None}
    mean = u / p
    variance = (v + u) / p - mean * mean
    return {'probability': str(p), 'first_moment': str(u),
            'second_factorial_moment': str(v),
            'mean_tags': str(mean), 'variance_tags': str(variance)}


def oracle_query(elements, one_query):
    """Full oracle report for one query, independent of the production code."""
    cutoff = one_query['cutoff']
    states, base = enumerate_joint(
        elements, one_query['formula'], cutoff)
    channels = one_query['channels']

    def assign(mass, tags):
        result = [Fraction(0)] * len(channels)
        firsts = [Fraction(0)] * len(channels)
        seconds = [Fraction(0)] * len(channels)
        for i, ch in enumerate(channels):
            if (ch['mass']['lower'] <= mass <= ch['mass']['upper']
                    and ch['tag_count']['lower'] <= tags
                    <= ch['tag_count']['upper']):
                result[i] += 1
                firsts[i] += tags
                seconds[i] += tags * (tags - 1)
                break
        return result, firsts, seconds

    probability = [Fraction(0)] * len(channels)
    firsts = [Fraction(0)] * len(channels)
    seconds = [Fraction(0)] * len(channels)
    envelope = {}
    for (mass, tags), p in states.items():
        envelope[mass] = envelope.get(mass, Fraction(0)) + p
        hits, f, s = assign(mass, tags)
        for i in range(len(channels)):
            probability[i] += hits[i] * p
            firsts[i] += f[i] * p
            seconds[i] += s[i] * p

    # Untruncated total, itself computed from the full assignment enumeration
    # rather than any marginal formula.
    from reaction_batch.chemistry import formula
    counts = formula(one_query['formula'])
    table = {entry['element']: entry['isotopes'] for entry in elements}
    atoms = []
    for token, number in counts.items():
        for _ in range(number):
            atoms.append([(item['tagged'], Fraction(item['abundance']))
                          for item in table[token]])
    total_p = Fraction(0)
    total_u = Fraction(0)
    total_v = Fraction(0)
    for choices in itertools.product(*atoms):
        prob = Fraction(1)
        for _t, p in choices:
            prob *= p
        tags = sum(1 for t, _p in choices if t)
        total_p += prob
        total_u += tags * prob
        total_v += tags * (tags - 1) * prob
    assert total_p == 1

    kept = sum(envelope.values(), Fraction(0))
    kept_u = sum((tags * p for (_m, tags), p in states.items()),
                 Fraction(0))
    kept_v = sum((tags * (tags - 1) * p
                  for (_m, tags), p in states.items()), Fraction(0))
    overflow = stats(1 - kept, total_u - kept_u, total_v - kept_v)

    assigned_p = sum(probability, Fraction(0))
    assigned_u = sum(firsts, Fraction(0))
    assigned_v = sum(seconds, Fraction(0))
    missed = stats(total_p - assigned_p, total_u - assigned_u,
                   total_v - assigned_v)

    return {
        'query_id': one_query['query_id'],
        'base_mass': base,
        'envelope': [{'mass': mass, 'probability': str(envelope[mass])}
                     for mass in sorted(envelope) if envelope[mass]],
        'overflow': overflow,
        'channels': [{'channel_id': ch['channel_id'],
                      'stats': stats(probability[i], firsts[i], seconds[i])}
                     for i, ch in enumerate(channels)],
        'missed': missed,
        'total': stats(total_p, total_u, total_v),
    }


def by_id(entries):
    return {entry['query_id']: entry for entry in entries['queries']}


def check_oracle(test, data):
    report = analyze(data)
    test.assertEqual([q['query_id'] for q in report['queries']],
                     [q['query_id'] for q in data['queries']])
    for produced, one_query in zip(report['queries'], data['queries']):
        expected = oracle_query(data['elements'], one_query)
        test.assertEqual(produced, expected)


class NumericAcceptance(unittest.TestCase):
    def test_a2_joint_distribution_and_total_moments(self):
        data = {'elements': a2_table(),
                'queries': [query('q', 'A2', 2, [])]}
        out = analyze(data)['queries'][0]
        self.assertEqual(out['base_mass'], 20)
        self.assertEqual(out['envelope'],
                         [{'mass': 20, 'probability': '1/4'},
                          {'mass': 21, 'probability': '1/2'},
                          {'mass': 22, 'probability': '1/4'}])
        self.assertEqual(out['total']['first_moment'], '1')
        self.assertEqual(out['total']['second_factorial_moment'], '1/2')
        self.assertEqual(out['total']['mean_tags'], '1')
        self.assertEqual(out['total']['variance_tags'], '1/2')
        self.assertEqual(out['missed'], out['total'])
        self.assertEqual(out['overflow']['probability'], '0')
        self.assertIsNone(out['overflow']['mean_tags'])

    def test_marginals_positive_but_joint_channel_zero(self):
        data = {'elements': a2_table(),
                'queries': [query('q', 'A2', 2,
                                  [channel('gate', (21, 21), (0, 0))])]}
        out = analyze(data)['queries'][0]
        # Mass 21 has probability 1/2 and tag 0 has probability 1/4, yet no
        # outcome is both: (21,0) never occurs.
        mass_marginal = Fraction(1, 2)
        tag_marginal = Fraction(1, 4)
        self.assertGreater(mass_marginal, 0)
        self.assertGreater(tag_marginal, 0)
        entry = out['channels'][0]['stats']
        self.assertEqual(entry['probability'], '0')
        self.assertIsNone(entry['mean_tags'])
        self.assertIsNone(entry['variance_tags'])
        self.assertEqual(out['missed'], out['total'])

    def test_linkage_prefix_shrink_changes_successor_conditionals(self):
        channels = [channel('c1', (20, 21), (0, 2)),
                    channel('c2', (20, 22), (0, 2))]
        base = {'elements': a2_table(),
                'queries': [query('q', 'A2', 2, channels)]}
        first = by_id(analyze(base))['q']
        self.assertEqual(first['channels'][1]['stats'],
                         {'probability': '1/4', 'first_moment': '1/2',
                          'second_factorial_moment': '1/2',
                          'mean_tags': '2', 'variance_tags': '0'})
        narrowed = [channel('c1', (20, 20), (0, 2)),
                    channel('c2', (20, 22), (0, 2))]
        second = by_id(analyze({'elements': a2_table(),
                                'queries': [query('q', 'A2', 2, narrowed)]}))['q']
        # c2 itself is unchanged but now also receives (21,1).
        self.assertEqual(second['channels'][1]['stats'],
                         {'probability': '3/4', 'first_moment': '1',
                          'second_factorial_moment': '1/2',
                          'mean_tags': '4/3', 'variance_tags': '2/9'})
        self.assertEqual(first['missed']['probability'], '0')
        self.assertEqual(second['missed']['probability'], '0')

    def test_channel_order_defines_priority_not_identifier(self):
        channels = [channel('z', (20, 22), (0, 2)),
                    channel('a', (20, 22), (0, 2))]
        out = by_id(analyze({'elements': a2_table(),
                             'queries': [query('q', 'A2', 2, channels)]}))['q']
        # The first array entry takes everything; sorting by channel_id would
        # hand the mass to 'a' instead.
        self.assertEqual(out['channels'][0]['channel_id'], 'z')
        self.assertEqual(out['channels'][0]['stats']['probability'], '1')
        self.assertEqual(out['channels'][1]['stats']['probability'], '0')

    def test_cutoff_widening_leaves_covered_channels_identical(self):
        channels = [channel('c1', (20, 20), (0, 2)),
                    channel('c2', (21, 21), (0, 2))]
        narrow = by_id(analyze({'elements': a2_table(),
                                'queries': [query('q', 'A2', 1, channels)]}))['q']
        wide = by_id(analyze({'elements': a2_table(),
                              'queries': [query('q', 'A2', 2, channels)]}))['q']
        self.assertEqual(narrow['channels'], wide['channels'])
        self.assertEqual(narrow['missed'], wide['missed'])
        self.assertEqual(narrow['overflow'],
                         {'probability': '1/4', 'first_moment': '1/2',
                          'second_factorial_moment': '1/2',
                          'mean_tags': '2', 'variance_tags': '0'})
        self.assertEqual(wide['overflow'],
                         {'probability': '0', 'first_moment': '0',
                          'second_factorial_moment': '0',
                          'mean_tags': None, 'variance_tags': None})

    def test_truncated_tail_carries_tags(self):
        data = {'elements': a2_table(),
                'queries': [query('q', 'A', 0, [])]}
        out = analyze(data)['queries'][0]
        self.assertEqual(out['envelope'],
                         [{'mass': 10, 'probability': '1/2'}])
        self.assertEqual(out['overflow'],
                         {'probability': '1/2', 'first_moment': '1/2',
                          'second_factorial_moment': '0',
                          'mean_tags': '1', 'variance_tags': '0'})

    def test_zero_probability_channel_is_retained(self):
        data = {'elements': a2_table(),
                'queries': [query('q', 'A2', 2,
                                  [channel('empty', (0, 19), (0, 2)),
                                   channel('full', (20, 22), (0, 2))])]}
        out = analyze(data)['queries'][0]
        self.assertEqual(len(out['channels']), 2)
        self.assertEqual(out['channels'][0]['stats']['probability'], '0')
        self.assertEqual(out['channels'][1]['stats']['probability'], '1')


class ScaleAcceptance(unittest.TestCase):
    def test_a80_binomial_envelope_and_channel_conservation(self):
        data = {'elements': a2_table(),
                'queries': [query('q', 'A80', 80, [])]}
        out = analyze(data)['queries'][0]
        self.assertEqual(len(out['envelope']), 81)
        for k in range(81):
            self.assertEqual(out['envelope'][k],
                             {'mass': 800 + k,
                              'probability': str(Fraction(comb(80, k), 2 ** 80))})
        self.assertEqual(out['total']['mean_tags'], '40')
        self.assertEqual(out['total']['variance_tags'], '20')
        self.assertEqual(out['total']['first_moment'], '40')
        self.assertEqual(out['overflow']['probability'], '0')

    def test_a80_selective_first_channel_and_conservation(self):
        channels = [channel('exact', (840, 840), (40, 40)),
                    channel('broad', (800, 880), (0, 80))]
        data = {'elements': a2_table(),
                'queries': [query('q', 'A80', 80, channels)]}
        out = analyze(data)['queries'][0]
        first = out['channels'][0]['stats']
        self.assertEqual(first['probability'], str(Fraction(comb(80, 40), 2 ** 80)))
        self.assertEqual(first['mean_tags'], '40')
        self.assertEqual(first['variance_tags'], '0')
        # Second channel gets every retained outcome from 800..880 except the
        # single mass 840 (which the first channel took).
        reference = Fraction(comb(80, 40), 2 ** 80)
        self.assertEqual(Fraction(out['channels'][1]['stats']['probability']),
                         1 - reference)
        self.assertEqual(out['missed']['probability'], '0')
        for key in ('probability', 'first_moment', 'second_factorial_moment'):
            total = sum(Fraction(item['stats'][key])
                        for item in out['channels'])
            total += Fraction(out['missed'][key])
            self.assertEqual(total, Fraction(out['total'][key]))
        self.assertEqual(out['envelope'][40]['mass'], 840)

    def test_multi_element_offset_differs_from_tag_count(self):
        elements = [
            element('A', [iso(10, '0.5', False), iso(12, '0.5', True)]),
            element('B', [iso(1, '1', True)]),
        ]
        data = {'elements': elements,
                'queries': [query('q', 'A2B', 4, [])]}
        check_oracle(self, data)
        out = analyze(data)['queries'][0]
        # T is not S: the two B atoms are always tagged, so offset 0 (both A
        # light) still carries T = 2, and the mass/tag pairs are correlated.
        self.assertEqual(out['base_mass'], 21)
        joint = {(e['mass'], e['probability']) for e in out['envelope']}
        self.assertIn((21, '1/4'), joint)
        self.assertEqual(out['total']['mean_tags'], '2')
        self.assertEqual(out['envelope'][0], {'mass': 21, 'probability': '1/4'})


class OracleCrossCheck(unittest.TestCase):
    def test_nested_formula_against_position_oracle(self):
        # Chemistry.quantity only accepts decimal abundances.
        elements = [
            element('A', [iso(10, '0.333333333', False),
                          iso(11, '0.666666667', True)]),
            element('B', [iso(20, '0.75', False), iso(22, '0.25', True)]),
        ]
        data = {'elements': elements,
                'queries': [query('q', '(AB2)2', 5,
                                  [channel('c1', (80, 83), (1, 3)),
                                   channel('c2', (80, 90), (0, 4))])]}
        check_oracle(self, data)

    def test_random_small_queries_match_oracle(self):
        import random
        random.seed(20261118)
        splits = [['1'], ['0.5', '0.5'], ['0.25', '0.75'], ['0.2', '0.8'],
                  ['0.5', '0.25', '0.25'], ['0.1', '0.2', '0.7']]
        checked = 0
        for _ in range(40):
            elements = []
            for e_index in range(random.randint(1, 2)):
                token = 'AB'[e_index]
                abundances = random.choice(splits)
                count = len(abundances)
                masses = random.sample(range(1, 12), count)
                isotopes = [iso(masses[k], abundances[k], bool(random.getrandbits(1)))
                            for k in range(count)]
                elements.append(element(token, isotopes))
            pieces = []
            for _p in range(random.randint(1, 2)):
                token = random.choice([e['element'] for e in elements])
                pieces.append(token)
            if not pieces:
                pieces = [elements[0]['element']]
            text = ''.join(token + (str(random.randint(1, 3)) if random.random() < 0.6 else '')
                           for token in pieces)
            from reaction_batch.chemistry import formula
            counts = formula(text)
            if sum(counts.values()) > 5 or any(t not in {e['element'] for e in elements} for t in counts):
                continue
            base = sum(n * min(item['mass_number']
                               for item in next(e for e in elements if e['element'] == t)['isotopes'])
                       for t, n in counts.items())
            channels = []
            for c in range(random.randint(0, 3)):
                lo = base + random.randint(-1, 4)
                hi = lo + random.randint(0, 4)
                tlo = random.randint(0, 3)
                thi = tlo + random.randint(0, 3)
                if hi > base + 6:
                    hi = base + 6
                if lo > hi:
                    lo = hi
                channels.append(channel('c%d' % c, (max(0, lo), hi), (tlo, thi)))
            data = {'elements': elements,
                    'queries': [query('q', text, 6, channels)]}
            try:
                check_oracle(self, data)
            except ValueError:
                continue
            checked += 1
        self.assertGreater(checked, 20)


class Reordering(unittest.TestCase):
    def test_element_and_query_order_do_not_change_results(self):
        first = {'elements': [
            element('A', [iso(10, '0.5', False), iso(11, '0.5', True)]),
            element('B', [iso(5, '0.25', True), iso(6, '0.75', False)]),
        ], 'queries': [
            query('one', 'AB', 3, [channel('c', (15, 18), (0, 2))]),
            query('two', 'B2', 4, []),
        ]}
        second = {'elements': [
            element('B', [iso(6, '0.75', False), iso(5, '0.25', True)]),
            element('A', [iso(11, '0.5', True), iso(10, '0.5', False)]),
        ], 'queries': [
            query('one', 'BA', 3, [channel('c', (15, 18), (0, 2))]),
            query('two', 'B2', 4, []),
        ]}
        left = by_id(analyze(first))
        right = by_id(analyze(second))
        self.assertEqual(left['one'], right['one'])
        self.assertEqual(left['two'], right['two'])

    def test_zero_abundance_light_isotope_keeps_empty_envelope(self):
        elements = [element('A', [iso(10, '0', False), iso(11, '1', True)])]
        data = {'elements': elements,
                'queries': [query('q', 'A', 0, [])]}
        out = analyze(data)['queries'][0]
        self.assertEqual(out['base_mass'], 10)
        self.assertEqual(out['envelope'], [])
        self.assertEqual(out['overflow'],
                         {'probability': '1', 'first_moment': '1',
                          'second_factorial_moment': '0',
                          'mean_tags': '1', 'variance_tags': '0'})

    def test_zero_abundance_entry_sets_base_mass(self):
        elements = [element('A', [iso(10, '0', False), iso(12, '1', True)])]
        data = {'elements': elements,
                'queries': [query('q', 'A2', 4, [])]}
        out = analyze(data)['queries'][0]
        self.assertEqual(out['base_mass'], 20)
        self.assertEqual(out['envelope'], [{'mass': 24, 'probability': '1'}])
        self.assertEqual(out['total']['mean_tags'], '2')


class EdgeCases(unittest.TestCase):
    def test_tiny_but_nonzero_probability_stays_exact(self):
        data = {'elements': a2_table(),
                'queries': [query('q', 'A80', 80, [])]}
        out = analyze(data)['queries'][0]
        tail = out['envelope'][0]
        self.assertEqual(tail['mass'], 800)
        self.assertEqual(Fraction(tail['probability']), Fraction(1, 2 ** 80))
        self.assertNotEqual(Fraction(tail['probability']), 0)

    def test_sparse_mass_gaps_are_omitted_not_zero_filled(self):
        elements = [element('A', [iso(10, '0.5', False), iso(12, '0.5', False)])]
        data = {'elements': elements,
                'queries': [query('q', 'A2', 4, [])]}
        out = analyze(data)['queries'][0]
        self.assertEqual([entry['mass'] for entry in out['envelope']],
                         [20, 22, 24])
        self.assertEqual([entry['probability'] for entry in out['envelope']],
                         ['1/4', '1/2', '1/4'])

    def test_cutoff_zero_and_maximum(self):
        zero = analyze({'elements': a2_table(),
                        'queries': [query('q', 'A2', 0, [])]})['queries'][0]
        self.assertEqual(zero['envelope'],
                         [{'mass': 20, 'probability': '1/4'}])
        self.assertEqual(zero['overflow']['probability'], '3/4')
        top = analyze({'elements': a2_table(),
                       'queries': [query('q', 'A80', 128, [])]})['queries'][0]
        self.assertEqual(len(top['envelope']), 81)
        self.assertEqual(top['overflow']['probability'], '0')

    def test_empty_gate_beyond_atom_count(self):
        data = {'elements': a2_table(),
                'queries': [query('q', 'A2', 2,
                                  [channel('wide', (20, 22), (3, 80))])]}
        out = analyze(data)['queries'][0]
        self.assertEqual(out['channels'][0]['stats']['probability'], '0')
        self.assertEqual(out['missed'], out['total'])

    def test_no_channels_missed_equals_total(self):
        data = {'elements': a2_table(),
                'queries': [query('q', 'A2', 2, [])]}
        out = analyze(data)['queries'][0]
        self.assertEqual(out['channels'], [])
        self.assertEqual(out['missed'], out['total'])
        self.assertEqual(Fraction(sum(Fraction(e['probability'])
                                      for e in out['envelope']))
                         + Fraction(out['overflow']['probability']), 1)

    def test_every_query_is_independent(self):
        data = {'elements': a2_table(),
                'queries': [query('q1', 'A2', 2, [channel('c', (20, 22), (0, 2))]),
                            query('q2', 'A2', 2, [])]}
        first, second = analyze(data)['queries']
        self.assertEqual(first['channels'][0]['stats']['probability'], '1')
        self.assertEqual(second['channels'], [])
        self.assertEqual(second['missed'], second['total'])
        self.assertEqual(first['envelope'], second['envelope'])

    def test_unused_element_table_is_legal(self):
        data = {'elements': a2_table() + [element('Z', [iso(5, '1', False)])],
                'queries': [query('q', 'A', 0, [])]}
        out = analyze(data)['queries'][0]
        self.assertEqual(out['base_mass'], 10)


class CoverageRegressions(unittest.TestCase):
    def test_all_isotopes_share_one_tag_value(self):
        elements = [element('A', [iso(10, '0.25', True), iso(12, '0.75', True)])]
        data = {'elements': elements,
                'queries': [query('q', 'A2', 4, [])]}
        check_oracle(self, data)
        out = analyze(data)['queries'][0]
        self.assertEqual(out['total']['mean_tags'], '2')
        self.assertEqual(out['total']['variance_tags'], '0')
        self.assertEqual([e['mass'] for e in out['envelope']], [20, 22, 24])
        elements = [element('A', [iso(10, '0.25', False), iso(12, '0.75', False)])]
        data = {'elements': elements, 'queries': [query('q', 'A2', 4, [])]}
        out = analyze(data)['queries'][0]
        self.assertEqual(out['total']['mean_tags'], '0')
        self.assertEqual(out['total']['variance_tags'], '0')

    def test_nested_and_repeated_channel_bounds(self):
        channels = [channel('outer', (20, 22), (0, 2)),
                    channel('inner', (21, 21), (1, 1)),
                    channel('dup', (20, 22), (0, 2))]
        data = {'elements': a2_table(),
                'queries': [query('q', 'A2', 2, channels)]}
        check_oracle(self, data)
        out = analyze(data)['queries'][0]
        # 'outer' takes everything in priority; 'inner' and 'dup' are empty.
        self.assertEqual([c['stats']['probability'] for c in out['channels']],
                         ['1', '0', '0'])

    def test_same_formula_and_channel_ids_across_queries(self):
        data = {'elements': a2_table(), 'queries': [
            query('q1', 'A2', 2, [channel('c', (20, 20), (0, 0))]),
            query('q2', 'A2', 2, [channel('c', (22, 22), (2, 2))])]}
        out = analyze(data)['queries']
        self.assertEqual(out[0]['channels'][0]['stats']['probability'], '1/4')
        self.assertEqual(out[1]['channels'][0]['stats']['probability'], '1/4')

    def test_single_atom_and_cutoff_128(self):
        data = {'elements': a2_table(),
                'queries': [query('q', 'A', 128, [])]}
        out = analyze(data)['queries'][0]
        self.assertEqual(out['envelope'],
                         [{'mass': 10, 'probability': '1/2'},
                          {'mass': 11, 'probability': '1/2'}])

    def test_unreachable_channel_lower_above_max_mass(self):
        # Mass upper must be <= B+cutoff, but a lower bound can sit above every
        # reachable mass, legitimately giving an empty gate.
        data = {'elements': a2_table(),
                'queries': [query('q', 'A2', 2,
                                  [channel('high', (22, 22), (0, 0)),
                                   channel('all', (20, 22), (0, 2))])]}
        out = analyze(data)['queries'][0]
        self.assertEqual(out['channels'][0]['stats']['probability'], '0')
        self.assertEqual(out['channels'][1]['stats']['probability'], '1')


class InputValidation(unittest.TestCase):
    def base(self):
        return {'elements': a2_table(),
                'queries': [query('q', 'A2', 2, [channel('c', (20, 22), (0, 2))])]}

    def assert_rejected(self, mutate):
        data = self.base()
        mutate(data)
        with self.assertRaises(ValueError):
            analyze(data)

    def test_unknown_and_missing_fields(self):
        cases = [
            lambda d: d.__setitem__('extra', 1),
            lambda d: d['elements'][0].__setitem__('extra', 1),
            lambda d: d['elements'][0]['isotopes'][0].__setitem__('extra', 1),
            lambda d: d['queries'][0].__setitem__('extra', 1),
            lambda d: d['queries'][0]['channels'][0].__setitem__('extra', 1),
            lambda d: d['queries'][0]['channels'][0]['mass'].__setitem__('extra', 1),
            lambda d: d['queries'][0]['channels'][0]['tag_count'].__setitem__('extra', 1),
            lambda d: d['queries'][0].pop('cutoff'),
        ]
        for index, mutate in enumerate(cases):
            with self.subTest(case=index):
                self.assert_rejected(mutate)

    def test_literal_errors(self):
        cases = [
            lambda d: d['elements'][0].__setitem__('element', 'a'),
            lambda d: d['elements'][0].__setitem__('element', 'Abc'),
            lambda d: d['elements'][0]['isotopes'][0].__setitem__('mass_number', True),
            lambda d: d['elements'][0]['isotopes'][0].__setitem__('mass_number', 0),
            lambda d: d['elements'][0]['isotopes'][0].__setitem__('mass_number', 301),
            lambda d: d['elements'][0]['isotopes'][0].__setitem__('mass_number', 10.0),
            lambda d: d['elements'][0]['isotopes'][0].__setitem__('abundance', 0.5),
            lambda d: d['elements'][0]['isotopes'][0].__setitem__('abundance', '1.5'),
            lambda d: d['elements'][0]['isotopes'][0].__setitem__('abundance', '1.'),
            lambda d: d['elements'][0]['isotopes'][0].__setitem__('tagged', 1),
            lambda d: d['elements'][0]['isotopes'][1].__setitem__('mass_number', 10),
            lambda d: d['elements'][0]['isotopes'][0].__setitem__('abundance', '0.4'),
            lambda d: d['queries'][0].__setitem__('cutoff', True),
            lambda d: d['queries'][0].__setitem__('cutoff', 129),
            lambda d: d['queries'][0]['channels'][0]['mass'].__setitem__('lower', True),
            lambda d: d['queries'][0]['channels'][0]['mass'].__setitem__('lower', 23),
            lambda d: d['queries'][0]['channels'][0]['mass'].__setitem__('upper', 10 ** 9),
            lambda d: d['queries'][0]['channels'][0]['mass'].__setitem__('upper', 23),
            lambda d: d['queries'][0]['channels'][0]['tag_count'].__setitem__('lower', 81),
            lambda d: d['queries'][0]['channels'][0]['tag_count'].__setitem__('lower', 3),
            lambda d: d['queries'][0]['channels'][0]['tag_count'].__setitem__('upper', True),
        ]
        for index, mutate in enumerate(cases):
            with self.subTest(case=index):
                self.assert_rejected(mutate)

    def test_structure_limits(self):
        cases = [
            lambda d: d.__setitem__('elements', []),
            lambda d: d.__setitem__('queries', []),
            lambda d: d['elements'][0].__setitem__('isotopes', []),
            lambda d: d['queries'][0].__setitem__('channels',
                                                  [channel('c%d' % i, (20, 22), (0, 2))
                                                   for i in range(17)]),
            lambda d: d['queries'][0].__setitem__('formula', 'A81'),
            lambda d: d['queries'][0].__setitem__('formula', 'AB'),
            lambda d: d['queries'][0].__setitem__('formula', 'A0'),
            lambda d: d['queries'].append(query('q', 'A', 0, [])),
            lambda d: d['queries'][0]['channels'].append(
                channel('c', (20, 22), (0, 2))),
            lambda d: d['elements'].append(element('A', [iso(1, '1', False)])),
            lambda d: d['elements'].append(element('B', [iso(1, '1', False),
                                                          iso(2, '1', False)])),
            lambda d: d['queries'].extend(
                query('q%d' % i, 'A', 0, []) for i in range(2, 6)),
        ]
        for index, mutate in enumerate(cases):
            with self.subTest(case=index):
                self.assert_rejected(mutate)


class CliRegression(unittest.TestCase):
    def run_cli(self, source, output, module='reaction_batch.isotopes_cli'):
        return subprocess.run(
            [sys.executable, '-m', module,
             '--source', str(source), '--output', str(output)],
            capture_output=True, cwd=ROOT,
            env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')})

    def write(self, tmp, payload):
        source = Path(tmp) / 'in.json'
        source.write_text(json.dumps(payload))
        return source

    def test_cli_round_trip_and_digest(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = self.write(tmp, {'elements': a2_table(),
                                      'queries': [query('q', 'A2', 2, [])]})
            output = Path(tmp) / 'out.json'
            completed = self.run_cli(source, output)
            self.assertEqual(completed.returncode, 0, completed.stderr[:400])
            report = json.loads(output.read_text())
            self.assertEqual(report['schema_version'], 1)
            self.assertEqual(report['source_sha256'],
                             hashlib.sha256(source.read_bytes()).hexdigest())
            self.assertEqual(set(report), {'schema_version', 'source_sha256',
                                           'result'})

    def test_rerun_byte_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = self.write(tmp, {'elements': a2_table(),
                                      'queries': [query('q', 'A80', 80, [])]})
            first, second = Path(tmp) / 'a.json', Path(tmp) / 'b.json'
            self.assertEqual(self.run_cli(source, first).returncode, 0)
            self.assertEqual(self.run_cli(source, second).returncode, 0)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_late_query_error_preserves_old_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = {'elements': a2_table(),
                       'queries': [query('ok', 'A', 1, []),
                                   query('bad', 'Z', 1, [])]}
            source = self.write(tmp, payload)
            output = Path(tmp) / 'out.json'
            output.write_text('old report')
            completed = self.run_cli(source, output)
            self.assertEqual(completed.returncode, 2)
            self.assertEqual(output.read_text(), 'old report')

    def test_output_guards(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = self.write(tmp, {'elements': a2_table(),
                                      'queries': [query('q', 'A', 1, [])]})
            link = Path(tmp) / 'link.json'
            target = Path(tmp) / 'target.json'
            target.write_text('keep')
            link.symlink_to(target)
            for destination in (source, link):
                completed = self.run_cli(source, destination)
                self.assertEqual(completed.returncode, 2)
            self.assertEqual(target.read_text(), 'keep')

    def test_duplicate_json_keys_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'dup.json'
            source.write_text('{"elements": [], "elements": [], "queries": []}')
            output = Path(tmp) / 'out.json'
            completed = self.run_cli(source, output)
            self.assertEqual(completed.returncode, 2)
            self.assertFalse(output.exists())

    def test_existing_commands_still_work(self):
        cases = [('reaction_batch.cli', 'input.json'),
                 ('reaction_batch.lots', 'lots-example.json'),
                 ('reaction_batch.assay_cli', 'assay-example1.json'),
                 ('reaction_batch.peaks_cli', 'peaks-example1.json')]
        with tempfile.TemporaryDirectory() as tmp:
            for index, (module, fixture) in enumerate(cases):
                with self.subTest(module=module):
                    source = ROOT / 'fixtures' / fixture
                    output = Path(tmp) / ('out%d.json' % index)
                    completed = self.run_cli(source, output, module=module)
                    self.assertEqual(completed.returncode, 0,
                                     completed.stderr[:400])
                    self.assertEqual(json.loads(output.read_text())['schema_version'], 1)


if __name__ == '__main__':
    unittest.main()
