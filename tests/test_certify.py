"""Tests for reaction-certify: documented acceptance, independent oracle
cross-checks with change-permutation sweeps, input-contract errors, CLI
regression and the full 300-batch / 6000-change scale boundary.
"""

import itertools
import json
import os
import random
import subprocess
import sys
import tempfile
import time
import unittest
from fractions import Fraction
from pathlib import Path

from reaction_batch.certify import analyze
from reaction_batch.cli import run_certify
from reference_certify import replay

ROOT = Path(__file__).resolve().parents[1]

# --- record builders -------------------------------------------------------


def cal(cal_id, instrument='i', start=0, end=20, points=None):
    if points is None:
        points = [('0', '0'), ('10', '10')]
    return {'id': cal_id, 'instrument': instrument, 'start': start, 'end': end,
            'points': [{'raw': raw, 'value': value} for raw, value in points]}


def obs(obs_id, instrument='i', time=10, formula_text='H2O', raw='4'):
    return {'id': obs_id, 'instrument': instrument, 'time': time,
            'formula': formula_text, 'raw': raw}


def cert(cert_id, formula_text='H2O', members=None, lower='0', upper='100',
         max_range='100'):
    if members is None:
        members = [('o', '1')]
    return {'id': cert_id, 'formula': formula_text,
            'members': [{'observation_id': m, 'weight': w} for m, w in members],
            'lower': lower, 'upper': upper, 'max_range': max_range}


def put(entity, record):
    return {'entity': entity, 'action': 'put', 'record': record}


def delete(entity, record_id):
    return {'entity': entity, 'action': 'delete', 'id': record_id}


def batch(batch_id, *changes):
    return {'batch_id': batch_id, 'changes': list(changes)}


def query(query_id, as_of):
    return {'query_id': query_id, 'as_of': as_of}


def data(batches, queries=None):
    return {'batches': batches,
            'queries': queries if queries is not None else [query('q', 1)]}


def statuses(result):
    return [(item['status'], item['reason']) for item in result['batches']]


# --- documented acceptance walkthrough -------------------------------------


class DocumentedAcceptance(unittest.TestCase):
    def test_identity_then_double_then_double_with_withdrawal(self):
        source = json.loads((ROOT / 'fixtures/certify-example.json').read_text())
        result = analyze(source)
        self.assertEqual(statuses(result),
                         [('accepted', None),
                          ('rejected', 'mean_out_of_range'),
                          ('accepted', None)])
        by_query = {item['query_id']: item for item in result['queries']}
        q1, q2, q3, q0 = (by_query[name] for name in ('q1', 'q2', 'q3', 'q0'))
        for view in (q1, q2):
            self.assertEqual(len(view['observations']), 1)
            observation = view['observations'][0]
            self.assertEqual(observation['value'], '4')
            self.assertEqual(observation['calibration_id'], 'c')
            self.assertEqual(view['certificates'][0]['members'][0]['value'], '4')
            self.assertEqual(view['certificates'][0]['mean'], '4')
        self.assertEqual(q3['observations'][0]['value'], '8')
        self.assertEqual(q3['certificates'], [])
        self.assertEqual(q0['observations'], [])
        self.assertEqual(q0['certificates'], [])

    def test_fixture_matches_cli_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'in.json'
            output = Path(tmp) / 'out.json'
            source.write_text((ROOT / 'fixtures/certify-example.json').read_text())
            completed = subprocess.run(
                [sys.executable, '-m', 'reaction_batch.certify_cli',
                 '--source', str(source), '--output', str(output)],
                capture_output=True, cwd=ROOT,
                env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')})
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = json.loads(output.read_text())
            self.assertEqual(report['schema_version'], 1)
            self.assertEqual(len(report['source_sha256']), 64)
            self.assertEqual(len(report['result']['batches']), 3)


# --- supplemental acceptance cases -----------------------------------------


class HalfOpenBoundaries(unittest.TestCase):
    def test_time_equal_to_start_chooses_later_interval(self):
        source = data([
            batch('b1', put('calibration', cal('c1', start=0, end=10)),
                  put('calibration', cal('c2', start=10, end=20)),
                  put('observation', obs('o', time=10))),
            batch('b2', delete('calibration', 'c2')),
        ], [query('q1', 1), query('q2', 2)])
        result = analyze(source)
        self.assertEqual(statuses(result), [('accepted', None),
                                            ('accepted', None)])
        views = {item['as_of']: item for item in result['queries']}
        self.assertEqual(views[1]['observations'][0]['calibration_id'], 'c2')
        gone = views[2]['observations'][0]
        self.assertEqual(gone['status'], 'unavailable')
        self.assertEqual(gone['reason'], 'no_calibration')
        self.assertIsNone(gone['calibration_id'])

    def test_removing_covering_calibration_rejected_while_certificate_active(self):
        source = data([
            batch('b1', put('calibration', cal('c1', start=0, end=10)),
                  put('calibration', cal('c2', start=10, end=20)),
                  put('observation', obs('o', time=10)),
                  put('certificate', cert('s', lower='4', upper='4',
                                          max_range='0'))),
            batch('b2', delete('calibration', 'c2')),
            batch('b3', delete('calibration', 'c2'),
                  delete('certificate', 's')),
        ], [query('q3', 3)])
        result = analyze(source)
        self.assertEqual(statuses(result),
                         [('accepted', None),
                          ('rejected', 'unavailable_member'),
                          ('accepted', None)])
        self.assertEqual(result['queries'][0]['certificates'], [])

    def test_end_boundary_excludes_first_interval(self):
        source = data([
            batch('b1', put('calibration', cal('c1', start=0, end=10)),
                  put('observation', obs('o', time=9, raw='10'))),
        ], [query('q', 1)])
        item = analyze(source)['queries'][0]['observations'][0]
        self.assertEqual(item['status'], 'available')
        self.assertEqual(item['value'], '10')
        outside = analyze(data([
            batch('b1', put('calibration', cal('c1', start=0, end=10)),
                  put('observation', obs('o', time=10))),
        ], [query('q', 1)]))['queries'][0]['observations'][0]
        self.assertEqual(outside['reason'], 'no_calibration')


class WeightedStatistics(unittest.TestCase):
    def test_weights_one_and_three_mean_and_range_priority(self):
        members = [('o0', '1'), ('o1', '3')]
        curve = put('calibration', cal('c', points=[('0', '0'), ('10', '10')]))
        observations = [
            put('observation', obs('o0', time=1, raw='0')),
            put('observation', obs('o1', time=2, raw='10'))]
        # mean = 15/2 sits inside [7,8]; only range = 10 > 9 fails.
        strict_range = analyze(data([
            batch('b1', curve, *observations,
                  put('certificate', cert('s', members=members, lower='7',
                                          upper='8', max_range='9'))),
        ], [query('q', 1)]))
        self.assertEqual(statuses(strict_range),
                         [('rejected', 'range_exceeded')])
        # Widening the range to exactly 10 lets the same certificate issue.
        accepted = analyze(data([
            batch('b1', curve, *observations,
                  put('certificate', cert('s', members=members, lower='7',
                                          upper='8', max_range='10'))),
        ], [query('q', 1)]))['queries'][0]['certificates'][0]
        self.assertEqual(accepted['mean'], '15/2')
        self.assertEqual(accepted['range'], '10')
        self.assertEqual([m['observation_id'] for m in accepted['members']],
                         ['o0', 'o1'])
        # A mean envelope excluding 15/2 with a loose range fails on mean.
        mean_bad = analyze(data([
            batch('b1', curve, *observations,
                  put('certificate', cert('s', members=members, lower='0',
                                          upper='7', max_range='10'))),
        ], [query('q', 1)]))
        self.assertEqual(statuses(mean_bad),
                         [('rejected', 'mean_out_of_range')])
        # Replacing a failing certificate in place stays rejected; the prior
        # accepted view is what the prefix query reports.
        replacement = analyze(data([
            batch('b1', curve, *observations,
                  put('certificate', cert('s', members=members, lower='7',
                                          upper='8', max_range='10'))),
            batch('b2', put('certificate', cert('s', members=members,
                                                lower='7', upper='8',
                                                max_range='9'))),
        ], [query('q1', 1), query('q2', 2)]))
        self.assertEqual(statuses(replacement),
                         [('accepted', None), ('rejected', 'range_exceeded')])
        self.assertEqual(
            replacement['queries'][1]['certificates'][0]['max_range'], '10')


class FormulaEquivalence(unittest.TestCase):
    def test_h2o_and_oh2_are_the_same_formula(self):
        result = analyze(data([
            batch('b1', put('calibration', cal('c')),
                  put('observation', obs('o', formula_text='H2O')),
                  put('certificate', cert('s', formula_text='OH2',
                                          lower='4', upper='4',
                                          max_range='0'))),
        ]))
        self.assertEqual(statuses(result), [('accepted', None)])

    def test_a_and_a2_never_reduce_to_each_other(self):
        result = analyze(data([
            batch('b1', put('calibration', cal('c')),
                  put('observation', obs('o', formula_text='A')),
                  put('certificate', cert('s', formula_text='A2',
                                          lower='4', upper='4',
                                          max_range='0'))),
        ]))
        self.assertEqual(statuses(result),
                         [('rejected', 'formula_mismatch')])


class CurveSemantics(unittest.TestCase):
    def test_non_monotonic_curve_and_internal_node_equality(self):
        source = data([
            batch('b1', put('calibration', cal(
                'c', points=[('0', '10'), ('5', '2'), ('10', '8')])),
                  put('observation', obs('at_node', time=1, raw='5')),
                  put('observation', obs('mid', time=2, raw='2.5')),
                  put('observation', obs('beyond', time=3, raw='11')),
                  put('observation', obs('before', time=4, raw='0.5'))),
        ], [query('q', 1)])
        items = {item['id']: item
                 for item in analyze(source)['queries'][0]['observations']}
        self.assertEqual(items['at_node']['value'], '2')
        # y0 + (2.5-0)*(2-10)/(5-0) = 10 - 4 = 6
        self.assertEqual(items['mid']['value'], '6')
        self.assertEqual(items['beyond']['reason'], 'outside_curve')
        self.assertEqual(items['beyond']['calibration_id'], 'c')
        self.assertIsNone(items['beyond']['value'])
        # 10 + 0.5 * (-8)/5 = 10 - 4/5 = 46/5
        self.assertEqual(items['before']['value'], '46/5')

    def test_exact_one_third_interpolation(self):
        result = analyze(data([
            batch('b1', put('calibration', cal(
                'c', points=[('0', '0'), ('3', '1')])),
                  put('observation', obs('o', raw='1'))),
        ], [query('q', 1)]))['queries'][0]['observations'][0]
        self.assertEqual(result['value'], '1/3')

    def test_endpoint_nodes_are_inside_domain(self):
        result = analyze(data([
            batch('b1', put('calibration', cal(
                'c', points=[('2', '7'), ('9', '21')])),
                  put('observation', obs('lo', raw='2')),
                  put('observation', obs('hi', raw='9'))),
        ], [query('q', 1)]))['queries'][0]['observations']
        self.assertEqual([item['value'] for item in result], ['21', '7'])


class ReplayAndReuse(unittest.TestCase):
    def test_replacement_then_correction_keeps_calibration_identity(self):
        # One change per (entity, id) per batch; replacement spans batches.
        source = data([
            batch('b1', put('calibration', cal('c')),
                  put('observation', obs('o')),
                  put('certificate', cert('s', lower='4', upper='4',
                                          max_range='0'))),
            batch('b2', put('calibration', cal(
                'c', points=[('0', '0'), ('10', '20')])),
                  put('certificate', cert('s', lower='8', upper='8',
                                          max_range='0'))),
            batch('b3', delete('calibration', 'c')),
        ], [query('q1', 1), query('q2', 2), query('q3', 3)])
        result = analyze(source)
        self.assertEqual(statuses(result),
                         [('accepted', None), ('accepted', None),
                          ('rejected', 'unavailable_member')])
        views = {item['as_of']: item for item in result['queries']}
        self.assertEqual(views[1]['observations'][0]['value'], '4')
        self.assertEqual(views[2]['observations'][0]['value'], '8')
        # The rejected deletion leaves the doubled curve in place.
        self.assertEqual(views[3]['observations'][0]['value'], '8')

    def test_split_interval_and_fill_gap_in_one_batch(self):
        source = data([
            batch('b1', put('calibration', cal(
                'wide', start=0, end=40,
                points=[('0', '0'), ('10', '10')])),
                  put('observation', obs('early', time=5, raw='4')),
                  put('observation', obs('late', time=35, raw='4'))),
            batch('b2',
                  delete('calibration', 'wide'),
                  put('calibration', cal(
                      'left', start=0, end=20,
                      points=[('0', '0'), ('10', '10')])),
                  put('calibration', cal(
                      'right', start=20, end=40,
                      points=[('0', '0'), ('10', '10')]))),
        ], [query('q2', 2)])
        result = analyze(source)
        self.assertEqual(statuses(result),
                         [('accepted', None), ('accepted', None)])
        items = {item['id']: item
                 for item in result['queries'][0]['observations']}
        self.assertEqual(items['early']['calibration_id'], 'left')
        self.assertEqual(items['late']['calibration_id'], 'right')
        # Filling a previously empty time window in the same batch is legal.
        filled = analyze(data([
            batch('b1', put('calibration', cal('only', start=0, end=10)),
                  put('observation', obs('o', time=30))),
            batch('b2', put('calibration', cal('fill', start=20, end=40))),
        ], [query('q', 2)]))['queries'][0]['observations'][0]
        self.assertEqual(filled['status'], 'available')
        self.assertEqual(filled['calibration_id'], 'fill')
    def test_same_batch_member_and_calibration_replacement_is_order_independent(self):
        base = [
            batch('b1', put('calibration', cal('c')),
                  put('observation', obs('o')),
                  put('certificate', cert('s', lower='4', upper='4',
                                          max_range='0'))),
        ]
        b2 = batch('b2',
                   put('calibration', cal(
                       'c', points=[('0', '0'), ('10', '20')])),
                   put('certificate', cert('s', lower='8', upper='8',
                                           max_range='0')))
        queries = [query('q1', 1), query('q2', 2)]
        first = analyze(data(base + [b2], queries))
        shuffled = dict(b2)
        shuffled['changes'] = list(reversed(b2['changes']))
        second = analyze(data(base + [shuffled], queries))
        self.assertEqual(first, second)
        self.assertEqual(statuses(second), [('accepted', None),
                                             ('accepted', None)])
        self.assertEqual(second['queries'][1]['certificates'][0]['mean'], '8')

    def test_deleted_id_rebuilt_keeps_old_query_isolated(self):
        source = data([
            batch('b1', put('observation', obs('o', raw='4'))),
            batch('b2', delete('observation', 'o')),
            batch('b3', put('observation', obs('o', raw='9'))),
        ], [query('q1', 1), query('q2', 2), query('q3', 3)])
        views = {item['as_of']: item for item in analyze(source)['queries']}
        self.assertEqual(views[1]['observations'][0]['raw'], '4')
        self.assertEqual(views[2]['observations'], [])
        self.assertEqual(views[3]['observations'][0]['raw'], '9')

    def test_rejected_batch_publishes_no_partial_change(self):
        source = data([
            batch('b1', put('calibration', cal('c')),
                  put('observation', obs('o')),
                  put('certificate', cert('s', lower='4', upper='4',
                                          max_range='0'))),
            batch('b2', put('observation', obs('o2', formula_text='A')),
                  put('certificate', cert('s2', formula_text='A2',
                                          lower='0', upper='0',
                                          max_range='0'))),
        ], [query('q', 2)])
        result = analyze(source)
        self.assertEqual(statuses(result)[1],
                         ('rejected', 'formula_mismatch'))
        self.assertEqual(
            [item['id'] for item in result['queries'][0]['observations']],
            ['o'])


class ReasonCoverage(unittest.TestCase):
    def _one_batch(self, changes):
        return analyze(data([batch('b1', *changes)]))['batches'][0]

    def test_each_business_reason(self):
        self.assertEqual(self._one_batch([delete('observation', 'ghost')]),
                         {'batch_id': 'b1', 'status': 'rejected',
                          'reason': 'missing_delete'})
        overlap = self._one_batch([
            put('calibration', cal('c1', start=0, end=12)),
            put('calibration', cal('c2', start=10, end=20))])
        self.assertEqual(overlap['reason'], 'overlap')
        missing_member = self._one_batch([
            put('certificate', cert('s', members=[('ghost', '1')]))])
        self.assertEqual(missing_member['reason'], 'missing_member')
        unavailable = self._one_batch([
            put('observation', obs('o', time=99)),
            put('certificate', cert('s'))])
        self.assertEqual(unavailable['reason'], 'unavailable_member')
        mean_bad = self._one_batch([
            put('calibration', cal('c')), put('observation', obs('o')),
            put('certificate', cert('s', lower='5', upper='6',
                                    max_range='100'))])
        self.assertEqual(mean_bad['reason'], 'mean_out_of_range')
        range_bad = self._one_batch([
            put('calibration', cal('c', points=[('0', '0'), ('10', '10')])),
            put('observation', obs('a', time=1, raw='0')),
            put('observation', obs('b', time=2, raw='10')),
            put('certificate', cert('s', members=[('a', '1'), ('b', '1')],
                                    lower='0', upper='100', max_range='9'))])
        self.assertEqual(range_bad['reason'], 'range_exceeded')

    def test_reason_priority_fixed_order(self):
        # One candidate carries several failures; the highest priority wins.
        result = self._one_batch([
            delete('observation', 'ghost'),
            put('calibration', cal('c1', start=0, end=12)),
            put('calibration', cal('c2', start=10, end=20)),
            put('certificate', cert('s', members=[('ghost', '1')],
                                    formula_text='A2', lower='9', upper='9',
                                    max_range='0')),
        ])
        self.assertEqual(result['reason'], 'missing_delete')
        result = self._one_batch([
            put('calibration', cal('c1', start=0, end=12)),
            put('calibration', cal('c2', start=10, end=20)),
            put('observation', obs('o', formula_text='A')),
            put('certificate', cert('s', formula_text='A2',
                                    lower='9', upper='9', max_range='0')),
        ])
        self.assertEqual(result['reason'], 'overlap')
        result = self._one_batch([
            put('observation', obs('o', formula_text='A')),
            put('certificate', cert('s', formula_text='A2',
                                    members=[('o', '1'), ('ghost', '1')],
                                    lower='9', upper='9', max_range='0')),
        ])
        self.assertEqual(result['reason'], 'missing_member')
        result = self._one_batch([
            put('observation', obs('o', formula_text='A', time=99)),
            put('certificate', cert('s', formula_text='A2', lower='9',
                                    upper='9', max_range='0')),
        ])
        self.assertEqual(result['reason'], 'formula_mismatch')
        result = self._one_batch([
            put('observation', obs('o', time=99)),
            put('certificate', cert('s', lower='9', upper='9',
                                    max_range='0')),
        ])
        self.assertEqual(result['reason'], 'unavailable_member')


class EmptyInput(unittest.TestCase):
    def test_empty_state(self):
        result = analyze({'batches': [], 'queries': [query('q0', 0)]})
        self.assertEqual(result['batches'], [])
        view = result['queries'][0]
        self.assertEqual(view['observations'], [])
        self.assertEqual(view['certificates'], [])
        self.assertEqual(view['as_of'], 0)


# --- independent oracle cross-check ----------------------------------------


FORMULAS = ['H2O', 'OH2', 'A', 'A2', 'CO2', 'C(O2)']


def _random_scenario(seed, batch_count):
    rng = random.Random(seed)
    sequences = {kind: [] for kind in ('calibration', 'observation',
                                       'certificate')}
    live = {kind: set() for kind in sequences}
    counter = 0

    def new_id(kind):
        nonlocal counter
        counter += 1
        return kind[0] + str(counter)

    batches = []
    for index in range(batch_count):
        changes = []
        touched = set()
        for _ in range(rng.randint(1, 4)):
            kind = rng.choice(('calibration', 'observation', 'certificate'))
            if live[kind] and rng.random() < 0.3:
                target = rng.choice(sorted(live[kind]))
                key = (kind, target)
                if key in touched:
                    continue
                touched.add(key)
                changes.append(delete(kind, target))
                live[kind].discard(target)
                continue
            record_id = new_id(kind)
            key = (kind, record_id)
            if key in touched:
                continue
            touched.add(key)
            live[kind].add(record_id)
            if kind == 'calibration':
                start = rng.randint(0, 40)
                end = start + rng.randint(1, 40)
                raws = sorted(rng.sample(range(0, 30), 3))
                points = [(str(raw), str(rng.randint(0, 25))) for raw in raws]
                record = cal(record_id, instrument=rng.choice(['i', 'j']),
                             start=start, end=end, points=points)
            elif kind == 'observation':
                record = obs(record_id,
                             instrument=rng.choice(['i', 'j']),
                             time=rng.randint(0, 80),
                             formula_text=rng.choice(FORMULAS),
                             raw=str(rng.randint(0, 30)))
            else:
                pool = sorted(live['observation'])
                if not pool:
                    live[kind].discard(record_id)
                    touched.discard(key)
                    continue
                chosen = rng.sample(pool, min(len(pool), rng.randint(1, 3)))
                record = cert(
                    record_id, formula_text=rng.choice(FORMULAS),
                    members=[(m, str(rng.randint(1, 5))) for m in chosen],
                    lower=str(rng.randint(0, 10)),
                    upper=str(rng.randint(10, 30)),
                    max_range=str(rng.randint(0, 30)))
            changes.append(put(kind, record))
        batches.append(batch('b' + str(index), *changes))
    queries = [query('q' + str(k), rng.randint(0, batch_count))
               for k in range(6)]
    return {'batches': batches, 'queries': queries}


class OracleCrossCheck(unittest.TestCase):
    def test_random_sequences_match_reference(self):
        for seed in range(40):
            source = _random_scenario(seed, (seed % 7) + 2)
            expected_statuses, expected_queries = replay(source)
            result = analyze(source)
            self.assertEqual(result['batches'], expected_statuses)
            self.assertEqual(result['queries'], expected_queries)

    def test_reference_covers_each_rejection_reason(self):
        # Hand-shaped interlocking sequences so the independent oracle
        # witnesses every business reason, including mean and range failures.
        curve = put('calibration', cal('c', points=[('0', '0'), ('10', '10')]))
        crafted = [
            # missing_delete on the very first batch.
            batch('d0', delete('observation', 'ghost')),
            # Overlap from two covering intervals.
            batch('d1', put('calibration', cal('x1', start=0, end=12)),
                  put('calibration', cal('x2', start=10, end=20))),
            # A healthy accepted base for the remaining reasons.
            batch('ok', curve,
                  put('observation', obs('a', time=1, raw='0')),
                  put('observation', obs('b', time=2, raw='10')),
                  put('observation', obs('fm', formula_text='A', time=1,
                                         raw='4')),
                  put('certificate', cert('good', members=[('a', '1'),
                                                            ('b', '1')],
                                          lower='0', upper='100',
                                          max_range='100'))),
            batch('missing', put('certificate', cert('m',
                                                     members=[('nope', '1')]))),
            batch('formula', put('certificate', cert(
                'f', formula_text='A2', members=[('fm', '1')]))),
            batch('unavailable', put('observation', obs('u', time=99)),
                  put('certificate', cert('u_c', members=[('u', '1')]))),
            batch('mean', put('certificate', cert(
                'mn', members=[('a', '1')], lower='9', upper='10',
                max_range='100'))),
            batch('range', put('certificate', cert(
                'rg', members=[('a', '1'), ('b', '1')], lower='0',
                upper='100', max_range='9'))),
        ]
        source = data(crafted, [query('q', len(crafted))])
        expected_statuses, expected_queries = replay(source)
        result = analyze(source)
        self.assertEqual(result['batches'], expected_statuses)
        self.assertEqual(result['queries'], expected_queries)
        observed = {item['reason'] for item in result['batches']
                    if item['reason']}
        self.assertEqual(observed, {
            'missing_delete', 'overlap', 'missing_member',
            'formula_mismatch', 'unavailable_member',
            'mean_out_of_range', 'range_exceeded'})

    def test_change_permutation_invariance_against_reference(self):
        for seed in range(12):
            source = _random_scenario(100 + seed, 5)
            multi = [batch for batch in source['batches']
                     if len(batch['changes']) > 1]
            if not multi:
                continue
            target = multi[0]
            permutations = list(itertools.permutations(target['changes']))
            baseline = analyze(source)
            for ordering in permutations:
                shuffled = json.loads(json.dumps(source))
                for candidate in shuffled['batches']:
                    if candidate['batch_id'] == target['batch_id']:
                        candidate['changes'] = [
                            json.loads(json.dumps(change))
                            for change in ordering]
                expected_statuses, expected_queries = replay(shuffled)
                result = analyze(shuffled)
                self.assertEqual(result['batches'], expected_statuses)
                self.assertEqual(result['queries'], expected_queries)
                self.assertEqual(result['batches'], baseline['batches'])

    def test_reference_tracks_rejected_prefixes(self):
        source = _random_scenario(7, 6)
        expected_statuses, _ = replay(source)
        result = analyze(source)
        rejected = [index for index, item in enumerate(expected_statuses)
                    if item['status'] == 'rejected']
        self.assertTrue(rejected)
        # Every prefix, including those ending in a rejected batch, must agree
        # with the reference view rendered at that prefix.
        for prefix in range(len(source['batches']) + 1):
            probe = json.loads(json.dumps(source))
            probe['queries'] = [query('probe', prefix)]
            _, expected = replay(probe)
            actual = analyze(probe)['queries']
            self.assertEqual(actual, expected)


# --- whole-file input errors ------------------------------------------------


class InputContract(unittest.TestCase):
    def assert_invalid(self, source):
        with self.assertRaises((ValueError, TypeError)):
            analyze(source)

    def test_unknown_and_missing_fields(self):
        self.assert_invalid({'batches': [], 'queries': [], 'extra': 1})
        self.assert_invalid({'batches': []})
        self.assert_invalid({'batches': {}, 'queries': []})
        self.assert_invalid(data([{'batch_id': 'b', 'changes': [],
                                   'extra': 1}]))
        self.assert_invalid(data([batch('b', {'entity': 'calibration',
                                              'action': 'put'})]))
        self.assert_invalid(data([batch('b', {'entity': 'ghost',
                                              'action': 'delete',
                                              'id': 'x'})]))
        self.assert_invalid(data([batch('b', {'entity': 'observation',
                                              'action': 'freeze',
                                              'id': 'x'})]))

    def test_literal_and_bound_errors(self):
        self.assert_invalid(data([batch('b', put('observation',
                                                 obs('o', raw='1/3')))]))
        self.assert_invalid(data([batch('b', put('observation',
                                                 obs('o', raw=4)))]))
        self.assert_invalid(data([batch('b', put('observation',
                                                 obs('o', time=True)))]))
        self.assert_invalid(data([batch('b', put('calibration',
                                                 cal('c', start=10, end=10)))]))
        self.assert_invalid(data([batch('b', put('calibration', cal(
            'c', points=[('5', '0'), ('5', '1')])))]))
        self.assert_invalid(data([batch('b', put('calibration', cal(
            'c', points=[('5', '0')])))]))
        self.assert_invalid(data([batch('b', put('certificate', cert(
            's', members=[('o', '0')])))]))
        self.assert_invalid(data([batch('b', put('certificate', cert(
            's', members=[('o', '-1')])))]))
        self.assert_invalid(data([batch('b', put('certificate', cert(
            's', lower='9', upper='8')))]))
        self.assert_invalid(data([batch('b', put('certificate', cert(
            's', members=[('o', '1'), ('o', '2')])))]))
        self.assert_invalid(data([batch('b', put('observation',
                                                 obs('o', formula_text='H0')))]))
        self.assert_invalid(data([batch('b', put('observation',
                                                 obs('', time=1)))]))

    def test_duplicate_and_as_of_errors(self):
        self.assert_invalid(data([batch('b'), batch('b')]))
        self.assert_invalid({'batches': [batch('b')],
                             'queries': [query('q', 2)]})
        self.assert_invalid({'batches': [], 'queries': [query('q', -1)]})
        self.assert_invalid({'batches': [], 'queries': [query('q', 0),
                                                         query('q', 0)]})
        self.assert_invalid(data([batch('b',
                                        put('observation', obs('o')),
                                        put('observation', obs('o')))]))

    def test_trailing_record_still_validated(self):
        valid = data([batch('good', put('observation', obs('o'))),
                      batch('bad', put('observation', obs('o2', raw='nope')))])
        self.assert_invalid(valid)

    def test_duplicate_json_key_rejected_at_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'dup.json'
            path.write_text('{"batches":[],"batches":[],"queries":[]}')
            with self.assertRaises(ValueError):
                run_certify(path)


# --- CLI transport ----------------------------------------------------------


class CliRegression(unittest.TestCase):
    def test_bad_input_exits_two_and_keeps_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'in.json'
            output = root / 'out.json'
            source.write_text('{}')
            output.write_text('previous')
            completed = subprocess.run(
                [sys.executable, '-m', 'reaction_batch.certify_cli',
                 '--source', str(source), '--output', str(output)],
                capture_output=True, cwd=ROOT,
                env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')})
            self.assertEqual(completed.returncode, 2)
            self.assertEqual(output.read_text(), 'previous')

    def test_same_path_and_symlink_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'in.json'
            source.write_text(json.dumps(data([])))
            same = subprocess.run(
                [sys.executable, '-m', 'reaction_batch.certify_cli',
                 '--source', str(source), '--output', str(source)],
                capture_output=True, cwd=ROOT,
                env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')})
            self.assertEqual(same.returncode, 2)
            link = root / 'link'
            link.symlink_to(source)
            linked = subprocess.run(
                [sys.executable, '-m', 'reaction_batch.certify_cli',
                 '--source', str(source), '--output', str(link)],
                capture_output=True, cwd=ROOT,
                env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')})
            self.assertEqual(linked.returncode, 2)


# --- full-scale boundary ----------------------------------------------------


def _scale_input():
    """Hit the batch/change/member/query scale boundaries in one legal input.

    300 batches with exactly 6000 cumulative changes (20 per batch), 2000
    cumulative certificate members (50 certificates of 40 members each) and
    60 shuffled historical queries.  Every batch puts 5 single-instrument
    calibrations whose validity windows tile the time axis without overlap
    (1500 calibrations, 4500 nodes).  Batches 0..99 add 15 observations each
    (1500 observations); batches 100..149 add one 40-member certificate over
    observations already present plus 14 observations; later batches add 15
    observations each.  Every curve covers raw 0..50, so all observations
    stay available and every certificate remains valid; all 300 batches
    accept.
    """
    batches = []
    changes = 0
    nodes = 0
    members = 0
    observation_index = 0
    for index in range(300):
        current = []
        for slot in range(5):
            points = [('0', '0'), ('25', '25'), ('50', '50')]
            current.append(put('calibration', cal(
                f'c{index}_{slot}', instrument=f'i{slot}',
                start=index * 1000, end=index * 1000 + 1000,
                points=points)))
            nodes += len(points)
        changes += 5
        if 100 <= index < 150:
            start_member = (index - 100) * 40
            chosen = [f'o{k}' for k in range(start_member, start_member + 40)]
            current.append(put('certificate', cert(
                f's{index}', members=[(name, '1') for name in chosen],
                lower='0', upper='19', max_range='19')))
            changes += 1
            members += 40
            extra = 14
        else:
            extra = 15
        for slot in range(extra):
            current.append(put(
                'observation',
                obs(f'o{observation_index}', instrument='i0',
                    time=index * 1000, raw=str(slot))))
            observation_index += 1
        changes += extra
        batches.append(batch(f'b{index}', *current))
    rng = random.Random(2026)
    prefixes = list(range(301))
    rng.shuffle(prefixes)
    # Reserve boundary prefixes, including the empty and final view.
    chosen = [0, 100, 125, 150, 299, 300]
    chosen += [prefix for prefix in prefixes if prefix not in chosen][:54]
    queries = [query(f'q{k}', chosen[k]) for k in range(60)]
    return {'batches': batches, 'queries': queries}, changes, nodes, members


class ScaleBoundary(unittest.TestCase):
    @staticmethod
    def _observations_before(prefix):
        # Fifteen observation puts per batch, except fourteen in b100..149.
        certificate_batches = max(0, min(prefix, 150) - 100)
        return 15 * prefix - certificate_batches

    def test_300_batches_6000_changes_2000_member_budget_and_60_queries(self):
        source, change_count, node_count, member_count = _scale_input()
        self.assertEqual(change_count, 6000)
        self.assertEqual(member_count, 2000)
        self.assertLessEqual(node_count, 20000)
        self.assertEqual(len(source['batches']), 300)
        started = time.monotonic()
        result = analyze(source)
        elapsed = time.monotonic() - started
        self.assertEqual(len(result['batches']), 300)
        self.assertTrue(all(item['status'] == 'accepted'
                            for item in result['batches']))
        self.assertEqual(len(result['queries']), 60)
        prefixes = {item['as_of'] for item in result['queries']}
        self.assertEqual(len(prefixes), 60)
        latest = max(result['queries'], key=lambda item: item['as_of'])
        self.assertEqual(latest['as_of'], 300)
        self.assertEqual(len(latest['observations']), 4450)
        self.assertEqual(len(latest['certificates']), 50)
        for view in result['queries']:
            self.assertEqual(len(view['observations']),
                             self._observations_before(view['as_of']))
            self.assertTrue(all(
                item['status'] == 'available'
                for item in view['observations']))
        self.assertLess(elapsed, 60)

    def test_scale_python_allocation_in_budget(self):
        # ru_maxrss reports cgroup-wide peaks in some containers, so measure
        # the interpreter's own Python object allocations with tracemalloc.
        import tracemalloc
        source, _, _, _ = _scale_input()
        tracemalloc.start()
        analyze(source)
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        self.assertLess(peak / (1024 ** 2), 192)

    def test_scale_input_through_cli(self):
        source, _, _, _ = _scale_input()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = root / 'scale.json'
            out = root / 'report.json'
            src.write_text(json.dumps(source))
            started = time.monotonic()
            completed = subprocess.run(
                [sys.executable, '-m', 'reaction_batch.certify_cli',
                 '--source', str(src), '--output', str(out)],
                capture_output=True, cwd=ROOT,
                env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')})
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertLess(time.monotonic() - started, 120)
            report = json.loads(out.read_text())
            self.assertEqual(len(report['result']['queries']), 60)
            self.assertTrue(all(
                item['status'] == 'accepted'
                for item in report['result']['batches']))


if __name__ == '__main__':
    unittest.main()
