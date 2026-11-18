"""Independent checks for the reaction-peaks shared-drift reconciliation.

Small cases are verified against a brute-force oracle that enumerates every
injective explanation (peaks may be left unmatched), intersects the closed
drift interval of each explanation, and aggregates by size; nothing in the
oracle calls the production matching code.  The 40x40 scale case would need
40! enumerations and is checked directly against the documented expectation.
"""

import itertools
import json
import os
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from reaction_batch.chemistry import formula
from reaction_batch.peaks import analyze

ROOT = Path(__file__).resolve().parents[1]


def reference(reference_id, text, lower, upper):
    return {'reference_id': reference_id, 'formula': text,
            'window': {'lower': lower, 'upper': upper}}


def peak(peak_id, text, time):
    return {'peak_id': peak_id, 'formula': text, 'time': time}


def run(run_id, lower, upper, peaks):
    return {'run_id': run_id,
            'drift': {'lower': lower, 'upper': upper},
            'peaks': peaks}


def batch(references, runs):
    return {'references': references, 'runs': runs}


def merge(intervals):
    """Merge a set of closed Fraction intervals, sorted by lower."""
    merged = []
    for lower, upper in sorted(intervals):
        if merged and lower <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], upper))
        else:
            merged.append((lower, upper))
    return merged


def oracle(references, one_run):
    """Brute-force K, optimal drift set, supports and freeable flags.

    Every injective partial assignment is an explanation; its feasible drift
    set is the intersection of its edge intervals with the run drift domain.
    """
    peaks = one_run['peaks']
    lower_d = Fraction(one_run['drift']['lower'])
    upper_d = Fraction(one_run['drift']['upper'])
    edges = {}
    for i, p in enumerate(peaks):
        for j, r in enumerate(references):
            if formula(p['formula']) != formula(r['formula']):
                continue
            lower = Fraction(p['time']) - Fraction(r['window']['upper'])
            upper = Fraction(p['time']) - Fraction(r['window']['lower'])
            lower = max(lower, lower_d)
            upper = min(upper, upper_d)
            if lower <= upper:
                edges[(i, j)] = (lower, upper)
    explanations = []

    def visit(i, used, lower, upper, chosen):
        if i == len(peaks):
            explanations.append((frozenset(chosen), lower, upper))
            return
        visit(i + 1, used, lower, upper, chosen)
        for j in range(len(references)):
            if j in used or (i, j) not in edges:
                continue
            a, b = edges[(i, j)]
            lo, hi = max(lower, a), min(upper, b)
            if lo <= hi:
                visit(i + 1, used | {j}, lo, hi, chosen + [(i, j)])

    visit(0, set(), lower_d, upper_d, [])
    best = max(len(chosen) for chosen, _lo, _hi in explanations)
    optimal = [(lo, hi) for chosen, lo, hi in explanations
               if len(chosen) == best]
    support = {i: set() for i in range(len(peaks))}
    freeable = {i: False for i in range(len(peaks))}
    for chosen, lo, hi in explanations:
        if len(chosen) != best:
            continue
        covered = {i for i, _j in chosen}
        for i, j in chosen:
            support[i].add(j)
        for i in range(len(peaks)):
            if i not in covered:
                freeable[i] = True
    return best, merge(optimal), support, freeable, explanations


def oracle_witness(references, one_run, explanations, best, drift):
    """Lexicographically smallest assignment vector among K-explanations."""
    peaks = one_run['peaks']
    reference_ids = [r['reference_id'] for r in references]
    order = sorted(range(len(peaks)), key=lambda i: peaks[i]['peak_id'])
    vectors = []
    for chosen, lo, hi in explanations:
        if len(chosen) != best or not lo <= drift <= hi:
            continue
        mapping = dict(chosen)
        vectors.append(tuple(mapping.get(i) for i in order))

    def key(vector):
        return tuple((1, '') if v is None else (0, reference_ids[v])
                     for v in vector)

    winner = min(vectors, key=key)
    return [(peaks[i]['peak_id'],
             None if v is None else reference_ids[v])
            for i, v in zip(order, winner)]


def check_against_oracle(test, references, one_run, report_run):
    best, optimal, support, freeable, explanations = oracle(references, one_run)
    test.assertEqual(report_run['maximum_matches'], best)
    produced = [(Fraction(item['lower']), Fraction(item['upper']))
                for item in report_run['best_drifts']]
    test.assertEqual(produced, optimal)
    peaks = one_run['peaks']
    order = sorted(range(len(peaks)), key=lambda i: peaks[i]['peak_id'])
    reference_ids = [r['reference_id'] for r in references]
    expected_peaks = [
        {'peak_id': peaks[i]['peak_id'],
         'reference_ids': sorted(reference_ids[j] for j in support[i]),
         'can_be_unmatched': freeable[i]}
        for i in order]
    test.assertEqual(report_run['peaks'], expected_peaks)
    drift = optimal[0][0] if optimal else Fraction(one_run['drift']['lower'])
    test.assertEqual(report_run['witness']['drift'], str(drift))
    expected_witness = oracle_witness(
        references, one_run, explanations, best, drift)
    produced_witness = [(item['peak_id'], item['reference_id'])
                        for item in report_run['witness']['assignments']]
    test.assertEqual(produced_witness, expected_witness)
    # The witness itself is a legal K-explanation at its drift.
    assigned = [(item['peak_id'], item['reference_id'])
                for item in report_run['witness']['assignments']
                if item['reference_id'] is not None]
    test.assertEqual(len(assigned), best)
    test.assertEqual(len({r for _p, r in assigned}), len(assigned))
    by_id = {r['reference_id']: r for r in references}
    moment = Fraction(report_run['witness']['drift'])
    for peak_id, reference_id in assigned:
        entry = by_id[reference_id]
        observed = next(p for p in peaks if p['peak_id'] == peak_id)
        test.assertEqual(formula(entry['formula']), formula(observed['formula']))
        test.assertLessEqual(
            Fraction(entry['window']['lower']) + moment,
            Fraction(observed['time']))
        test.assertLessEqual(
            Fraction(observed['time']),
            Fraction(entry['window']['upper']) + moment)


class AcceptanceExamples(unittest.TestCase):
    def test_single_point_drift_holds_both_edges(self):
        data = batch(
            [reference('a', 'A', '0', '1'), reference('b', 'B', '0', '1')],
            [run('r', '0', '2', [peak('p', 'A', '1'), peak('q', 'B', '2')])])
        check_against_oracle(self, data['references'], data['runs'][0],
                             analyze(data)['runs'][0])
        out = analyze(data)['runs'][0]
        self.assertEqual(out['maximum_matches'], 2)
        self.assertEqual(out['best_drifts'], [{'lower': '1', 'upper': '1'}])
        self.assertEqual(out['witness']['drift'], '1')
        self.assertEqual(
            [(item['peak_id'], item['reference_id'])
             for item in out['witness']['assignments']],
            [('p', 'a'), ('q', 'b')])

    def test_two_separated_optimal_intervals(self):
        data = batch(
            [reference('a', 'A', '0', '1'), reference('b', 'A', '3', '4')],
            [run('r', '-1', '3', [peak('p', 'A', '3')])])
        out = analyze(data)['runs'][0]
        check_against_oracle(self, data['references'], data['runs'][0], out)
        self.assertEqual(out['maximum_matches'], 1)
        self.assertEqual(out['best_drifts'],
                         [{'lower': '-1', 'upper': '0'},
                          {'lower': '2', 'upper': '3'}])
        self.assertEqual(out['peaks'],
                         [{'peak_id': 'p', 'reference_ids': ['a', 'b'],
                           'can_be_unmatched': False}])
        self.assertEqual(out['witness']['drift'], '-1')
        self.assertEqual(out['witness']['assignments'],
                         [{'peak_id': 'p', 'reference_id': 'b'}])

    def test_contending_edge_is_never_supported(self):
        data = batch(
            [reference('a', 'A', '0', '1'), reference('b', 'A', '1', '1')],
            [run('r', '0', '0', [peak('p', 'A', '0'), peak('q', 'A', '1')])])
        out = analyze(data)['runs'][0]
        check_against_oracle(self, data['references'], data['runs'][0], out)
        self.assertEqual(out['maximum_matches'], 2)
        # q/a satisfies the window on its own but can never join a 2-match.
        self.assertEqual(out['peaks'],
                         [{'peak_id': 'p', 'reference_ids': ['a'],
                           'can_be_unmatched': False},
                          {'peak_id': 'q', 'reference_ids': ['b'],
                           'can_be_unmatched': False}])

    def test_one_reference_two_peaks_both_freeable(self):
        data = batch(
            [reference('a', 'A', '0', '0')],
            [run('r', '0', '2', [peak('p', 'A', '0'), peak('q', 'A', '0')])])
        out = analyze(data)['runs'][0]
        check_against_oracle(self, data['references'], data['runs'][0], out)
        self.assertEqual(out['maximum_matches'], 1)
        self.assertEqual(out['best_drifts'], [{'lower': '0', 'upper': '0'}])
        self.assertEqual(out['peaks'],
                         [{'peak_id': 'p', 'reference_ids': ['a'],
                           'can_be_unmatched': True},
                          {'peak_id': 'q', 'reference_ids': ['a'],
                           'can_be_unmatched': True}])
        self.assertEqual(out['witness']['assignments'],
                         [{'peak_id': 'p', 'reference_id': 'a'},
                          {'peak_id': 'q', 'reference_id': None}])

    def test_window_change_moves_other_group_support(self):
        references = [reference('a', 'A', '0', '0'), reference('b', 'A', '1', '1'),
                      reference('c', 'B', '0', '0')]
        peaks = [peak('p', 'A', '1'), peak('q', 'B', '1')]
        first = analyze(batch(references, [run('r', '0', '1', peaks)]))['runs'][0]
        self.assertEqual(first['maximum_matches'], 2)
        self.assertEqual(first['best_drifts'], [{'lower': '1', 'upper': '1'}])
        self.assertEqual(first['peaks'][0],
                         {'peak_id': 'p', 'reference_ids': ['a'],
                          'can_be_unmatched': False})
        widened = [references[0], references[1], reference('c', 'B', '0', '1')]
        second_data = batch(widened, [run('r', '0', '1', peaks)])
        second = analyze(second_data)['runs'][0]
        check_against_oracle(self, widened, second_data['runs'][0], second)
        self.assertEqual(second['maximum_matches'], 2)
        self.assertEqual(second['best_drifts'],
                         [{'lower': '0', 'upper': '0'},
                          {'lower': '1', 'upper': '1'}])
        self.assertEqual(second['peaks'][0],
                         {'peak_id': 'p', 'reference_ids': ['a', 'b'],
                          'can_be_unmatched': False})


class ScaleAndCompleteness(unittest.TestCase):
    def test_forty_by_forty_full_support(self):
        data = json.loads((ROOT / 'fixtures' / 'peaks-scale-40.json').read_text())
        out = analyze(data)['runs'][0]
        self.assertEqual(out['maximum_matches'], 40)
        self.assertEqual(out['best_drifts'], [{'lower': '0', 'upper': '1'}])
        all_refs = sorted(r['reference_id'] for r in data['references'])
        for entry in out['peaks']:
            self.assertEqual(entry['reference_ids'], all_refs)
            self.assertIs(entry['can_be_unmatched'], False)
        assignments = out['witness']['assignments']
        self.assertEqual([item['peak_id'] for item in assignments],
                         sorted(p['peak_id'] for p in data['runs'][0]['peaks']))
        self.assertEqual([item['reference_id'] for item in assignments],
                         sorted(r['reference_id'] for r in data['references']))

    def test_many_equivalent_regions_stay_exact(self):
        # Distinct breakpoints per edge: the sweep visits every region.
        references = [reference('r%d' % j, 'A', str(j), str(j + 1))
                      for j in range(6)]
        peaks = [peak('p%d' % i, 'A', str(i + 1)) for i in range(6)]
        data = batch(references, [run('r', '-6', '7', peaks)])
        out = analyze(data)['runs'][0]
        check_against_oracle(self, references, data['runs'][0], out)


class FormulaAndIdentity(unittest.TestCase):
    def test_equivalent_formula_writings_are_compatible(self):
        data = batch(
            [reference('w1', 'H2O', '0', '1'), reference('w2', 'OH2', '0', '1')],
            [run('r', '0', '1', [peak('x', 'O(H)2', '1')])])
        out = analyze(data)['runs'][0]
        check_against_oracle(self, data['references'], data['runs'][0], out)
        self.assertEqual(out['peaks'][0]['reference_ids'], ['w1', 'w2'])

    def test_duplicate_formula_texts_keep_distinct_identities(self):
        data = batch(
            [reference('a', 'H2O', '0', '1'), reference('b', 'H2O', '0', '1')],
            [run('r', '0', '1', [peak('x', 'H2O', '1'), peak('y', 'H2O', '1')])])
        out = analyze(data)['runs'][0]
        check_against_oracle(self, data['references'], data['runs'][0], out)
        self.assertEqual(out['maximum_matches'], 2)
        for entry in out['peaks']:
            self.assertEqual(entry['reference_ids'], ['a', 'b'])
            self.assertIs(entry['can_be_unmatched'], False)

    def test_incompatible_elements_never_match(self):
        data = batch(
            [reference('a', 'C', '0', '10')],
            [run('r', '0', '10', [peak('x', 'N', '5')])])
        out = analyze(data)['runs'][0]
        self.assertEqual(out['maximum_matches'], 0)
        self.assertEqual(out['peaks'][0]['reference_ids'], [])
        self.assertIs(out['peaks'][0]['can_be_unmatched'], True)

    def test_nested_formula_counts(self):
        data = batch(
            [reference('a', 'Al2(SO4)3', '2', '3')],
            [run('r', '0', '1', [peak('x', '(SO4)3Al2', '3.5')])])
        out = analyze(data)['runs'][0]
        check_against_oracle(self, data['references'], data['runs'][0], out)
        self.assertEqual(out['maximum_matches'], 1)


class BoundaryConditions(unittest.TestCase):
    def test_no_compatible_edge_gives_full_domain(self):
        data = batch(
            [reference('a', 'A', '0', '1')],
            [run('r', '-2', '5', [peak('x', 'B', '0'), peak('y', 'A', '9')])])
        out = analyze(data)['runs'][0]
        self.assertEqual(out['maximum_matches'], 0)
        self.assertEqual(out['best_drifts'], [{'lower': '-2', 'upper': '5'}])
        self.assertEqual(out['peaks'],
                         [{'peak_id': 'x', 'reference_ids': [],
                           'can_be_unmatched': True},
                          {'peak_id': 'y', 'reference_ids': [],
                           'can_be_unmatched': True}])
        self.assertEqual(out['witness'],
                         {'drift': '-2',
                          'assignments': [{'peak_id': 'x', 'reference_id': None},
                                          {'peak_id': 'y', 'reference_id': None}]})

    def test_empty_peaks(self):
        data = batch([reference('a', 'A', '0', '1')],
                     [run('r', '-1', '2', [])])
        out = analyze(data)['runs'][0]
        self.assertEqual(out['maximum_matches'], 0)
        self.assertEqual(out['best_drifts'], [{'lower': '-1', 'upper': '2'}])
        self.assertEqual(out['peaks'], [])
        self.assertEqual(out['witness'], {'drift': '-1', 'assignments': []})

    def test_negative_zero_drift_is_normalised(self):
        data = batch([reference('a', 'A', '0', '0')],
                     [run('r', '-0', '1', [peak('x', 'A', '0')])])
        out = analyze(data)['runs'][0]
        self.assertEqual(out['best_drifts'], [{'lower': '0', 'upper': '0'}])
        self.assertEqual(out['witness']['drift'], '0')
        for item in out['best_drifts'] + [out['witness']]:
            for value in item.values() if isinstance(item, dict) else ():
                self.assertNotEqual(value, '-0')

    def test_fractional_drift_bounds_are_reduced(self):
        data = batch([reference('a', 'A', '0.3', '0.3')],
                     [run('r', '0', '1', [peak('x', 'A', '0.5')])])
        out = analyze(data)['runs'][0]
        self.assertEqual(out['best_drifts'], [{'lower': '1/5', 'upper': '1/5'}])
        self.assertEqual(out['witness']['drift'], '1/5')

    def test_zero_width_window_and_touching_endpoint(self):
        data = batch(
            [reference('a', 'A', '1', '1'), reference('b', 'A', '2', '2')],
            [run('r', '0', '2', [peak('x', 'A', '2')])])
        out = analyze(data)['runs'][0]
        check_against_oracle(self, data['references'], data['runs'][0], out)
        # x/a needs drift 1; x/b needs drift 0; both singletons optimal.
        self.assertEqual(out['best_drifts'],
                         [{'lower': '0', 'upper': '0'},
                          {'lower': '1', 'upper': '1'}])
        self.assertEqual(out['peaks'][0]['reference_ids'], ['a', 'b'])

    def test_closed_endpoint_beats_open_neighbourhood(self):
        # Only the exact drift 1 aligns both pairs; either side matches one.
        data = batch(
            [reference('a', 'A', '0', '1'), reference('b', 'B', '0', '1')],
            [run('r', '0', '2', [peak('p', 'A', '1'), peak('q', 'B', '2')])])
        out = analyze(data)['runs'][0]
        self.assertEqual(out['best_drifts'], [{'lower': '1', 'upper': '1'}])

    def test_output_types_are_json_scalars(self):
        data = json.loads((ROOT / 'fixtures' / 'peaks-example2.json').read_text())
        out = analyze(data)['runs'][0]
        self.assertIsInstance(out['maximum_matches'], int)
        for interval in out['best_drifts']:
            self.assertIsInstance(interval['lower'], str)
            self.assertIsInstance(interval['upper'], str)
        for entry in out['peaks']:
            self.assertIsInstance(entry['can_be_unmatched'], bool)
            for reference_id in entry['reference_ids']:
                self.assertIsInstance(reference_id, str)
        self.assertIsInstance(out['witness']['drift'], str)


class Reordering(unittest.TestCase):
    def test_permuted_source_arrays_keep_identity_results(self):
        references = [reference('a', 'A', '0', '1'), reference('b', 'A', '3', '4'),
                      reference('c', 'B', '0', '2')]
        runs = [run('r1', '-1', '3', [peak('p', 'A', '3'), peak('q', 'B', '1')]),
                run('r2', '0', '1', [peak('p', 'A', '4'), peak('z', 'B', '2')])]
        base = analyze(batch(references, runs))
        shuffled = batch(
            [references[2], references[0], references[1]],
            [run('r2', '0', '1', [runs[1]['peaks'][1], runs[1]['peaks'][0]]),
             run('r1', '-1', '3', [runs[0]['peaks'][1], runs[0]['peaks'][0]])])
        second = analyze(shuffled)
        self.assertEqual([r['run_id'] for r in second['runs']], ['r2', 'r1'])
        by_id = {r['run_id']: r for r in base['runs']}
        for entry in second['runs']:
            self.assertEqual(entry, by_id[entry['run_id']])

    def test_peak_ids_may_repeat_across_runs(self):
        data = batch(
            [reference('a', 'A', '0', '1')],
            [run('r1', '0', '1', [peak('p', 'A', '1')]),
             run('r2', '0', '1', [peak('p', 'A', '0')])])
        out = analyze(data)
        self.assertEqual([r['maximum_matches'] for r in out['runs']], [1, 1])


class RandomisedCrossCheck(unittest.TestCase):
    def test_random_small_batches_against_enumeration(self):
        import random
        random.seed(20261118)
        formulas = ['A', 'B', 'H2O', 'OH2', 'C']
        checked = 0
        for _ in range(120):
            references = []
            for j in range(random.randint(1, 4)):
                low = random.randint(0, 4)
                references.append(reference(
                    'r%d' % j, random.choice(formulas),
                    str(low), str(low + random.randint(0, 3))))
            runs = []
            for q in range(random.randint(1, 2)):
                low = random.randint(-3, 2)
                peaks = []
                for i in range(random.randint(0, 4)):
                    peaks.append(peak(
                        'p%d' % i, random.choice(formulas),
                        str(random.randint(0, 6))))
                runs.append(run('q%d' % q, str(low),
                                str(low + random.randint(0, 4)), peaks))
            data = batch(references, runs)
            report = analyze(data)
            self.assertEqual([r['run_id'] for r in report['runs']],
                             [r['run_id'] for r in runs])
            for produced, one_run in zip(report['runs'], runs):
                check_against_oracle(self, references, one_run, produced)
                checked += 1
        self.assertGreater(checked, 150)


class InputValidation(unittest.TestCase):
    def base(self):
        return batch([reference('a', 'A', '0', '1')],
                     [run('r', '0', '1', [peak('p', 'A', '1')])])

    def assert_rejected(self, mutate):
        data = self.base()
        mutate(data)
        with self.assertRaises(ValueError):
            analyze(data)

    def test_unknown_fields(self):
        def top(data):
            data['extra'] = 1
        def ref(data):
            data['references'][0]['extra'] = 1
        def window(data):
            data['references'][0]['window']['extra'] = 1
        def run_level(data):
            data['runs'][0]['extra'] = 1
        def drift(data):
            data['runs'][0]['drift']['extra'] = 1
        def peak_level(data):
            data['runs'][0]['peaks'][0]['extra'] = 1
        def missing(data):
            del data['runs'][0]['peaks'][0]['time']
        for mutate in (top, ref, window, run_level, drift, peak_level, missing):
            with self.subTest(mutate=mutate.__name__):
                self.assert_rejected(mutate)

    def test_duplicate_identifiers(self):
        def ref(data):
            data['references'].append(reference('a', 'B', '0', '1'))
        def run_id(data):
            data['runs'].append(run('r', '0', '1', []))
        def peak_id(data):
            data['runs'][0]['peaks'].append(peak('p', 'A', '0'))
        for mutate in (ref, run_id, peak_id):
            with self.subTest(mutate=mutate.__name__):
                self.assert_rejected(mutate)

    def test_invalid_literals(self):
        cases = [
            lambda d: d['references'][0].__setitem__('formula', 'H0'),
            lambda d: d['references'][0].__setitem__('formula', 'h2o'),
            lambda d: d['references'][0].__setitem__('formula', 12),
            lambda d: d['references'][0]['window'].__setitem__('lower', '02'),
            lambda d: d['references'][0]['window'].__setitem__('lower', 0),
            lambda d: d['references'][0]['window'].__setitem__('lower', True),
            lambda d: d['references'][0]['window'].__setitem__('lower', '-1'),
            lambda d: d['references'][0]['window'].__setitem__('lower', '1e3'),
            lambda d: d['references'][0]['window'].__setitem__('lower', ' 1'),
            lambda d: d['references'][0]['window'].__setitem__('lower', '2'),
            lambda d: d['runs'][0]['drift'].__setitem__('lower', '+1'),
            lambda d: d['runs'][0]['drift'].__setitem__('lower', '--1'),
            lambda d: d['runs'][0]['drift'].__setitem__('lower', '-'),
            lambda d: d['runs'][0]['drift'].__setitem__('lower', 0),
            lambda d: d['runs'][0]['drift'].__setitem__('lower', '2'),
            lambda d: d['runs'][0]['peaks'][0].__setitem__('time', '-1'),
            lambda d: d['runs'][0]['peaks'][0].__setitem__('time', 1),
            lambda d: d['runs'][0]['peaks'][0].__setitem__('peak_id', ' x'),
            lambda d: d['runs'][0]['peaks'][0].__setitem__('peak_id', ''),
        ]
        for index, mutate in enumerate(cases):
            with self.subTest(case=index):
                self.assert_rejected(mutate)

    def test_size_limits(self):
        def no_refs(data):
            data['references'] = []
        def many_refs(data):
            data['references'] = [reference('x%d' % i, 'A', '0', '1')
                                  for i in range(41)]
        def no_runs(data):
            data['runs'] = []
        def many_runs(data):
            data['runs'] = [run('q%d' % i, '0', '1', []) for i in range(9)]
        def many_peaks(data):
            data['runs'][0]['peaks'] = [peak('p%d' % i, 'A', '1')
                                        for i in range(41)]
        for mutate in (no_refs, many_refs, no_runs, many_runs, many_peaks):
            with self.subTest(mutate=mutate.__name__):
                self.assert_rejected(mutate)

    def test_late_run_error_rejects_whole_batch(self):
        data = self.base()
        data['runs'].append(run('second', '0', '1', [peak('p', 'A', '1')]))
        data['runs'][1]['peaks'][0]['time'] = 'not-a-number'
        with self.assertRaises(ValueError):
            analyze(data)


class CliRegression(unittest.TestCase):
    def run_cli(self, source, output, module='reaction_batch.peaks_cli'):
        return subprocess.run(
            [sys.executable, '-m', module,
             '--source', str(source), '--output', str(output)],
            capture_output=True, cwd=ROOT,
            env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')})

    def test_fixtures_end_to_end(self):
        import hashlib
        for name in ('peaks-example1.json', 'peaks-example2.json',
                     'peaks-linkage.json', 'peaks-scale-40.json'):
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as tmp:
                    output = Path(tmp) / 'report.json'
                    source = ROOT / 'fixtures' / name
                    completed = self.run_cli(source, output)
                    self.assertEqual(completed.returncode, 0,
                                     completed.stderr[:400])
                    report = json.loads(output.read_text())
                    self.assertEqual(report['schema_version'], 1)
                    self.assertEqual(
                        report['source_sha256'],
                        hashlib.sha256(source.read_bytes()).hexdigest())
                    self.assertEqual(set(report), {'schema_version',
                                                   'source_sha256', 'result'})
                    self.assertEqual(set(report['result']), {'runs'})
                    data = json.loads(source.read_text())
                    self.assertEqual([r['run_id'] for r in report['result']['runs']],
                                     [r['run_id'] for r in data['runs']])
                    for produced, one_run in zip(report['result']['runs'],
                                                 data['runs']):
                        self.assertEqual(
                            set(produced), {'run_id', 'maximum_matches',
                                            'best_drifts', 'peaks', 'witness'})
                        if name != 'peaks-scale-40.json':
                            # The enumeration oracle is exponential; the 40x40
                            # fixture is checked directly in ScaleAndCompleteness.
                            check_against_oracle(
                                self, data['references'], one_run, produced)

    def test_rerun_is_byte_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            first = Path(tmp) / 'a.json'
            second = Path(tmp) / 'b.json'
            source = ROOT / 'fixtures' / 'peaks-scale-40.json'
            self.assertEqual(self.run_cli(source, first).returncode, 0)
            self.assertEqual(self.run_cli(source, second).returncode, 0)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_invalid_input_preserves_old_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'bad.json'
            output = Path(tmp) / 'out.json'
            source.write_text('{"references": [], "runs": []}')
            output.write_text('old report')
            completed = self.run_cli(source, output)
            self.assertEqual(completed.returncode, 2)
            self.assertEqual(output.read_text(), 'old report')

    def test_late_run_error_preserves_old_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'bad.json'
            output = Path(tmp) / 'out.json'
            source.write_text(json.dumps({
                'references': [reference('a', 'A', '0', '1')],
                'runs': [run('ok', '0', '1', [peak('p', 'A', '1')]),
                         run('bad', '2', '1', [])]}))
            output.write_text('old report')
            completed = self.run_cli(source, output)
            self.assertEqual(completed.returncode, 2)
            self.assertEqual(output.read_text(), 'old report')

    def test_output_guards(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = ROOT / 'fixtures' / 'peaks-example1.json'
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
            source.write_text(
                '{"references": [], "references": [], "runs": []}')
            output = Path(tmp) / 'out.json'
            completed = self.run_cli(source, output)
            self.assertEqual(completed.returncode, 2)
            self.assertFalse(output.exists())

    def test_existing_commands_still_work(self):
        cases = [('reaction_batch.cli', 'input.json'),
                 ('reaction_batch.lots', 'lots-example.json'),
                 ('reaction_batch.assay_cli', 'assay-example1.json')]
        with tempfile.TemporaryDirectory() as tmp:
            for index, (module, fixture) in enumerate(cases):
                with self.subTest(module=module):
                    source = ROOT / 'fixtures' / fixture
                    output = Path(tmp) / ('out%d.json' % index)
                    completed = self.run_cli(source, output, module=module)
                    self.assertEqual(completed.returncode, 0,
                                     completed.stderr[:400])
                    report = json.loads(output.read_text())
                    self.assertEqual(report['schema_version'], 1)


if __name__ == '__main__':
    unittest.main()
