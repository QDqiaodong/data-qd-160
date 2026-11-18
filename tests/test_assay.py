"""Independent checks for the assay report.

Every certificate is verified straight from the documented canonical system
``A x <= b``, rebuilt here from the input formulas and the contract's row
order.  Nothing in this file calls the simplex; the expected values are either
hand-derived or obtained by enumerating the vertices of tiny polyhedra.
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

from reaction_batch.assay import analyze
from reaction_batch.chemistry import formula

ROOT = Path(__file__).resolve().parents[1]


def species_pair(species_id, text):
    return {'species_id': species_id, 'formula': text}


def measurement(measurement_id, element, lower, upper):
    return {'measurement_id': measurement_id, 'element': element,
            'lower': lower, 'upper': upper}


def sample(sample_id, lower, upper, measurements):
    return {'sample_id': sample_id,
            'total_moles': {'lower': lower, 'upper': upper},
            'measurements': measurements}


def canonical(species, one_sample):
    """Rebuild ``A`` and ``b`` exactly as the contract defines them."""
    n = len(species)
    compositions = [formula(item['formula']) for item in species]
    rows = []
    bounds = []
    for j in range(n):
        rows.append([Fraction(-1) if k == j else Fraction(0) for k in range(n)])
        bounds.append(Fraction(0))
    rows.append([Fraction(1)] * n)
    bounds.append(Fraction(one_sample['total_moles']['upper']))
    rows.append([Fraction(-1)] * n)
    bounds.append(-Fraction(one_sample['total_moles']['lower']))
    for item in one_sample['measurements']:
        coefficients = [Fraction(composition.get(item['element'], 0))
                        for composition in compositions]
        rows.append(list(coefficients))
        bounds.append(Fraction(item['upper']))
        rows.append([-value for value in coefficients])
        bounds.append(-Fraction(item['lower']))
    return rows, bounds


def dot(left, right):
    return sum((a * b for a, b in zip(left, right)), Fraction(0))


def check_report(test, data, report):
    """Verify every numeric claim of a report against the canonical system."""
    species = data['species']
    n = len(species)
    test.assertEqual(report['result']['species_order'],
                     [item['species_id'] for item in species])
    samples = report['result']['samples']
    test.assertEqual([item['sample_id'] for item in samples],
                     [item['sample_id'] for item in data['samples']])
    for out, one_sample in zip(samples, data['samples']):
        rows, bounds = canonical(species, one_sample)
        row_count = len(rows)
        test.assertEqual(row_count, n + 2 + 2 * len(one_sample['measurements']))
        if out['status'] == 'inconsistent':
            test.assertNotIn('ranges', out)
            y = [Fraction(value) for value in out['contradiction']]
            test.assertEqual(len(y), row_count)
            test.assertTrue(all(value >= 0 for value in y))
            for j in range(n):
                column = [rows[i][j] for i in range(row_count)]
                test.assertEqual(dot(column, y), 0, 'Farkas column %d' % j)
            test.assertEqual(dot(bounds, y), -1)
            continue
        test.assertEqual(out['status'], 'consistent')
        test.assertNotIn('contradiction', out)
        test.assertEqual(len(out['ranges']), n)
        for index, entry in enumerate(out['ranges']):
            test.assertEqual(entry['species_id'], species[index]['species_id'])
            for label, sign in (('minimum', -1), ('maximum', 1)):
                block = entry[label]
                composition = [Fraction(value) for value in block['composition']]
                dual = [Fraction(value) for value in block['dual']]
                value = Fraction(block['value'])
                test.assertEqual(len(composition), n)
                test.assertEqual(len(dual), row_count)
                test.assertTrue(all(item >= 0 for item in composition))
                test.assertTrue(all(item >= 0 for item in dual))
                for i in range(row_count):
                    test.assertLessEqual(dot(rows[i], composition), bounds[i])
                test.assertEqual(composition[index], value)
                for j in range(n):
                    column = [rows[i][j] for i in range(row_count)]
                    target = Fraction(sign) if j == index else Fraction(0)
                    test.assertEqual(dot(column, dual), target)
                test.assertEqual(dot(bounds, dual), sign * value)


def solve_square(matrix, rhs):
    size = len(matrix)
    augmented = [list(matrix[i]) + [rhs[i]] for i in range(size)]
    for column in range(size):
        pivot = next((i for i in range(column, size) if augmented[i][column]), None)
        if pivot is None:
            return None
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        factor = augmented[column][column]
        augmented[column] = [x / factor for x in augmented[column]]
        for i in range(size):
            if i != column and augmented[i][column]:
                ratio = augmented[i][column]
                augmented[i] = [a - ratio * b
                                for a, b in zip(augmented[i], augmented[column])]
    return [augmented[i][size] for i in range(size)]


def vertex_ranges(rows, bounds):
    """Exact min/max of every coordinate over basic feasible solutions.

    Used only for tiny polyhedra as an independent oracle for optimality.
    """
    row_count = len(rows)
    n = len(rows[0])
    minimum = [None] * n
    maximum = [None] * n
    feasible = False
    for subset in itertools.combinations(range(row_count), n):
        point = solve_square([rows[i] for i in subset], [bounds[i] for i in subset])
        if point is None or any(value < 0 for value in point):
            continue
        if any(dot(rows[i], point) > bounds[i] for i in range(row_count)):
            continue
        feasible = True
        for j in range(n):
            if minimum[j] is None or point[j] < minimum[j]:
                minimum[j] = point[j]
            if maximum[j] is None or point[j] > maximum[j]:
                maximum[j] = point[j]
    return feasible, minimum, maximum


class AcceptanceExamples(unittest.TestCase):
    def test_water_only_is_unique(self):
        data = {'species': [species_pair('h', 'H2'), species_pair('o', 'O2'),
                            species_pair('w', 'H2O')],
                'samples': [sample('s1', '2', '2', [
                    measurement('mH', 'H', '3', '3'),
                    measurement('mO', 'O', '1.5', '1.5')])]}
        report = analyze(data)
        check_report(self, data, {'result': report})
        ranges = report['samples'][0]['ranges']
        self.assertEqual(
            [(entry['species_id'], entry['minimum']['value'],
              entry['maximum']['value']) for entry in ranges],
            [('h', '1', '1'), ('o', '1/2', '1/2'), ('w', '1/2', '1/2')])
        for entry in ranges:
            self.assertEqual(entry['minimum']['composition'],
                             ['1', '1/2', '1/2'])
            self.assertEqual(entry['maximum']['composition'],
                             ['1', '1/2', '1/2'])

    def test_two_formula_writings_each_endpoint_has_own_witness(self):
        data = {'species': [species_pair('a', 'H2O'),
                            species_pair('b', 'OH2')],
                'samples': [sample('s2', '1', '1', [
                    measurement('mH', 'H', '2', '2'),
                    measurement('mO', 'O', '1', '1')])]}
        report = analyze(data)
        check_report(self, data, {'result': report})
        ranges = report['samples'][0]['ranges']
        for entry in ranges:
            self.assertEqual(entry['minimum']['value'], '0')
            self.assertEqual(entry['maximum']['value'], '1')
        self.assertEqual(ranges[0]['minimum']['composition'], ['0', '1'])
        self.assertEqual(ranges[0]['maximum']['composition'], ['1', '0'])
        self.assertEqual(ranges[1]['minimum']['composition'], ['1', '0'])
        self.assertEqual(ranges[1]['maximum']['composition'], ['0', '1'])

    def test_nanomole_mismatch_is_inconsistent(self):
        data = {'species': [species_pair('w', 'H2O')],
                'samples': [sample('s3', '1', '1', [
                    measurement('mH', 'H', '2', '2'),
                    measurement('mO', 'O', '1.000000001', '1.000000001')])]}
        report = analyze(data)
        check_report(self, data, {'result': report})
        entry = report['samples'][0]
        self.assertEqual(entry['status'], 'inconsistent')
        self.assertNotIn('ranges', entry)
        self.assertTrue(all(Fraction(value) >= 0
                            for value in entry['contradiction']))


class EdgeCases(unittest.TestCase):
    def test_total_upper_zero_feasible(self):
        data = {'species': [species_pair('a', 'H2'), species_pair('b', 'O2')],
                'samples': [sample('z', '0', '0', [])]}
        report = analyze(data)
        check_report(self, data, {'result': report})
        self.assertEqual(report['samples'][0]['status'], 'consistent')
        for entry in report['samples'][0]['ranges']:
            self.assertEqual(entry['minimum']['value'], '0')
            self.assertEqual(entry['maximum']['value'], '0')

    def test_total_upper_zero_infeasible(self):
        data = {'species': [species_pair('a', 'H2')],
                'samples': [sample('z', '0', '0', [
                    measurement('m', 'H', '1', '1')])]}
        report = analyze(data)
        check_report(self, data, {'result': report})
        self.assertEqual(report['samples'][0]['status'], 'inconsistent')

    def test_empty_measurements_leave_total_range(self):
        data = {'species': [species_pair('a', 'H2')],
                'samples': [sample('z', '1', '2', [])]}
        report = analyze(data)
        check_report(self, data, {'result': report})
        entry = report['samples'][0]['ranges'][0]
        self.assertEqual(entry['minimum']['value'], '1')
        self.assertEqual(entry['maximum']['value'], '2')

    def test_element_absent_from_every_species(self):
        zero = {'species': [species_pair('a', 'H2')],
                'samples': [sample('z', '1', '1', [
                    measurement('m', 'Xe', '0', '0')])]}
        check_report(self, zero, {'result': analyze(zero)})
        self.assertEqual(analyze(zero)['samples'][0]['status'], 'consistent')
        nonzero = {'species': [species_pair('a', 'H2')],
                   'samples': [sample('z', '1', '1', [
                       measurement('m', 'Xe', '1', '1')])]}
        self.assertEqual(analyze(nonzero)['samples'][0]['status'],
                         'inconsistent')
        check_report(self, nonzero, {'result': analyze(nonzero)})

    def test_repeated_compatible_intervals(self):
        data = {'species': [species_pair('a', 'H2')],
                'samples': [sample('z', '1', '1', [
                    measurement('m1', 'H', '1', '3'),
                    measurement('m2', 'H', '2', '2')])]}
        report = analyze(data)
        check_report(self, data, {'result': report})
        entry = report['samples'][0]['ranges'][0]
        self.assertEqual(entry['minimum']['value'], '1')
        self.assertEqual(entry['maximum']['value'], '1')

    def test_repeated_contradictory_intervals(self):
        data = {'species': [species_pair('a', 'H2')],
                'samples': [sample('z', '1', '1', [
                    measurement('m1', 'H', '1', '1'),
                    measurement('m2', 'H', '3', '3')])]}
        report = analyze(data)
        check_report(self, data, {'result': report})
        self.assertEqual(report['samples'][0]['status'], 'inconsistent')

    def test_multiple_degrees_of_freedom(self):
        data = {'species': [species_pair('a', 'H2O'),
                            species_pair('b', 'OH2'),
                            species_pair('c', 'H2')],
                'samples': [sample('z', '0', '2', [
                    measurement('mH', 'H', '0', '4')])]}
        report = analyze(data)
        check_report(self, data, {'result': report})
        ranges = report['samples'][0]['ranges']
        self.assertEqual([entry['maximum']['value'] for entry in ranges],
                         ['2', '2', '2'])


class VertexOracle(unittest.TestCase):
    """Confirm bounds are optimal, not merely feasible, on tiny systems."""

    CASES = [
        ([species_pair('h', 'H2'), species_pair('o', 'O2'),
          species_pair('w', 'H2O')], sample('s', '1', '3', [
            measurement('a', 'H', '2', '4'),
            measurement('b', 'O', '1', '2'),
            measurement('c', 'H', '3', '3')])),
        ([species_pair('a', 'H2O'), species_pair('b', 'OH2'),
          species_pair('c', 'H2')], sample('s', '0', '2', [
            measurement('a', 'H', '0', '4'),
            measurement('b', 'O', '0', '3')])),
        ([species_pair('a', 'CH4'), species_pair('b', 'C2H6'),
          species_pair('c', 'O2')], sample('s', '0', '5', [
            measurement('a', 'C', '0', '6'),
            measurement('b', 'H', '0', '20')])),
    ]

    def test_vertices_match_reported_bounds(self):
        for species, one_sample in self.CASES:
            with self.subTest(species=species):
                data = {'species': species, 'samples': [one_sample]}
                report = analyze(data)
                rows, bounds = canonical(species, one_sample)
                feasible, minimum, maximum = vertex_ranges(rows, bounds)
                if not feasible:
                    self.assertEqual(report['samples'][0]['status'],
                                     'inconsistent')
                    continue
                ranges = report['samples'][0]['ranges']
                for j, entry in enumerate(ranges):
                    self.assertEqual(Fraction(entry['minimum']['value']),
                                     minimum[j])
                    self.assertEqual(Fraction(entry['maximum']['value']),
                                     maximum[j])


class RandomisedCrossCheck(unittest.TestCase):
    def test_random_small_systems_against_vertex_enumeration(self):
        import random
        random.seed(4242)
        elements = ['H', 'C', 'N', 'O']
        formulas = ['H2', 'O2', 'H2O', 'CO', 'CO2', 'CH4', 'N2', 'NO',
                    'NH3', 'C2H6', 'H2SO4', 'SO2']
        checked = 0
        for _ in range(60):
            n = random.randint(1, 3)
            chosen = random.sample(formulas, n)
            species = [species_pair('s%d' % i, text)
                       for i, text in enumerate(chosen)]
            measurements = []
            for q in range(random.randint(0, 3)):
                low = random.randint(0, 3)
                measurements.append(measurement(
                    'm%d' % q, random.choice(elements),
                    str(low), str(random.randint(low, 6))))
            low = random.randint(0, 3)
            one_sample = sample('s', str(low), str(random.randint(low, 5)),
                                measurements)
            data = {'species': species, 'samples': [one_sample]}
            report = analyze(data)
            check_report(self, data, {'result': report})
            rows, bounds = canonical(species, one_sample)
            feasible, minimum, maximum = vertex_ranges(rows, bounds)
            if not feasible:
                self.assertEqual(report['samples'][0]['status'],
                                 'inconsistent')
                continue
            self.assertEqual(report['samples'][0]['status'], 'consistent')
            for j, entry in enumerate(report['samples'][0]['ranges']):
                self.assertEqual(Fraction(entry['minimum']['value']),
                                 minimum[j])
                self.assertEqual(Fraction(entry['maximum']['value']),
                                 maximum[j])
            checked += 1
        self.assertGreater(checked, 20)


class NonUniqueComposition(unittest.TestCase):
    def test_same_composition_different_writings(self):
        data = {'species': [species_pair('a', 'H2O'),
                            species_pair('b', 'OH2')],
                'samples': [sample('s', '0', '2', [
                    measurement('mH', 'H', '0', '4')])]}
        report = analyze(data)
        check_report(self, data, {'result': report})
        first, second = report['samples'][0]['ranges']
        for entry in (first, second):
            self.assertEqual(entry['minimum']['value'], '0')
            self.assertEqual(entry['maximum']['value'], '2')
        self.assertNotEqual(first['minimum']['composition'],
                            first['maximum']['composition'])

    def test_distinct_formula_texts_are_not_merged(self):
        data = {'species': [species_pair('a', 'H2O'),
                            species_pair('b', 'OH2')],
                'samples': [sample('s', '1', '1', [
                    measurement('mH', 'H', '2', '2'),
                    measurement('mO', 'O', '1', '1')])]}
        self.assertEqual(len(analyze(data)['species_order']), 2)


class InputValidation(unittest.TestCase):
    def base(self):
        return {'species': [species_pair('a', 'H2')],
                'samples': [sample('s', '0', '1', [])]}

    def assert_rejected(self, mutate):
        data = self.base()
        mutate(data)
        with self.assertRaises(ValueError):
            analyze(data)

    def test_unknown_fields(self):
        for path in ('species', 'sample', 'total', 'measurement'):
            def mutate(data, path=path):
                if path == 'species':
                    data['species'][0]['extra'] = 1
                elif path == 'sample':
                    data['samples'][0]['extra'] = 1
                elif path == 'total':
                    data['samples'][0]['total_moles']['extra'] = 1
                else:
                    data['samples'][0]['measurements'] = [
                        measurement('m', 'H', '0', '1')]
                    data['samples'][0]['measurements'][0]['extra'] = 1
            with self.subTest(path=path):
                self.assert_rejected(mutate)

    def test_duplicate_identifiers(self):
        def duplicate_species(data):
            data['species'].append(species_pair('a', 'O2'))
        def duplicate_formula(data):
            data['species'].append(species_pair('b', 'H2'))
        def duplicate_sample(data):
            data['samples'].append(sample('s', '0', '1', []))
        def duplicate_measurement(data):
            data['samples'][0]['measurements'] = [
                measurement('m', 'H', '0', '1'),
                measurement('m', 'O', '0', '1')]
        for mutate in (duplicate_species, duplicate_formula,
                       duplicate_sample, duplicate_measurement):
            with self.subTest(mutate=mutate.__name__):
                self.assert_rejected(mutate)

    def test_invalid_literals(self):
        cases = [
            lambda data: data['species'].__setitem__(
                0, species_pair('a', 'H0')),
            lambda data: data['species'].__setitem__(
                0, species_pair('a', 'h2o')),
            lambda data: data['samples'][0].__setitem__(
                'measurements', [measurement('m', 'h', '0', '1')]),
            lambda data: data['samples'][0].__setitem__(
                'measurements', [measurement('m', 'H1', '0', '1')]),
            lambda data: data['samples'][0]['total_moles'].__setitem__(
                'lower', '1.5'),
            lambda data: data['samples'][0]['total_moles'].__setitem__(
                'lower', 0),
            lambda data: data['samples'][0]['total_moles'].__setitem__(
                'lower', '02'),
            lambda data: data['samples'][0]['total_moles'].__setitem__(
                'lower', '-1'),
            lambda data: data['samples'][0]['measurements'].append(
                measurement('m', 'H', '2', '1')),
            lambda data: data['samples'][0]['measurements'].append(
                measurement('m', 'H', True, '1')),
        ]
        for index, mutate in enumerate(cases):
            with self.subTest(case=index):
                self.assert_rejected(mutate)

    def test_size_limits(self):
        def too_many_species(data):
            data['species'] = [species_pair('x%d' % i, 'H2')
                               for i in range(61)]
        def too_many_samples(data):
            data['samples'] = [sample('x%d' % i, '0', '1', [])
                               for i in range(9)]
        def too_many_measurements(data):
            data['samples'][0]['measurements'] = [
                measurement('m%d' % i, 'H', '0', '1') for i in range(121)]
        def empty_species(data):
            data['species'] = []
        def empty_samples(data):
            data['samples'] = []
        for mutate in (too_many_species, too_many_samples,
                       too_many_measurements, empty_species, empty_samples):
            with self.subTest(mutate=mutate.__name__):
                self.assert_rejected(mutate)


class CliRegression(unittest.TestCase):
    def run_cli(self, source, output):
        return subprocess.run(
            [sys.executable, '-m', 'reaction_batch.assay_cli',
             '--source', str(source), '--output', str(output)],
            capture_output=True, cwd=ROOT,
            env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')})

    def test_fixtures_and_sha256(self):
        for name in ('assay-example1.json', 'assay-example2.json',
                     'assay-example3.json', 'assay-mixed-48.json'):
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as tmp:
                    output = Path(tmp) / 'report.json'
                    source = ROOT / 'fixtures' / name
                    completed = self.run_cli(source, output)
                    self.assertEqual(completed.returncode, 0,
                                     completed.stderr[:400])
                    report = json.loads(output.read_text())
                    import hashlib
                    self.assertEqual(report['source_sha256'],
                                     hashlib.sha256(source.read_bytes()).hexdigest())
                    data = json.loads(source.read_text())
                    check_report(self, data, report)
                    self.assertEqual(report['schema_version'], 1)

    def test_rerun_is_byte_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            first = Path(tmp) / 'a.json'
            second = Path(tmp) / 'b.json'
            source = ROOT / 'fixtures' / 'assay-mixed-48.json'
            self.assertEqual(self.run_cli(source, first).returncode, 0)
            self.assertEqual(self.run_cli(source, second).returncode, 0)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_invalid_input_preserves_old_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'bad.json'
            output = Path(tmp) / 'out.json'
            source.write_text('{"species": [], "samples": []}')
            output.write_text('old report')
            completed = self.run_cli(source, output)
            self.assertEqual(completed.returncode, 2)
            self.assertEqual(output.read_text(), 'old report')

    def test_output_guards(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = ROOT / 'fixtures' / 'assay-example1.json'
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
                '{"species": [], "species": [], "samples": []}')
            output = Path(tmp) / 'out.json'
            completed = self.run_cli(source, output)
            self.assertEqual(completed.returncode, 2)
            self.assertFalse(output.exists())


class ScaleAndMixedBatch(unittest.TestCase):
    def test_forty_eight_species_multiple_degrees_of_freedom(self):
        data = json.loads((ROOT / 'fixtures' / 'assay-mixed-48.json').read_text())
        report = analyze(data)
        check_report(self, data, {'result': report})
        self.assertEqual(len(report['species_order']), 48)
        entries = report['samples'][0]['ranges']
        self.assertTrue(any(entry['minimum']['value'] != entry['maximum']['value']
                            for entry in entries))
        # The marker species contains the only Se atoms, so its amount is
        # pinned to the Se window and has a strictly positive lower bound.
        marker = next(entry for entry in entries
                      if entry['species_id'] == 's47')
        self.assertNotEqual(marker['minimum']['value'], '0')
        self.assertLess(Fraction(marker['minimum']['value']),
                        Fraction(marker['maximum']['value']))

    def test_mixed_batch_reports_each_sample(self):
        data = json.loads((ROOT / 'fixtures' / 'assay-mixed-48.json').read_text())
        report = analyze(data)
        statuses = [(item['sample_id'], item['status'])
                    for item in report['samples']]
        self.assertIn(('contradictory', 'inconsistent'), statuses)
        self.assertIn(('zero-nonempty', 'inconsistent'), statuses)
        self.assertIn(('forty-eight-dof', 'consistent'), statuses)
        self.assertIn(('empty-zero', 'consistent'), statuses)

    def test_reordering_input_rebuilds_certificates(self):
        base = {'species': [species_pair('h', 'H2'), species_pair('o', 'O2'),
                            species_pair('w', 'H2O')],
                'samples': [sample('s', '1', '3', [
                    measurement('m1', 'H', '2', '4'),
                    measurement('m2', 'O', '1', '2'),
                    measurement('m3', 'H', '3', '3')])]}
        report = analyze(base)
        check_report(self, base, {'result': report})
        reordered = {
            'species': [base['species'][2], base['species'][0],
                        base['species'][1]],
            'samples': [sample('s', '1', '3', [
                measurement('m3', 'H', '3', '3'),
                measurement('m1', 'H', '2', '4'),
                measurement('m2', 'O', '1', '2')])],
        }
        second = analyze(reordered)
        check_report(self, reordered, {'result': second})
        self.assertEqual(second['species_order'], ['w', 'h', 'o'])
        by_id = {entry['species_id']: entry for entry in report['samples'][0]['ranges']}
        for entry in second['samples'][0]['ranges']:
            original = by_id[entry['species_id']]
            self.assertEqual(entry['minimum']['value'], original['minimum']['value'])
            self.assertEqual(entry['maximum']['value'], original['maximum']['value'])

    def test_maximum_scale_within_seconds(self):
        import random, time
        random.seed(808)
        elements = ['H', 'C', 'N', 'O', 'S', 'P', 'Cl', 'Na', 'Fe', 'K',
                    'Mg', 'Ca']
        texts = []
        seen = set()
        while len(texts) < 60:
            chosen = random.sample(elements, random.randint(2, 4))
            text = ''.join(element + (str(random.randint(2, 4))
                                      if random.random() < 0.5 else '')
                           for element in chosen)
            if text not in seen:
                seen.add(text)
                texts.append(text)
        species = [species_pair('x%d' % i, text)
                   for i, text in enumerate(texts)]
        measurements = []
        for q in range(120):
            low = random.randint(0, 100)
            measurements.append(measurement(
                'm%d' % q, random.choice(elements),
                str(low), str(random.randint(low, 300))))
        data = {'species': species,
                'samples': [sample('big', '0', '80', measurements)]}
        started = time.monotonic()
        report = analyze(data)
        elapsed = time.monotonic() - started
        self.assertLess(elapsed, 60)
        check_report(self, data, {'result': report})


if __name__ == '__main__':
    unittest.main()
