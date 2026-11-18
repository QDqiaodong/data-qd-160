"""Leakage-safe split tests: independent oracle, regressions and scale.

Small inputs are checked against an oracle that builds the leakage components
by *directly testing every record pair* for equal element counts or equal
experiment id, traverses them with a plain breadth-first search, and then
enumerates every whole-component assignment.  The oracle never uses the
production disjoint-set structure, the production suffix reachability table or
the production greedy reconstruction; it compares the complete 0/1 record-id
vectors directly to resolve ties.  Scale cases use hand-derived expectations.
"""

import itertools
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from reaction_batch.chemistry import formula
from reaction_batch.split import analyze

ROOT = Path(__file__).resolve().parents[1]


def record(record_id, formula_text, experiment_id, label, pin='free'):
    return {'record_id': record_id, 'formula': formula_text,
            'experiment_id': experiment_id, 'label': label, 'pin': pin}


def quota(lower, upper, target):
    return {'lower': lower, 'upper': upper, 'target': target}


def quotas(positive, negative):
    return {'positive': positive, 'negative': negative}


# ---------------------------------------------------------------------------
# Independent oracle


def oracle_groups(data):
    """Components via pairwise relation checks + BFS, independent of the DSU."""
    records = data['records']
    counts = [tuple(sorted(formula(entry['formula']).items()))
              for entry in records]
    size = len(records)
    adjacency = [set() for _ in range(size)]
    for left in range(size):
        for right in range(left + 1, size):
            if (counts[left] == counts[right]
                    or records[left]['experiment_id']
                    == records[right]['experiment_id']):
                adjacency[left].add(right)
                adjacency[right].add(left)
    seen = set()
    components = []
    for start in range(size):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        members = []
        while stack:
            current = stack.pop()
            members.append(current)
            for neighbour in adjacency[current]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        components.append(members)
    result = []
    for members in components:
        identifiers = sorted(records[index]['record_id'] for index in members)
        pins = {records[index]['pin'] for index in members}
        if 'train' in pins and 'test' in pins:
            forced = 'conflict'
        elif 'train' in pins:
            forced = 'train'
        elif 'test' in pins:
            forced = 'test'
        else:
            forced = 'free'
        result.append({
            'group_id': identifiers[0],
            'record_ids': identifiers,
            'counts': {'positive': sum(1 for index in members
                                       if records[index]['label'] == 'positive'),
                       'negative': sum(1 for index in members
                                       if records[index]['label'] == 'negative')},
            'forced': forced,
            'members': members,
        })
    result.sort(key=lambda entry: entry['group_id'])
    return result


def oracle(data):
    records = data['records']
    groups = oracle_groups(data)
    plain = [{key: group[key] for key in
              ('group_id', 'record_ids', 'counts', 'forced')}
             for group in groups]
    conflicts = [group['group_id'] for group in groups
                 if group['forced'] == 'conflict']
    if conflicts:
        return {'groups': plain, 'status': 'infeasible',
                'reason': 'pin_conflict', 'conflict_groups': conflicts}

    positive_quota = data['test_quotas']['positive']
    negative_quota = data['test_quotas']['negative']
    order = sorted(range(len(records)),
                   key=lambda index: records[index]['record_id'])
    best = None
    for choice in itertools.product((0, 1), repeat=len(groups)):
        sides = {}
        for selected, group in zip(choice, groups):
            for index in group['members']:
                sides[index] = selected
        valid = True
        for index, entry in enumerate(records):
            if entry['pin'] == 'train' and sides[index] == 1:
                valid = False
            if entry['pin'] == 'test' and sides[index] == 0:
                valid = False
        if not valid:
            continue
        test_positive = sum(1 for index in order
                            if sides[index] and records[index]['label'] == 'positive')
        test_negative = sum(1 for index in order
                            if sides[index] and records[index]['label'] == 'negative')
        if not positive_quota['lower'] <= test_positive <= positive_quota['upper']:
            continue
        if not negative_quota['lower'] <= test_negative <= negative_quota['upper']:
            continue
        deviation = (abs(test_positive - positive_quota['target'])
                     + abs(test_negative - negative_quota['target']))
        vector = tuple(sides[index] for index in order)
        key = (deviation, vector)
        if best is None or key < best[0]:
            best = (key, sides, test_positive, test_negative)

    if best is None:
        return {'groups': plain, 'status': 'infeasible', 'reason': 'quotas'}
    _key, sides, test_positive, test_negative = best
    train = [records[index]['record_id'] for index in order
             if not sides[index]]
    test = [records[index]['record_id'] for index in order
            if sides[index]]
    return {'groups': plain, 'status': 'feasible',
            'split': {'train': train, 'test': test,
                      'test_counts': {'positive': test_positive,
                                      'negative': test_negative},
                      'deviation': (abs(test_positive - positive_quota['target'])
                                    + abs(test_negative - negative_quota['target']))}}


# ---------------------------------------------------------------------------


class OracleComparison(unittest.TestCase):
    def check(self, data):
        expected = oracle(data)
        produced = analyze(data)
        self.assertEqual(produced, expected)
        return produced

    def test_bridge_free(self):
        data = {'records': [record('a', 'H2O', 'e1', 'positive'),
                            record('b', 'CO2', 'e2', 'negative')],
                'test_quotas': quotas(quota(1, 1, 1), quota(0, 0, 0))}
        report = self.check(data)
        self.assertEqual(report['split']['test'], ['a'])
        self.assertEqual(report['split']['train'], ['b'])

    def test_bridge_shared_experiment_quotas_infeasible(self):
        data = {'records': [record('a', 'H2O', 'e1', 'positive'),
                            record('b', 'CO2', 'e1', 'negative')],
                'test_quotas': quotas(quota(1, 1, 1), quota(0, 0, 0))}
        report = self.check(data)
        self.assertEqual(report['status'], 'infeasible')
        self.assertEqual(report['reason'], 'quotas')
        self.assertNotIn('split', report)

    def test_bridge_pins_become_conflict(self):
        data = {'records': [record('a', 'H2O', 'e1', 'positive', 'test'),
                            record('b', 'CO2', 'e1', 'negative', 'train')],
                'test_quotas': quotas(quota(1, 1, 1), quota(0, 0, 0))}
        report = self.check(data)
        self.assertEqual(report['reason'], 'pin_conflict')
        self.assertEqual(report['conflict_groups'], ['a'])

    def test_bridge_pins_separate_are_feasible(self):
        data = {'records': [record('a', 'H2O', 'e1', 'positive', 'test'),
                            record('b', 'CO2', 'e2', 'negative', 'train')],
                'test_quotas': quotas(quota(1, 1, 1), quota(0, 0, 0))}
        self.assertEqual(self.check(data)['status'], 'feasible')

    def test_pin_conflict_beats_quota_infeasibility(self):
        data = {'records': [record('a', 'A', 'e1', 'positive', 'test'),
                            record('b', 'B', 'e1', 'negative', 'train'),
                            record('c', 'C', 'e2', 'negative')],
                'test_quotas': quotas(quota(0, 0, 0), quota(0, 0, 0))}
        report = self.check(data)
        self.assertEqual(report['status'], 'infeasible')
        self.assertEqual(report['reason'], 'pin_conflict')
        self.assertEqual(report['conflict_groups'], ['a'])

    def test_forced_test_group_can_break_quota_upper(self):
        data = {'records': [record('a', 'A', 'e1', 'positive', 'test'),
                            record('b', 'B', 'e1', 'positive', 'free')],
                'test_quotas': quotas(quota(0, 1, 1), quota(0, 0, 0))}
        report = self.check(data)
        self.assertEqual(report['reason'], 'quotas')

    def test_separate_conflict_groups_sorted(self):
        data = {'records': [record('m', 'A', 'e1', 'positive', 'test'),
                            record('n', 'B', 'e1', 'negative', 'train'),
                            record('a', 'C', 'e2', 'positive', 'test'),
                            record('b', 'D', 'e2', 'negative', 'train')],
                'test_quotas': quotas(quota(0, 4, 1), quota(0, 4, 1))}
        report = self.check(data)
        self.assertEqual(report['conflict_groups'], ['a', 'm'])

    def test_transitive_closure_multi_hop(self):
        data = {'records': [record('a', 'H2O', 'e1', 'positive'),
                            record('b', 'OH2', 'e2', 'negative'),
                            record('c', 'CO2', 'e2', 'positive'),
                            record('d', 'O2', 'e3', 'negative')],
                'test_quotas': quotas(quota(1, 2, 1), quota(0, 2, 1))}
        report = self.check(data)
        self.assertEqual([group['record_ids'] for group in report['groups']],
                         [['a', 'b', 'c'], ['d']])
        self.assertEqual(report['groups'][0]['counts'],
                         {'positive': 2, 'negative': 1})

    def test_cross_count_tie_breaks_on_assignment(self):
        data = {'records': [record('a', 'A', 'e1', 'positive'),
                            record('b', 'A', 'e2', 'negative'),
                            record('c', 'A', 'e3', 'negative'),
                            record('d', 'B', 'e4', 'positive'),
                            record('e', 'B', 'e5', 'positive'),
                            record('f', 'B', 'e6', 'negative')],
                'test_quotas': quotas(quota(1, 2, 1), quota(1, 2, 1))}
        report = self.check(data)
        self.assertEqual(report['split']['test'], ['d', 'e', 'f'])
        self.assertEqual(report['split']['train'], ['a', 'b', 'c'])
        self.assertEqual(report['split']['deviation'], 1)

    def test_same_count_tie_breaks_on_assignment(self):
        data = {'records': [record('a', 'A', 'e1', 'positive'),
                            record('b', 'B', 'e2', 'positive'),
                            record('c', 'C', 'e3', 'negative'),
                            record('d', 'D', 'e4', 'negative')],
                'test_quotas': quotas(quota(1, 1, 1), quota(1, 1, 1))}
        report = self.check(data)
        self.assertEqual(report['split']['test'], ['b', 'd'])
        self.assertEqual(report['split']['train'], ['a', 'c'])

    def test_h4o2_does_not_merge_by_empirical_formula(self):
        data = {'records': [record('a', 'H2O', 'e1', 'positive'),
                            record('b', 'H4O2', 'e2', 'negative')],
                'test_quotas': quotas(quota(0, 1, 0), quota(0, 1, 0))}
        report = self.check(data)
        self.assertEqual(len(report['groups']), 2)

    def test_empty_records(self):
        data = {'records': [],
                'test_quotas': quotas(quota(0, 0, 3), quota(0, 0, 2))}
        report = self.check(data)
        self.assertEqual(report['split'],
                         {'train': [], 'test': [],
                          'test_counts': {'positive': 0, 'negative': 0},
                          'deviation': 5})

    def test_empty_records_positive_lower_is_infeasible(self):
        data = {'records': [],
                'test_quotas': quotas(quota(1, 2, 1), quota(0, 0, 0))}
        report = self.check(data)
        self.assertEqual(report['reason'], 'quotas')

    def test_all_pinned(self):
        data = {'records': [record('a', 'A', 'e1', 'positive', 'test'),
                            record('b', 'B', 'e2', 'negative', 'train'),
                            record('c', 'C', 'e3', 'positive', 'train'),
                            record('d', 'D', 'e4', 'negative', 'test')],
                'test_quotas': quotas(quota(1, 2, 1), quota(1, 2, 1))}
        self.check(data)

    def test_repeated_same_direction_pins(self):
        data = {'records': [record('a', 'A', 'e1', 'positive', 'test'),
                            record('b', 'B', 'e1', 'negative', 'test'),
                            record('c', 'C', 'e2', 'positive', 'test')],
                'test_quotas': quotas(quota(1, 3, 2), quota(0, 3, 1))}
        report = self.check(data)
        self.assertEqual(report['groups'][0]['forced'], 'test')

    def test_single_label_only(self):
        data = {'records': [record('a', 'A', 'e1', 'positive'),
                            record('b', 'B', 'e2', 'positive'),
                            record('c', 'C', 'e3', 'positive')],
                'test_quotas': quotas(quota(1, 2, 1), quota(0, 0, 0))}
        self.check(data)

    def test_quota_upper_above_sample_count(self):
        data = {'records': [record('a', 'A', 'e1', 'positive'),
                            record('b', 'B', 'e2', 'negative')],
                'test_quotas': quotas(quota(0, 300, 1), quota(0, 300, 1))}
        self.check(data)

    def test_target_outside_quota(self):
        data = {'records': [record('a', 'A', 'e1', 'positive'),
                            record('b', 'B', 'e2', 'negative')],
                'test_quotas': quotas(quota(0, 0, 5), quota(0, 0, 5))}
        report = self.check(data)
        self.assertEqual(report['split']['deviation'], 10)
        self.assertEqual(report['split']['test'], [])

    def test_count_hole_from_large_component(self):
        huge = [record(f'h{index}', 'H2O', 'shared', 'positive')
                for index in range(5)]
        data = {'records': huge + [record('z', 'CO2', 'e9', 'negative')],
                'test_quotas': quotas(quota(0, 300, 2), quota(0, 300, 1))}
        report = self.check(data)
        # The positive count can only be 0 or 5; target 2 picks 0 (D = 2+...)
        self.assertEqual(report['split']['test_counts']['positive'], 0)

    def test_non_ascii_record_ids(self):
        data = {'records': [record('ä', 'A', 'e1', 'positive'),
                            record('Z', 'B', 'e2', 'positive'),
                            record('é', 'C', 'e3', 'negative'),
                            record('a', 'D', 'e4', 'negative')],
                'test_quotas': quotas(quota(1, 2, 1), quota(1, 2, 1))}
        report = self.check(data)
        # Unicode order is Z, a, ä, é; train=0 sorts first, so both early
        # ids stay in train and the two later ids carry the quotas.
        self.assertEqual(report['split']['test'], ['ä', 'é'])

    def test_long_alternating_chain(self):
        records = []
        for index in range(12):
            if index % 2 == 0:
                records.append(record(f'r{index:02d}', 'H2O' if index % 4 == 0
                                      else 'OH2', f'e{index}', 'positive'))
            else:
                records.append(record(f'r{index:02d}', 'CO2' if index % 3 else 'O2C',
                                      f'e{index - 1}', 'negative'))
        data = {'records': records,
                'test_quotas': quotas(quota(0, 12, 3), quota(0, 12, 3))}
        report = self.check(data)
        self.assertEqual(len(report['groups']), 1)

    def test_random_small_inputs(self):
        rng = random.Random(20261118)
        pool = ['H2O', 'OH2', 'H4O2', 'CO2', 'O2C', 'CH4', 'A', 'A2',
                '(B2C)3', 'Zz9', 'C2H6O']
        for _case in range(180):
            size = rng.randrange(0, 9)
            records = []
            for index in range(size):
                records.append(record(
                    rng.choice(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h',
                                'ä', 'Ω', 'r1', 'r2']) + str(index),
                    rng.choice(pool),
                    'e' + str(rng.randrange(0, 3)),
                    rng.choice(['positive', 'negative']),
                    rng.choice(['free', 'free', 'free', 'train', 'test'])))
            data = {'records': records,
                    'test_quotas': quotas(
                        quota(rng.randrange(0, size + 2),
                              rng.randrange(0, size + 2),
                              rng.randrange(0, 6)),
                        quota(rng.randrange(0, size + 2),
                              rng.randrange(0, size + 2),
                              rng.randrange(0, 6)))}
            for label in ('positive', 'negative'):
                bound = data['test_quotas'][label]
                if bound['lower'] > bound['upper']:
                    bound['lower'], bound['upper'] = (bound['upper'],
                                                      bound['lower'])
            try:
                expected = oracle(data)
            except ValueError:
                continue
            produced = analyze(data)
            self.assertEqual(produced, expected,
                             msg=json.dumps(data, ensure_ascii=False))
            # Source array order must only affect nothing: reversing the
            # records leaves the business report identical.
            shuffled = {'records': list(reversed(data['records'])),
                        'test_quotas': json.loads(
                            json.dumps(data['test_quotas']))}
            self.assertEqual(analyze(shuffled), produced,
                             msg=json.dumps(shuffled, ensure_ascii=False))


# ---------------------------------------------------------------------------
# Scale


def scale_simple():
    records = [record(f'r{index:03d}', f'A{index + 1}', f'e{index:03d}',
                      'positive' if index < 300 else 'negative')
               for index in range(600)]
    return {'records': records,
            'test_quotas': quotas(quota(150, 150, 150),
                                  quota(150, 150, 150))}


def scale_mixed():
    """Mixed-label groups, both pin directions and several optimal cells.

    Free groups are 120 one-per-label pairs ``m000..m119`` followed by 60
    two-per-label pairs ``z000..z059``.  Forced-test singles contribute 30
    positive and 30 negative records; forced-train singles another 30 + 30.
    Targets 150/120 with quotas [100,200] make every ``P == N`` in [120,150]
    optimal (D = 30), so the canonical plan must still move to the
    lexicographically smallest record vector: all ``m`` groups train, the
    first fifteen ``z`` groups train, the rest test.
    """
    records = []
    counter = 0
    serial = 0

    def add(record_id, experiment_id, label, pin='free'):
        nonlocal serial
        records.append(record(record_id, f'A{serial + 1}', experiment_id,
                              label, pin))
        serial += 1

    for index in range(120):
        add(f'm{index:03d}n', f'em{index:03d}', 'negative')
        add(f'm{index:03d}p', f'em{index:03d}', 'positive')
    for index in range(60):
        for suffix in ('a', 'b'):
            add(f'z{index:03d}{suffix}0', f'ez{index:03d}', 'negative')
            add(f'z{index:03d}{suffix}1', f'ez{index:03d}', 'positive')
    for index in range(30):
        add(f'ftn{index:03d}', f'eftn{index:03d}', 'negative', 'test')
        add(f'ftp{index:03d}', f'eftp{index:03d}', 'positive', 'test')
    for index in range(30):
        add(f'frn{index:03d}', f'efrn{index:03d}', 'negative', 'train')
        add(f'frp{index:03d}', f'efrp{index:03d}', 'positive', 'train')
    assert len(records) == 600
    return {'records': records,
            'test_quotas': quotas(quota(100, 200, 150),
                                  quota(100, 200, 120))}


MEMORY_PROBE = """
import os, resource, subprocess, sys
env = {**os.environ, 'PYTHONPATH': sys.argv[1]}
completed = subprocess.run(
    [sys.executable, '-m', 'reaction_batch.split_cli', '--source', sys.argv[2],
     '--output', sys.argv[3]],
    capture_output=True, env=env)
sys.stderr.write(completed.stderr.decode()[:500])
sys.exit(completed.returncode or 0) if completed.returncode else None
print(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
"""


class Scale(unittest.TestCase):
    def test_600_independent_components(self):
        report = analyze(scale_simple())
        self.assertEqual(report['status'], 'feasible')
        self.assertEqual(len(report['groups']), 600)
        self.assertEqual(report['split']['test_counts'],
                         {'positive': 150, 'negative': 150})
        self.assertEqual(report['split']['deviation'], 0)
        self.assertEqual(report['split']['test'],
                         [f'r{index:03d}' for index in range(150, 300)]
                         + [f'r{index:03d}' for index in range(450, 600)])
        self.assertEqual(report['split']['train'],
                         [f'r{index:03d}' for index in range(150)]
                         + [f'r{index:03d}' for index in range(300, 450)])

    def test_600_mixed_groups_pins_and_multiple_optimal_cells(self):
        data = scale_mixed()
        report = analyze(data)
        self.assertEqual(report['status'], 'feasible')
        self.assertEqual(len(report['groups']), 120 + 60 + 60 + 60)
        split = report['split']
        self.assertEqual(split['test_counts'], {'positive': 120, 'negative': 120})
        self.assertEqual(split['deviation'], 30)
        expected_test = []
        for index in range(30):
            expected_test.append(f'ftn{index:03d}')
            expected_test.append(f'ftp{index:03d}')
        for index in range(15, 60):
            for suffix in ('a', 'b'):
                expected_test.append(f'z{index:03d}{suffix}0')
                expected_test.append(f'z{index:03d}{suffix}1')
        self.assertEqual(split['test'], sorted(expected_test))
        expected_train = []
        for index in range(120):
            expected_train.append(f'm{index:03d}n')
            expected_train.append(f'm{index:03d}p')
        for index in range(15):
            for suffix in ('a', 'b'):
                expected_train.append(f'z{index:03d}{suffix}0')
                expected_train.append(f'z{index:03d}{suffix}1')
        for index in range(30):
            expected_train.append(f'frn{index:03d}')
            expected_train.append(f'frp{index:03d}')
        self.assertEqual(split['train'], sorted(expected_train))
        # Every pin honoured and every group on one side.
        side = {record_id: 'test' for record_id in split['test']}
        side.update({record_id: 'train' for record_id in split['train']})
        for entry in data['records']:
            if entry['pin'] != 'free':
                self.assertEqual(side[entry['record_id']], entry['pin'])
        for group in report['groups']:
            self.assertEqual(len({side[record_id]
                                  for record_id in group['record_ids']}), 1)

    def test_cli_peak_rss_within_budget(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / 'out.json'
            for name, data in (('simple', scale_simple()),
                               ('mixed', scale_mixed())):
                source = root / (name + '.json')
                source.write_text(json.dumps(data))
                probe = subprocess.run(
                    [sys.executable, '-c', MEMORY_PROBE, str(ROOT / 'src'),
                     str(source), str(output)],
                    capture_output=True, cwd=ROOT)
                self.assertEqual(probe.returncode, 0, probe.stderr[:800])
                self.assertTrue(output.exists())
                peak_mib = int(probe.stdout.strip()) / 1024
                self.assertLess(peak_mib, 192,
                                f'{name} peak RSS {peak_mib:.1f} MiB')
                result = json.loads(output.read_text())['result']
                self.assertEqual(result['status'], 'feasible')


# ---------------------------------------------------------------------------
# Structural validation


def valid_sample():
    return {'records': [record('a', 'H2O', 'e1', 'positive'),
                        record('b', 'CO2', 'e2', 'negative')],
            'test_quotas': quotas(quota(1, 1, 1), quota(0, 0, 0))}


class Validation(unittest.TestCase):
    def rejects(self, data):
        with self.assertRaises(ValueError):
            analyze(data)

    def test_unknown_fields(self):
        data = valid_sample()
        data['extra'] = 1
        self.rejects(data)
        data = valid_sample()
        data['test_quotas']['positive']['extra'] = 1
        self.rejects(data)
        data = valid_sample()
        data['records'][0]['extra'] = 1
        self.rejects(data)

    def test_duplicate_record_id(self):
        data = valid_sample()
        data['records'][1]['record_id'] = 'a'
        self.rejects(data)

    def test_invalid_formula_pin_label(self):
        for field, value in (('formula', 'H2O)'), ('formula', ''),
                             ('formula', '2H'), ('pin', 'hold'),
                             ('label', 'Positive'), ('label', 'other')):
            data = valid_sample()
            data['records'][0][field] = value
            self.rejects(data)

    def test_boolean_quota_values(self):
        for field in ('lower', 'upper', 'target'):
            data = valid_sample()
            data['test_quotas']['positive'][field] = True
            self.rejects(data)

    def test_lower_above_upper(self):
        data = valid_sample()
        data['test_quotas']['negative'] = quota(2, 1, 1)
        self.rejects(data)

    def test_quota_out_of_range(self):
        data = valid_sample()
        data['test_quotas']['positive'] = quota(0, 301, 1)
        self.rejects(data)

    def test_record_and_label_limits(self):
        data = {'records': [record(f'r{index:04d}', f'A{index + 1}',
                                   f'e{index:04d}', 'positive')
                            for index in range(601)],
                'test_quotas': quotas(quota(0, 0, 0), quota(0, 0, 0))}
        self.rejects(data)
        data['records'] = data['records'][:600]
        self.rejects(data)  # 600 positives exceed the per-label limit of 300
        data['records'][300:] = [record(f'r{index:04d}', f'A{index + 1}',
                                        f'e{index:04d}', 'negative')
                                 for index in range(300, 600)]
        analyze(data)  # exactly 300 per label is allowed

    def test_records_must_be_array(self):
        data = valid_sample()
        data['records'] = {}
        self.rejects(data)

    def test_missing_top_level_fields(self):
        self.rejects({'records': []})
        self.rejects({'test_quotas': quotas(quota(0, 0, 0),
                                            quota(0, 0, 0))})


# ---------------------------------------------------------------------------
# CLI transport


class Transport(unittest.TestCase):
    def run_cli(self, source, output):
        env = {**os.environ, 'PYTHONPATH': str(ROOT / 'src')}
        return subprocess.run(
            [sys.executable, '-m', 'reaction_batch.split_cli',
             '--source', str(source), '--output', str(output)],
            capture_output=True, env=env, cwd=ROOT)

    def test_cli_success_and_determinism(self):
        data = valid_sample()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'in.json'
            output = root / 'out.json'
            source.write_text(json.dumps(data))
            completed = self.run_cli(source, output)
            self.assertEqual(completed.returncode, 0, completed.stderr[:500])
            first = json.loads(output.read_text())
            self.assertEqual(set(first), {'schema_version', 'source_sha256',
                                          'result'})
            self.assertEqual(first['schema_version'], 1)
            self.assertEqual(first['result'], oracle(data))
            self.assertEqual(len(first['source_sha256']), 64)
            completed = self.run_cli(source, output)
            self.assertEqual(completed.returncode, 0)
            self.assertEqual(json.loads(output.read_text()), first)

    def test_cli_reordered_records_keep_business_report(self):
        data = valid_sample()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'in.json'
            output = root / 'out.json'
            source.write_text(json.dumps(data))
            self.assertEqual(self.run_cli(source, output).returncode, 0)
            first = json.loads(output.read_text())
            data['records'].reverse()
            source.write_text(json.dumps(data))
            self.assertEqual(self.run_cli(source, output).returncode, 0)
            second = json.loads(output.read_text())
            self.assertEqual(first['result'], second['result'])
            self.assertNotEqual(first['source_sha256'], second['source_sha256'])

    def test_cli_input_errors_exit_2_and_keep_old_report(self):
        cases = []
        data = valid_sample()
        data['extra'] = 1
        cases.append(data)
        data = valid_sample()
        data['records'].append(record('a', 'H2O', 'e3', 'positive'))
        cases.append(data)
        data = valid_sample()
        data['records'][0]['formula'] = 'H2O('
        cases.append(data)
        data = valid_sample()
        data['records'][0]['pin'] = 'keep'
        cases.append(data)
        data = valid_sample()
        data['test_quotas']['positive']['target'] = False
        cases.append(data)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'in.json'
            output = root / 'out.json'
            output.write_text('old report')
            for data in cases:
                source.write_text(json.dumps(data))
                completed = self.run_cli(source, output)
                self.assertEqual(completed.returncode, 2)
                self.assertEqual(output.read_text(), 'old report')

    def test_cli_duplicate_json_and_nonfinite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / 'out.json'
            output.write_text('old')
            source = root / 'in.json'
            source.write_text('{"records": [], "records": []}')
            self.assertEqual(self.run_cli(source, output).returncode, 2)
            self.assertEqual(output.read_text(), 'old')
            source.write_text('{"records": [], "test_quotas": '
                              '{"positive": {"lower": 0, "upper": 0, '
                              '"target": NaN}, "negative": {"lower": 0, '
                              '"upper": 0, "target": 0}}}')
            self.assertEqual(self.run_cli(source, output).returncode, 2)
            self.assertEqual(output.read_text(), 'old')

    def test_cli_path_guards(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'in.json'
            source.write_text(json.dumps(valid_sample()))
            original = source.read_bytes()
            link = root / 'link.json'
            link.symlink_to(source)
            for destination in (source, link):
                completed = self.run_cli(source, destination)
                self.assertEqual(completed.returncode, 2)
                self.assertEqual(source.read_bytes(), original)

    def test_cli_matches_fixtures(self):
        expectations = {
            'split-example1.json': ('feasible', 1),
            'split-example2.json': ('feasible', 1),
            'split-example3.json': ('infeasible', None),
            'split-scale-600.json': ('feasible', 0),
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, (status, deviation) in expectations.items():
                source = ROOT / 'fixtures' / name
                output = root / (name + '.out')
                completed = self.run_cli(source, output)
                self.assertEqual(completed.returncode, 0, name)
                result = json.loads(output.read_text())['result']
                self.assertEqual(result['status'], status, name)
                if deviation is not None:
                    self.assertEqual(result['split']['deviation'], deviation)
                else:
                    self.assertEqual(result['reason'], 'pin_conflict')
                    self.assertEqual(result['conflict_groups'], ['a'])

    def test_cli_infeasible_is_a_successful_report(self):
        data = {'records': [record('a', 'H2O', 'e1', 'positive'),
                            record('b', 'CO2', 'e1', 'negative')],
                'test_quotas': quotas(quota(1, 1, 1), quota(0, 0, 0))}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'in.json'
            output = root / 'out.json'
            source.write_text(json.dumps(data))
            completed = self.run_cli(source, output)
            self.assertEqual(completed.returncode, 0)
            report = json.loads(output.read_text())
            self.assertEqual(report['result']['status'], 'infeasible')
            self.assertEqual(report['result']['reason'], 'quotas')


if __name__ == '__main__':
    unittest.main()
