import json
import os
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from reaction_batch.chemistry import balance, formula
from reaction_batch.lot_model import analyze

ROOT = Path(__file__).resolve().parents[1]


def stock(op_id, lot, formula_text, moles):
    return {'op_id': op_id, 'type': 'stock', 'lot': lot, 'formula': formula_text, 'moles': moles}


def run(op_id, inputs, products):
    return {'op_id': op_id, 'type': 'run', 'inputs': inputs, 'products': products}


def hold(op_id, lot):
    return {'op_id': op_id, 'type': 'hold', 'lot': lot}


def release(op_id, lot):
    return {'op_id': op_id, 'type': 'release', 'lot': lot}


def inp(lot, moles):
    return {'lot': lot, 'moles': moles}


def prod(lot, formula_text):
    return {'lot': lot, 'formula': formula_text}


def report(operations):
    return analyze({'operations': operations})


def lots_by_id(result):
    return {item['lot']: item for item in result['lots']}


def statuses(result):
    return {item['op_id']: item['status'] for item in result['operations']}


class OrderingAndFractions(unittest.TestCase):
    def test_documented_example(self):
        result = report([
            stock('s1', 'h1', 'H2', '0.1'),
            stock('s2', 'h2', 'H2', '0.2'),
            stock('s3', 'o', 'O2', '0.15'),
            run('r1', [inp('h1', '0.1'), inp('h2', '0.2'), inp('o', '0.15')], [prod('w', 'H2O')]),
            hold('h1op', 'h1'),
        ])
        lots = lots_by_id(result)
        self.assertEqual(lots['w']['remaining_moles'], '3/10')
        self.assertTrue(lots['w']['blocked'])
        self.assertTrue(lots['h1']['blocked'])
        self.assertFalse(lots['h2']['blocked'])
        for lot in ('h1', 'h2', 'o'):
            self.assertEqual(lots[lot]['remaining_moles'], '0')
        self.assertEqual(result['active_holds'], ['h1'])

    def test_species_input_order_does_not_change_result(self):
        base = [
            stock('s1', 'h1', 'H2', '3'),
            stock('s2', 'o', 'O2', '3'),
            stock('s3', 'n', 'N2', '3'),
        ]
        forward = base + [run('r', [inp('h1', '3'), inp('o', '1.5'), inp('n', '1')], [prod('w', 'H2O'), prod('a', 'NH3')])]
        backward = base + [run('r', [inp('n', '1'), inp('o', '1.5'), inp('h1', '3')], [prod('a', 'NH3'), prod('w', 'H2O')])]
        self.assertEqual(report(forward), report(backward))

    def test_exact_fraction_and_atom_conservation(self):
        operations = [
            stock('s1', 'a', 'Al', '10'),
            stock('s2', 'd', 'HCl', '21'),
            run('r', [inp('a', '3'), inp('d', '9')], [prod('salt', 'AlCl3'), prod('gas', 'H2')]),
        ]
        result = report(operations)
        lots = lots_by_id(result)
        consumed = {'Al': Fraction(3), 'HCl': Fraction(9)}
        produced = {'AlCl3': Fraction(lots['salt']['remaining_moles']), 'H2': Fraction(lots['gas']['remaining_moles'])}
        for element in ('Al', 'H', 'Cl'):
            before = sum(formula(name).get(element, 0) * amount for name, amount in consumed.items())
            after = sum(formula(name).get(element, 0) * amount for name, amount in produced.items())
            self.assertEqual(before, after)

    def test_same_formula_parents_are_merged_not_aliased(self):
        operations = [
            stock('s1', 'h1', 'H2', '1'),
            stock('s2', 'h2', 'H2', '2'),
            stock('s3', 'o', 'O2', '2'),
            run('r', [inp('h1', '1'), inp('h2', '2'), inp('o', '1.5')], [prod('w', 'H2O')]),
        ]
        lots = lots_by_id(report(operations))
        self.assertEqual(lots['w']['remaining_moles'], '3')
        self.assertEqual(lots['h1']['remaining_moles'], '0')
        self.assertEqual(lots['h2']['remaining_moles'], '0')
        self.assertEqual(lots['o']['remaining_moles'], '1/2')

    def test_ratio_mismatch_rolls_back_everything(self):
        operations = [
            stock('s1', 'h', 'H2', '2'),
            stock('s2', 'o', 'O2', '2'),
            run('r', [inp('h', '2'), inp('o', '2')], [prod('w', 'H2O')]),
        ]
        result = report(operations)
        lots = lots_by_id(result)
        self.assertEqual(statuses(result)['r'], 'rejected')
        self.assertNotIn('w', lots)
        self.assertEqual(lots['h']['remaining_moles'], '2')
        self.assertEqual(lots['o']['remaining_moles'], '2')

    def test_unbalanceable_is_business_rejection(self):
        operations = [
            stock('s1', 'c', 'C', '1'),
            stock('s2', 'o', 'O2', '1'),
            run('r', [inp('c', '1'), inp('o', '1')], [prod('co', 'CO'), prod('co2', 'CO2')]),
        ]
        result = report(operations)
        self.assertEqual(statuses(result)['r'], 'rejected')
        self.assertEqual(len(result['lots']), 2)

    def test_duplicate_or_reactant_product_formula_is_rejected(self):
        operations = [
            stock('s1', 'a', 'A', '2'),
            stock('s2', 'b', 'B', '2'),
            run('dup_formula', [inp('a', '1'), inp('b', '1')], [prod('p1', 'AB'), prod('p2', 'AB')]),
            run('reactant_formula', [inp('a', '1'), inp('b', '1')], [prod('p3', 'A')]),
        ]
        result = report(operations)
        self.assertEqual(statuses(result)['dup_formula'], 'rejected')
        self.assertEqual(statuses(result)['reactant_formula'], 'rejected')
        self.assertEqual(lots_by_id(result)['a']['remaining_moles'], '2')

    def test_failed_run_leaves_no_lineage_edge(self):
        operations = [
            stock('s1', 'a', 'A', '2'),
            stock('s2', 'b', 'B', '1'),
            run('failed', [inp('a', '1'), inp('b', '2')], [prod('P', 'AB')]),
            hold('h', 'a'),
            stock('s3', 'P', 'AB', '1'),
        ]
        result = report(operations)
        lots = lots_by_id(result)
        self.assertEqual(statuses(result)['failed'], 'rejected')
        self.assertFalse(lots['P']['blocked'])
        self.assertTrue(lots['a']['blocked'])

    def test_failed_run_does_not_reserve_lot_names(self):
        operations = [
            stock('s1', 'h', 'H2', '2'),
            stock('s2', 'o', 'O2', '2'),
            run('bad', [inp('h', '2')], [prod('w', 'H2O')]),
            run('good', [inp('h', '2'), inp('o', '1')], [prod('w', 'H2O')]),
        ]
        result = report(operations)
        lots = lots_by_id(result)
        self.assertEqual(statuses(result)['bad'], 'rejected')
        self.assertEqual(statuses(result)['good'], 'accepted')
        self.assertEqual(lots['w']['remaining_moles'], '2')
        self.assertFalse(lots['w']['blocked'])


class Isolation(unittest.TestCase):
    def diamond(self):
        return [
            stock('s', 'S', 'A2', '8'),
            run('b', [inp('S', '1')], [prod('B', 'A')]),
            run('c', [inp('S', '1')], [prod('C', 'A')]),
            run('d', [inp('B', '2'), inp('C', '2')], [prod('D', 'A2')]),
        ]

    def test_one_of_two_diamond_roots_still_blocks(self):
        operations = self.diamond() + [hold('hS', 'S'), hold('hB', 'B'), release('rS', 'S')]
        result = report(operations)
        lots = lots_by_id(result)
        self.assertNotIn('S', result['active_holds'])
        self.assertEqual(result['active_holds'], ['B'])
        self.assertFalse(lots['S']['blocked'])
        self.assertFalse(lots['C']['blocked'])
        self.assertTrue(lots['B']['blocked'])
        self.assertTrue(lots['D']['blocked'])

    def test_release_only_when_all_roots_release(self):
        operations = self.diamond() + [hold('hS', 'S'), hold('hB', 'B'), release('rS', 'S'), release('rB', 'B')]
        lots = lots_by_id(report(operations))
        self.assertEqual([item['blocked'] for item in lots.values()], [False] * 4)

    def test_depleted_ancestor_hold_blocks_existing_product(self):
        operations = [
            stock('s', 'S', 'A', '2'),
            run('r', [inp('S', '2')], [prod('P', 'A2')]),
            hold('h', 'S'),
        ]
        lots = lots_by_id(report(operations))
        self.assertEqual(lots['S']['remaining_moles'], '0')
        self.assertTrue(lots['P']['blocked'])

    def test_blocked_input_rejects_without_changes(self):
        operations = [
            stock('s', 'S', 'A', '4'),
            hold('h', 'S'),
            run('r', [inp('S', '1')], [prod('P', 'A')]),
        ]
        result = report(operations)
        lots = lots_by_id(result)
        self.assertEqual(statuses(result)['r'], 'rejected')
        self.assertEqual(lots['S']['remaining_moles'], '4')
        self.assertNotIn('P', lots)

    def test_cannot_feed_until_all_covering_roots_release(self):
        operations = self.diamond() + [
            hold('hS', 'S'), hold('hB', 'B'), release('rS', 'S'),
            run('blocked', [inp('D', '1')], [prod('E', 'A')]),
            release('rB', 'B'),
            run('allowed', [inp('D', '1')], [prod('E', 'A')]),
        ]
        result = report(operations)
        self.assertEqual(statuses(result)['blocked'], 'rejected')
        self.assertEqual(statuses(result)['allowed'], 'accepted')
        self.assertEqual(lots_by_id(result)['E']['remaining_moles'], '2')

    def test_failed_run_product_name_is_not_affected_by_ancestor_hold(self):
        operations = [
            stock('s', 'S', 'A', '1'),
            hold('h', 'S'),
            run('r', [inp('S', '1')], [prod('P', 'A')]),
        ]
        result = report(operations)
        self.assertNotIn('P', lots_by_id(result))
        self.assertFalse(any(item['blocked'] for item in result['lots'] if item['lot'] == 'P'))

    def test_hold_already_blocked_lot_is_accepted_root(self):
        operations = self.diamond() + [hold('hS', 'S'), hold('hB', 'B')]
        result = report(operations)
        self.assertEqual(statuses(result)['hB'], 'accepted')
        self.assertEqual(result['active_holds'], ['B', 'S'])

    def test_duplicate_active_hold_rejected_and_nonroot_release_rejected(self):
        operations = [stock('s', 'S', 'A', '1'), hold('h1', 'S'), hold('h2', 'S'), release('r1', 'S'), release('r2', 'S')]
        result = report(operations)
        self.assertEqual(statuses(result), {'s': 'accepted', 'h1': 'accepted', 'h2': 'rejected', 'r1': 'accepted', 'r2': 'rejected'})
        self.assertEqual(result['active_holds'], [])

    def test_rejected_hold_release_still_count_towards_limit(self):
        operations = [hold('h1', 'missing') for _ in range(501)]
        operations = [dict(item, op_id=item['op_id'] + str(index)) for index, item in enumerate(operations)]
        with self.assertRaises(ValueError):
            report(operations)


class InputContract(unittest.TestCase):
    def test_unknown_lot_is_rejection(self):
        result = report([run('r', [inp('ghost', '1')], [prod('p', 'A')])])
        self.assertEqual(statuses(result)['r'], 'rejected')
        self.assertEqual(result['lots'], [])

    def test_unknown_lot_rejects_whole_run_without_partial_deduction(self):
        operations = [
            stock('s1', 'a', 'A', '5'),
            stock('s2', 'b', 'B', '5'),
            run('r', [inp('a', '1'), inp('ghost', '1'), inp('b', '2')], [prod('p', 'AB')]),
        ]
        result = report(operations)
        lots = lots_by_id(result)
        self.assertEqual(statuses(result)['r'], 'rejected')
        self.assertEqual((lots['a']['remaining_moles'], lots['b']['remaining_moles']), ('5', '5'))
        self.assertNotIn('p', lots)

    def test_maximum_inputs_and_products_are_accepted(self):
        # 40 distinct same-formula lots merge to one reactant species.
        operations = [stock(f's{index}', f'l{index}', 'A2', '1') for index in range(40)]
        operations.append(run('r', [inp(f'l{index}', '1') for index in range(40)], [prod('p', 'A')]))
        result = report(operations)
        self.assertEqual(statuses(result)['r'], 'accepted')
        self.assertEqual(lots_by_id(result)['p']['remaining_moles'], '80')

        # One compound of 20 elements balanced into exactly 20 products.
        elements = [chr(ord('A') + index) for index in range(20)]
        operations = [
            stock('s', 'whole', ''.join(elements), '3'),
            run('split', [inp('whole', '3')], [prod(f'p{index}', element) for index, element in enumerate(elements)]),
        ]
        result = report(operations)
        self.assertEqual(statuses(result)['split'], 'accepted')
        self.assertEqual(len(result['lots']), 21)
        self.assertTrue(all(lots_by_id(result)[f'p{index}']['remaining_moles'] == '3' for index in range(20)))

    def test_hold_release_limit_counts_rejected_records(self):
        operations = [hold(f'h{index}', 'missing') for index in range(500)]
        for index, operation in enumerate(operations):
            operation['op_id'] = f'h{index}'
        self.assertEqual(statuses(report(operations)), {f'h{index}': 'rejected' for index in range(500)})
        operations.append(release('extra', 'missing'))
        with self.assertRaises(ValueError):
            report(operations)

    def test_duplicate_op_id_is_input_error(self):
        with self.assertRaises(ValueError):
            report([stock('same', 'a', 'A', '1'), stock('same', 'b', 'B', '1')])

    def test_stock_identity_reuse_is_rejection(self):
        result = report([stock('a', 'x', 'A', '1'), stock('b', 'x', 'B', '1')])
        self.assertEqual(statuses(result), {'a': 'accepted', 'b': 'rejected'})
        self.assertEqual(lots_by_id(result)['x']['formula'], 'A')

    def test_run_cannot_reuse_existing_product_identity(self):
        operations = [
            stock('s1', 'x', 'A', '1'),
            stock('s2', 'y', 'A', '1'),
            run('r', [inp('x', '1')], [prod('y', 'A2')]),
        ]
        result = report(operations)
        self.assertEqual(statuses(result)['r'], 'rejected')
        self.assertEqual(lots_by_id(result)['y']['remaining_moles'], '1')
        self.assertEqual(lots_by_id(result)['x']['remaining_moles'], '1')

    def test_empty_operations(self):
        self.assertEqual(report([]), {'operations': [], 'lots': [], 'active_holds': []})

    def test_structurally_invalid_inputs_raise(self):
        base = [stock('s1', 'x', 'A', '1')]
        bad = [
            [{'op_id': 'x', 'type': 'stock', 'lot': 'y', 'formula': 'A'}],
            [dict(stock('s2', 'y', 'A', '1'), extra=1)],
            [stock('s1', 'y', 'A', '1')],
            [stock('s1', 'x', 'A', '0')],
            [stock('s1', 'x', 'A', True)],
            [run('r', [], [prod('p', 'A')])],
            [run('r', [inp('x', '1'), inp('x', '2')], [prod('p', 'A')])],
            [run('r', [inp('x', '1')], [prod('p', 'A'), prod('p', 'A')])],
            [run('r', [inp('x', '1')] * 41, [prod('p', 'A')])],
            [run('r', [inp('x', '1')], [prod(f'p{i}', 'A') for i in range(21)])],
            [{'op_id': 'x', 'type': 'hold', 'lot': 'y', 'moles': '1'}],
        ]
        for operations in bad:
            with self.subTest(operations=operations), self.assertRaises(ValueError):
                report(base + operations)

    def test_merged_species_ceiling_is_input_error(self):
        elements = [chr(ord('A') + index) for index in range(21)]
        operations = [stock(f's{index}', f'l{index}', element, '1') for index, element in enumerate(elements)]
        operations.append(run('r', [inp(f'l{index}', '1') for index in range(21)], [prod('p', 'A2')]))
        with self.assertRaises(ValueError):
            report(operations)


class CliIntegration(unittest.TestCase):
    def invoke(self, source, output):
        return subprocess.run(
            [sys.executable, '-m', 'reaction_batch.lots', '--source', str(source), '--output', str(output)],
            capture_output=True,
            env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')},
        )

    def test_lots_cli_writes_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'in.json'
            output = root / 'out.json'
            source.write_text(json.dumps({'operations': [stock('s', 'x', 'A', '1')]}))
            completed = self.invoke(source, output)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(json.loads(output.read_text())['schema_version'], 1)

    def test_batch_cli_still_accepts_old_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            completed = subprocess.run(
                [sys.executable, '-m', 'reaction_batch.cli', '--source', str(ROOT / 'fixtures/input.json'), '--output', str(Path(tmp) / 'batch.json')],
                capture_output=True,
                cwd=ROOT,
                env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')},
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(json.loads((Path(tmp) / 'batch.json').read_text())['result']['coefficients'], {'H2': 2, 'O2': 1, 'H2O': 2})

    def test_output_symlink_and_same_path_are_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'in.json'
            source.write_text(json.dumps({'operations': [stock('s', 'x', 'A', '1')]}))
            previous = source.read_bytes()
            link = root / 'link'
            link.symlink_to(source)
            for target in (source, link):
                completed = self.invoke(source, target)
                self.assertEqual(completed.returncode, 2)
                self.assertEqual(source.read_bytes(), previous)

    def test_same_bytes_same_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'in.json'
            source.write_text(json.dumps({'operations': [stock('s', 'x', 'A', '1')]}))
            first = root / 'a.json'
            second = root / 'b.json'
            self.assertEqual(self.invoke(source, first).returncode, 0)
            self.assertEqual(self.invoke(source, second).returncode, 0)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_invalid_batch_input_keeps_old_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'in.json'
            output = root / 'out.json'
            source.write_text('{"operations": [{"op_id": "x", "type": "nope"}]}')
            output.write_text('old')
            completed = self.invoke(source, output)
            self.assertEqual(completed.returncode, 2)
            self.assertEqual(output.read_text(), 'old')
