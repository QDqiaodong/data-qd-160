import json
import os
import resource
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

from reaction_batch.lot_model import analyze

ROOT = Path(__file__).resolve().parents[1]


def chain_input(length):
    """One stock then ``length - 1`` single-parent runs alternating A/A2."""
    operations = [{'op_id': 'o0', 'type': 'stock', 'lot': 'L0', 'formula': 'A2', 'moles': '2'}]
    for index in range(1, length):
        parent, child = f'L{index - 1}', f'L{index}'
        if index % 2 == 1:
            operations.append({'op_id': f'o{index}', 'type': 'run', 'inputs': [{'lot': parent, 'moles': '1'}], 'products': [{'lot': child, 'formula': 'A'}]})
        else:
            operations.append({'op_id': f'o{index}', 'type': 'run', 'inputs': [{'lot': parent, 'moles': '2'}], 'products': [{'lot': child, 'formula': 'A2'}]})
    return operations


def convergence_input(runs):
    """Two parents per run and two products per run, within all limits."""
    operations = [
        {'op_id': 's0', 'type': 'stock', 'lot': 'L0', 'formula': 'CH4', 'moles': '1'},
        {'op_id': 's1', 'type': 'stock', 'lot': 'L1', 'formula': 'O2', 'moles': '2'},
    ]
    previous = ('L0', 'L1')
    for index in range(runs):
        left, right = previous
        if index % 2 == 0:
            products = ('CO2', 'H2O')
        else:
            products = ('CH4', 'O2')
        operations.append({
            'op_id': f'r{index}',
            'type': 'run',
            'inputs': [{'lot': left, 'moles': '1'}, {'lot': right, 'moles': '2'}],
            'products': [{'lot': f'P{index}a', 'formula': products[0]}, {'lot': f'P{index}b', 'formula': products[1]}],
        })
        previous = (f'P{index}a', f'P{index}b')
    return operations


MEMORY_PROBE = """
import json, resource, sys
from reaction_batch.lot_model import analyze
analyze(json.load(open(sys.argv[1])))
print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
"""


class Scale(unittest.TestCase):
    def test_20000_lot_chain_walks_hold_state_without_recursion(self):
        operations = chain_input(20000)
        started = time.monotonic()

        held = analyze({'operations': operations + [
            {'op_id': 'h', 'type': 'hold', 'lot': 'L0'},
        ]})
        self.assertEqual(len(held['lots']), 20000)
        self.assertTrue(all(item['blocked'] for item in held['lots']))
        self.assertEqual(held['lots'], sorted(held['lots'], key=lambda item: item['lot']))
        self.assertEqual({item['lot'] for item in held['lots']}, {f'L{index}' for index in range(20000)})

        released = analyze({'operations': operations + [
            {'op_id': 'h', 'type': 'hold', 'lot': 'L0'},
            {'op_id': 'r', 'type': 'release', 'lot': 'L0'},
        ]})
        self.assertFalse(any(item['blocked'] for item in released['lots']))
        self.assertEqual(released['active_holds'], [])
        self.assertLess(time.monotonic() - started, 30)

    def test_tail_lot_toggles_between_hold_and_release(self):
        operations = chain_input(2000)
        by_id = {item['lot']: item for item in analyze({'operations': operations + [
            {'op_id': 'h', 'type': 'hold', 'lot': 'L500'},
        ]})['lots']}
        self.assertTrue(by_id['L500']['blocked'])
        self.assertTrue(by_id['L1999']['blocked'])
        self.assertFalse(by_id['L499']['blocked'])
        self.assertFalse(by_id['L0']['blocked'])
        after = analyze({'operations': operations + [
            {'op_id': 'h', 'type': 'hold', 'lot': 'L500'},
            {'op_id': 'r', 'type': 'release', 'lot': 'L500'},
        ]})
        self.assertFalse(any(item['blocked'] for item in after['lots']))

    def test_convergence_graph_cli_within_memory_budget(self):
        operations = convergence_input(14000) + [{'op_id': 'h', 'type': 'hold', 'lot': 'L0'}]
        self.assertEqual(len(operations), 14003)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'in.json'
            output = root / 'out.json'
            source.write_text(json.dumps({'operations': operations}))
            completed = subprocess.run(
                [sys.executable, '-m', 'reaction_batch.lots', '--source', str(source), '--output', str(output)],
                capture_output=True,
                cwd=ROOT,
                env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')},
            )
            self.assertEqual(completed.returncode, 0, completed.stderr[:500])
            report = json.loads(output.read_text())
            self.assertEqual(len(report['result']['lots']), 28002)
            self.assertEqual(report['result']['active_holds'], ['L0'])
            self.assertTrue(all(item['blocked'] for item in report['result']['lots'] if item['lot'].startswith('P')))

            probe = subprocess.run(
                [sys.executable, '-c', MEMORY_PROBE, str(source)],
                capture_output=True,
                cwd=ROOT,
                env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')},
            )
            self.assertEqual(probe.returncode, 0, probe.stderr[:500])
            peak_mib = int(probe.stdout.strip()) / 1024
            self.assertLess(peak_mib, 192)


if __name__ == '__main__':
    unittest.main()
