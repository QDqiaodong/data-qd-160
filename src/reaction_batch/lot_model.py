"""Ordered lot accounting with exact fractions and lineage-based isolation."""
from functools import lru_cache

from .chemistry import balance, formula, quantity
from .common import ident, require, rows

_OPERATION_LIMIT = 30_000
_LOT_LIMIT = 30_000
_EDGE_LIMIT = 60_000
_HOLD_LIMIT = 500
_INPUT_LIMIT = 40
_PRODUCT_LIMIT = 20
_SPECIES_LIMIT = 20

@lru_cache(maxsize=8192)
def _parsed_formula(text):
    return formula(text)


@lru_cache(maxsize=4096)
def _cached_balance(reactants, products):
    return balance(list(reactants), list(products))


_TYPES = ('stock', 'run', 'hold', 'release')
_FIELDS = {
    'stock': {'op_id', 'lot', 'formula', 'moles', 'type'},
    'run': {'op_id', 'inputs', 'products', 'type'},
    'hold': {'op_id', 'lot', 'type'},
    'release': {'op_id', 'lot', 'type'},
}


def _positive_quantity(value):
    amount = quantity(value)
    if amount <= 0:
        raise ValueError('moles must be positive')
    return amount


def _object(value, fields):
    if not isinstance(value, dict) or set(value) != set(fields):
        raise ValueError('expected object fields: ' + ' '.join(sorted(fields)))
    return value


def _bounded_array(value, minimum, maximum):
    if not isinstance(value, list) or not minimum <= len(value) <= maximum:
        raise ValueError('array length outside allowed bounds')
    return value


def _validate_operations(data):
    """Collect every structural and literal error before business accounting."""
    require(data, 'operations')
    operations = rows(data['operations'], _OPERATION_LIMIT)
    op_ids = set()
    isolation_records = 0

    for operation in operations:
        if not isinstance(operation, dict):
            raise ValueError('operation must be an object')
        kind = operation.get('type')
        if kind not in _TYPES:
            raise ValueError('unknown operation type')
        _object(operation, _FIELDS[kind])

        op_id = ident(operation['op_id'])
        if op_id in op_ids:
            raise ValueError('duplicate op_id: ' + op_id)
        op_ids.add(op_id)

        if kind == 'stock':
            ident(operation['lot'])
            _parsed_formula(operation['formula'])
            _positive_quantity(operation['moles'])
        elif kind == 'run':
            inputs = _bounded_array(operation['inputs'], 1, _INPUT_LIMIT)
            input_lots = set()
            for item in inputs:
                _object(item, ('lot', 'moles'))
                lot = ident(item['lot'])
                if lot in input_lots:
                    raise ValueError('duplicate input lot: ' + lot)
                input_lots.add(lot)
                _positive_quantity(item['moles'])

            products = _bounded_array(operation['products'], 1, _PRODUCT_LIMIT)
            product_lots = set()
            for item in products:
                _object(item, ('lot', 'formula'))
                lot = ident(item['lot'])
                if lot in product_lots:
                    raise ValueError('duplicate product lot: ' + lot)
                product_lots.add(lot)
                _parsed_formula(item['formula'])
        else:
            ident(operation['lot'])
            # Rejected isolation records still consume the 500-record budget.
            isolation_records += 1

    if isolation_records > _HOLD_LIMIT:
        raise ValueError('too many hold/release records')
    return operations


class _Ledger:
    """Compact, non-recursive lot ledger.

    Lots live in append-only parallel arrays keyed by permanent integer
    indexes.  Lineage is stored once as child adjacency; blocked state is
    derived from the active roots instead of copying ancestor sets.
    """

    def __init__(self):
        self.ids = []
        self.formulas = []
        self.remaining = []
        self.children = []
        self.blocked = []
        self.index = {}
        self.roots = set()
        self.edge_count = 0

    def _resource_check(self, new_lots, new_edges):
        if len(self.ids) + new_lots > _LOT_LIMIT:
            raise ValueError('too many created lots')
        if self.edge_count + new_edges > _EDGE_LIMIT:
            raise ValueError('too many direct lineage edges')

    def _add_lot(self, lot, formula_text, amount, blocked=False):
        index = len(self.ids)
        self.ids.append(lot)
        self.formulas.append(formula_text)
        self.remaining.append(amount)
        self.children.append([])
        self.blocked.append(blocked)
        self.index[lot] = index
        return index

    def _recompute_blocked(self):
        """Mark exactly the descendants of every active root, iteratively."""
        blocked = [False] * len(self.ids)
        pending = []
        for root in self.roots:
            blocked[root] = True
            pending.append(root)
        while pending:
            node = pending.pop()
            for child in self.children[node]:
                if not blocked[child]:
                    blocked[child] = True
                    pending.append(child)
        self.blocked = blocked

    def stock(self, operation):
        lot = operation['lot']
        if lot in self.index:
            return 'rejected'
        self._resource_check(1, 0)
        self._add_lot(lot, operation['formula'], _positive_quantity(operation['moles']))
        return 'accepted'

    def run(self, operation):
        parent_indexes = []
        merged_amounts = {}
        merged_order = []

        for item in operation['inputs']:
            index = self.index.get(item['lot'])
            if index is None:
                return 'rejected'
            text = self.formulas[index]
            amount = _positive_quantity(item['moles'])
            parent_indexes.append(index)
            if text in merged_amounts:
                merged_amounts[text] += amount
            else:
                merged_amounts[text] = amount
                merged_order.append(text)

        # Exceeding the per-reaction species ceiling is an input-contract
        # violation (the supplemental rules classify it as a batch error),
        # even though the merged count is only known from current lot state.
        if len(merged_order) > _SPECIES_LIMIT:
            raise ValueError('too many merged reactant species')

        product_formulas = [item['formula'] for item in operation['products']]
        if len(set(product_formulas)) != len(product_formulas):
            return 'rejected'
        if set(product_formulas) & set(merged_order):
            return 'rejected'

        # Isolation and inventory are checked against the complete current
        # state before anything is deducted or created.
        if any(self.blocked[index] for index in parent_indexes):
            return 'rejected'
        if any(
            self.remaining[index] < _positive_quantity(item['moles'])
            for index, item in zip(parent_indexes, operation['inputs'])
        ):
            return 'rejected'

        try:
            coefficients, _ = _cached_balance(
                tuple(merged_order), tuple(product_formulas)
            )
        except ValueError:
            return 'rejected'

        extent = merged_amounts[merged_order[0]] / coefficients[0]
        if extent <= 0 or any(
            merged_amounts[text] != extent * coefficient
            for text, coefficient in zip(merged_order, coefficients)
        ):
            return 'rejected'

        product_amounts = [
            extent * coefficient
            for coefficient in coefficients[len(merged_order):]
        ]

        # A structurally valid but rejected run must not reserve the proposed
        # lot identities.  Collisions are therefore checked only after every
        # other reason to reject has been ruled out.
        if any(item['lot'] in self.index for item in operation['products']):
            return 'rejected'

        new_lots = len(operation['products'])
        new_edges = len(parent_indexes) * new_lots
        self._resource_check(new_lots, new_edges)

        # Commit as one unit.  Nothing below can fail after this point.
        for index, item in zip(parent_indexes, operation['inputs']):
            self.remaining[index] -= _positive_quantity(item['moles'])

        product_indexes = [
            self._add_lot(item['lot'], item['formula'], amount)
            for item, amount in zip(operation['products'], product_amounts)
        ]
        for parent in parent_indexes:
            self.children[parent].extend(product_indexes)
        self.edge_count += new_edges
        return 'accepted'

    def hold(self, operation):
        index = self.index.get(operation['lot'])
        if index is None or index in self.roots:
            return 'rejected'
        self.roots.add(index)
        self._recompute_blocked()
        return 'accepted'

    def release(self, operation):
        index = self.index.get(operation['lot'])
        if index is None or index not in self.roots:
            return 'rejected'
        self.roots.remove(index)
        self._recompute_blocked()
        return 'accepted'

    def result(self):
        lots = [
            {
                'lot': self.ids[index],
                'formula': self.formulas[index],
                'remaining_moles': str(self.remaining[index]),
                'blocked': self.blocked[index],
            }
            for index in range(len(self.ids))
        ]
        lots.sort(key=lambda item: item['lot'])
        return {
            'operations': [],
            'lots': lots,
            'active_holds': sorted(self.ids[index] for index in self.roots),
        }


def analyze(data):
    operations = _validate_operations(data)
    ledger = _Ledger()
    handlers = {
        'stock': ledger.stock,
        'run': ledger.run,
        'hold': ledger.hold,
        'release': ledger.release,
    }
    statuses = []
    for operation in operations:
        status = handlers[operation['type']](operation)
        statuses.append({'op_id': operation['op_id'], 'status': status})
    result = ledger.result()
    result['operations'] = statuses
    return result
