"""Sample composition consistency under element-wise molar measurements.

Every sample asks whether its measured element totals can be explained by one
nonnegative mixture of the candidate species, and if so which mole amount each
species can take.  The question is a polyhedron

    P = { x in Q^n : A x <= b,  x >= 0 }

over the ``n`` candidate species, with the canonical row order

    0 .. n-1            -x_j <= 0
    n                   sum_j x_j <= total_moles.upper
    n+1                 -sum_j x_j <= -total_moles.lower
    n+2+2k              a_e x <= upper_k
    n+3+2k              -a_e x <= -lower_k          (k in input order)

where ``a_ej`` is the number of atoms of element ``e`` in species ``j`` and the
measured values are element *atom* moles, not species moles.  ``m`` is the
number of measurements of that sample, so ``R = n + 2 + 2m``.

``P`` is always contained in the bounded set ``sum_j x_j <= total.upper``, so
each coordinate has a finite minimum and maximum.  All values are exact
:class:`fractions.Fraction`; measurement midpoints are never used, and the two
bounds of every species are computed jointly over the whole feasible region.

Parallel rows (several measurements of one element) encode the same quantity,
so for *solving* each element keeps only the tightest lower and upper rows;
the resulting multipliers are mapped back onto the original canonical rows so
that every certificate checks out against the exact ``A``/``b`` from the input.

Certificates
------------
For a maximum the objective is ``c = +e_j``: the report carries a feasible
composition ``x`` (``A x <= b``, ``x >= 0``) with ``x_j = value`` and a
nonnegative ``dual`` with ``A^T dual = c`` and ``b . dual = value``.  For a
minimum the objective is ``c = -e_j``, the composition attains ``x_j = value``,
and the dual satisfies ``A^T dual = -e_j`` and ``b . dual = -value``.  An
inconsistent sample instead carries a nonnegative ``contradiction`` with
``A^T contradiction = 0`` and ``b . contradiction = -1``.

The gap left by the dual simplex on ``-x_j <= 0`` rows is closed exactly by
adding it to those rows: their right-hand side is zero, so ``b . dual`` is
unchanged, and row ``j`` only touches column ``j``.
"""

from fractions import Fraction

from .assay_lp import AssayLP
from .chemistry import formula, quantity
from .common import ident, rows

_SPECIES_LIMIT = 60
_SAMPLE_LIMIT = 8
_MEASUREMENT_LIMIT = 120

ZERO = Fraction(0)
ONE = Fraction(1)


class _Measurement:
    __slots__ = ('measurement_id', 'element', 'lower', 'upper', 'index')

    def __init__(self, measurement_id, element, lower, upper, index):
        self.measurement_id = measurement_id
        self.element = element
        self.lower = lower
        self.upper = upper
        self.index = index


class _Sample:
    __slots__ = ('sample_id', 'lower_total', 'upper_total', 'measurements')

    def __init__(self, sample_id, lower_total, upper_total, measurements):
        self.sample_id = sample_id
        self.lower_total = lower_total
        self.upper_total = upper_total
        self.measurements = measurements


def _object(value, fields):
    if not isinstance(value, dict) or set(value) != set(fields):
        raise ValueError('expected object fields: ' + ' '.join(fields))
    return value


def _array(value, limit):
    if not isinstance(value, list) or len(value) > limit:
        raise ValueError('array length outside allowed bounds')
    return value


def _element(value):
    if not isinstance(value, str) or len(value) not in (1, 2):
        raise ValueError('invalid element token')
    if not 'A' <= value[0] <= 'Z':
        raise ValueError('invalid element token')
    if len(value) == 2 and not 'a' <= value[1] <= 'z':
        raise ValueError('invalid element token')
    return value


def _interval(value):
    _object(value, ('lower', 'upper'))
    low = quantity(value['lower'])
    high = quantity(value['upper'])
    if low > high:
        raise ValueError('interval lower must not exceed upper')
    return low, high


def _parse(data):
    """Validate the whole batch before any computation or publication."""
    _object(data, ('species', 'samples'))
    species_rows = _array(data['species'], _SPECIES_LIMIT)
    if not species_rows:
        raise ValueError('at least one species is required')
    species_ids = []
    compositions = []
    seen_ids = set()
    seen_formulas = set()
    for item in species_rows:
        _object(item, ('species_id', 'formula'))
        species_id = ident(item['species_id'])
        if species_id in seen_ids:
            raise ValueError('duplicate species_id: ' + species_id)
        seen_ids.add(species_id)
        text = item['formula']
        if not isinstance(text, str):
            raise ValueError('formula must be a string')
        if text in seen_formulas:
            raise ValueError('duplicate formula text: ' + text)
        seen_formulas.add(text)
        species_ids.append(species_id)
        compositions.append(formula(text))
    sample_rows = _array(data['samples'], _SAMPLE_LIMIT)
    if not sample_rows:
        raise ValueError('at least one sample is required')
    samples = []
    seen_samples = set()
    for item in sample_rows:
        _object(item, ('sample_id', 'total_moles', 'measurements'))
        sample_id = ident(item['sample_id'])
        if sample_id in seen_samples:
            raise ValueError('duplicate sample_id: ' + sample_id)
        seen_samples.add(sample_id)
        lower_total, upper_total = _interval(item['total_moles'])
        measurement_list = _array(item['measurements'], _MEASUREMENT_LIMIT)
        measurement_ids = set()
        measurements = []
        for index, raw in enumerate(measurement_list):
            _object(raw, ('measurement_id', 'element', 'lower', 'upper'))
            measurement_id = ident(raw['measurement_id'])
            if measurement_id in measurement_ids:
                raise ValueError('duplicate measurement_id: ' + measurement_id)
            measurement_ids.add(measurement_id)
            element = _element(raw['element'])
            low = quantity(raw['lower'])
            high = quantity(raw['upper'])
            if low > high:
                raise ValueError('measurement interval lower exceeds upper')
            measurements.append(
                _Measurement(measurement_id, element, low, high, index))
        samples.append(_Sample(sample_id, lower_total, upper_total, measurements))
    return species_ids, compositions, samples


def _element_order(sample):
    ordered = []
    seen = set()
    for measurement in sample.measurements:
        if measurement.element not in seen:
            seen.add(measurement.element)
            ordered.append(measurement.element)
    return ordered


def _reduced_system(compositions, sample):
    """Tighten parallel measurement rows; return rows, bounds, canonical index.

    Rows equal to another row are redundant.  Keeping the tightest bound per
    element leaves the *set* of primal optima unchanged, and the dropped rows
    receive zero multipliers when the certificate is spread back, so the
    canonical ``A``/``b`` are still satisfied exactly.
    """
    n = len(compositions)
    vector = []
    for j in range(n):
        vector.append([-ONE if k == j else ZERO for k in range(n)])
    bounds = [ZERO] * n
    canonical = list(range(n))
    vector.append([ONE] * n)
    bounds.append(sample.upper_total)
    canonical.append(n)
    vector.append([-ONE] * n)
    bounds.append(-sample.lower_total)
    canonical.append(n + 1)
    # Tightest lower and upper per element; ties keep the earliest measurement.
    tightest = {}
    for measurement in sample.measurements:
        entry = tightest.get(measurement.element)
        if entry is None:
            tightest[measurement.element] = {
                'lower': measurement.lower, 'lower_index': measurement.index,
                'upper': measurement.upper, 'upper_index': measurement.index,
            }
            continue
        if (measurement.lower > entry['lower']
                or (measurement.lower == entry['lower']
                    and measurement.index < entry['lower_index'])):
            entry['lower'] = measurement.lower
            entry['lower_index'] = measurement.index
        if (measurement.upper < entry['upper']
                or (measurement.upper == entry['upper']
                    and measurement.index < entry['upper_index'])):
            entry['upper'] = measurement.upper
            entry['upper_index'] = measurement.index
    for element in _element_order(sample):
        entry = tightest[element]
        coefficients = [Fraction(composition.get(element, 0))
                        for composition in compositions]
        vector.append(list(coefficients))
        bounds.append(entry['upper'])
        canonical.append(n + 2 + 2 * entry['upper_index'])
        vector.append([-value for value in coefficients])
        bounds.append(-entry['lower'])
        canonical.append(n + 3 + 2 * entry['lower_index'])
    return vector, bounds, canonical


def _spread(multipliers, canonical, row_count):
    result = [ZERO] * row_count
    for value, row in zip(multipliers, canonical):
        result[row] += value
    return result


def _column_sums(matrix, multipliers):
    width = len(matrix[0])
    sums = [ZERO] * width
    for coefficients, multiplier in zip(matrix, multipliers):
        if multiplier:
            for k in range(width):
                if coefficients[k]:
                    sums[k] += coefficients[k] * multiplier
    return sums


def _close_dual(multipliers, matrix, canonical, objective, row_count):
    """Make ``A^T dual = objective`` exact using the ``-x_j <= 0`` rows.

    Those rows have right-hand side zero and row ``j`` is ``-e_j``, so raising
    multiplier ``j`` by the gap on column ``j`` touches nothing else and does
    not change ``b . dual``.
    """
    gaps = _column_sums(matrix, multipliers)
    combined = list(multipliers)
    for j in range(len(objective)):
        gap = gaps[j] - objective[j]
        if gap < 0:
            raise RuntimeError('dual violated a column constraint')
        combined[j] += gap
    return _spread(combined, canonical, row_count)


def _canonical_bounds(sample, n):
    bounds = [ZERO] * n
    bounds.append(sample.upper_total)
    bounds.append(-sample.lower_total)
    for measurement in sample.measurements:
        bounds.append(measurement.upper)
        bounds.append(-measurement.lower)
    return bounds


def _farkas(multipliers, matrix, canonical, row_count, canonical_bounds):
    """Normalise and close a Farkas ray to ``A^T c = 0``, ``b . c = -1``."""
    gaps = _column_sums(matrix, multipliers)
    combined = list(multipliers)
    for j, gap in enumerate(gaps):
        if gap < 0:
            raise RuntimeError('Farkas ray violated a column sign')
        combined[j] += gap
    spread = _spread(combined, canonical, row_count)
    product = sum((value * bound for value, bound
                   in zip(spread, canonical_bounds)), ZERO)
    if product >= 0:
        raise RuntimeError('Farkas certificate does not refute the bounds')
    return [value / (-product) for value in spread]


def _solve_sample(compositions, sample):
    n = len(compositions)
    row_count = n + 2 + 2 * len(sample.measurements)
    matrix, bounds, canonical = _reduced_system(compositions, sample)
    program = AssayLP(matrix, bounds)
    if program.infeasible:
        contradiction = _farkas(
            program.farkas, matrix, canonical, row_count,
            _canonical_bounds(sample, n))
        return {'sample_id': sample.sample_id, 'status': 'inconsistent',
                'contradiction': [str(value) for value in contradiction]}
    ranges = []
    for j in range(n):
        maximum = program.bound(j, 1)
        minimum = program.bound(j, -1)
        upper_x, upper_value, upper_y = maximum
        lower_x, lower_value, lower_y = minimum
        upper_objective = [ONE if k == j else ZERO for k in range(n)]
        lower_objective = [-ONE if k == j else ZERO for k in range(n)]
        upper_dual = _close_dual(
            upper_y, matrix, canonical, upper_objective, row_count)
        lower_dual = _close_dual(
            lower_y, matrix, canonical, lower_objective, row_count)
        if upper_x[j] != upper_value:
            raise RuntimeError('maximising witness does not attain its value')
        if lower_x[j] != -lower_value:
            raise RuntimeError('minimising witness does not attain its value')
        ranges.append({
            'species_id': None,
            'minimum': {
                'value': str(-lower_value),
                'composition': [str(value) for value in lower_x],
                'dual': [str(value) for value in lower_dual],
            },
            'maximum': {
                'value': str(upper_value),
                'composition': [str(value) for value in upper_x],
                'dual': [str(value) for value in upper_dual],
            },
        })
    return {'sample_id': sample.sample_id, 'status': 'consistent',
            'ranges': ranges}


def analyze(data):
    species_ids, compositions, samples = _parse(data)
    results = []
    for sample in samples:
        item = _solve_sample(compositions, sample)
        if item['status'] == 'consistent':
            for position, entry in enumerate(item['ranges']):
                entry['species_id'] = species_ids[position]
        results.append(item)
    return {'species_order': list(species_ids), 'samples': results}
