"""Exact rational LP for the measurement-consistency assay.

A sample defines the closed polyhedron

    P = { x : A x <= b,  x >= 0 }

with ``n`` species columns and ``R`` canonical rows (``R`` grows with the number
of measurements).  The assay needs the exact minimum and maximum of every
column plus a certificate.  Solving the primal would need an ``R``-row
tableau, so this module solves the LP dual instead:

    maximise  -b.y   subject to   A^T y - s = c,   y >= 0, s >= 0.

Its row count is ``n`` no matter how many measurements the sample carries, and
the surplus block ``-I`` makes the coefficient matrix full rank, so phase one
can always eliminate every artificial variable.

Tableau (``n`` rows, max convention)::

    [ A^T (R columns) | -I (surplus, n) | I (artificial, n) ]

Rows are stored as integer numerators with one positive denominator per row,
so pivoting is fraction-free integer arithmetic with a periodic gcd reduction.
The artificial block starts as the identity and receives exactly the same row
operations, so it is ``B^-1`` at all times.  That yields the equality
multipliers ``pi = c_B B^-1`` and lets a new right-hand side be loaded as
``B^-1 d`` without refactorising.

Facts used by the callers:

* For ``c = e_j`` the dual optimum is ``max x_j``; for ``c = -e_j`` it is
  ``-min x_j``, and the primal witness is ``x = -pi``.
* The canonical rows include ``sum_j x_j <= total.upper`` and ``x >= 0``, so
  the recession cone of ``P`` is trivial.  By Farkas ``D(c)`` is then feasible
  for every ``c``.  If ``P`` is empty, ``D(c)`` is unbounded for every ``c``;
  the ray lies on ``y`` and ``s`` only (artificials never enter after phase
  one), so ``A^T y = s >= 0`` and ``b.y < 0``: a Farkas certificate.
* When ``P`` is nonempty every ``D(c)`` is bounded, so after the first optimal
  basis the remaining bounds use the dual simplex with the same cost vector
  and a different right-hand side.

The optimal ``y`` satisfies ``A^T y >= c`` and can exceed ``c`` on rows whose
surplus is basic.  Callers close that gap on the ``-x_j <= 0`` rows, whose
right-hand side is zero, so the reported certificate has ``A^T y = c`` exactly
with the same ``b.y``.

No floating point is used anywhere.  Dantzig pricing is used until a run of
degenerate pivots suggests stalling, then Bland's rule takes over so the
simplex cannot cycle.
"""

from fractions import Fraction
from math import gcd

ZERO = Fraction(0)
ONE = Fraction(1)

_BLAND_AFTER = 12


class _Unbounded(Exception):
    def __init__(self, ray):
        super().__init__('unbounded')
        self.ray = ray


class _Limit(RuntimeError):
    pass


def _common(numbers):
    divisor = 0
    for value in numbers:
        if value:
            divisor = gcd(divisor, value)
            if divisor == 1:
                return 1
    return divisor


class AssayLP:
    """Exact LP dual for ``A x <= b, x >= 0`` with warm starts.

    ``A`` is a list of ``R`` rows with ``n`` columns; ``b`` has ``R`` entries.
    Construction solves ``c = e_0`` and keeps the resulting optimal basis.

    If ``P`` is empty, ``infeasible`` is true and ``farkas`` holds
    nonnegative multipliers of the canonical rows with ``A^T farkas >= 0`` and
    ``b.farkas < 0`` (not yet normalised or closed).  Otherwise
    :meth:`bound` solves one coordinate.
    """

    __slots__ = (
        'n', 'r', 'width', 'art_start', 'T', 'den', 'rhs', 'basis',
        'cost', 'w', 'bounds', 'infeasible', 'farkas',
    )

    def __init__(self, A, b):
        n = len(A[0])
        r = len(A)
        self.n = n
        self.r = r
        self.art_start = r + n
        self.width = r + 2 * n
        T = []
        for i in range(n):
            row = []
            for k in range(r):
                value = Fraction(A[k][i])
                if value.denominator != 1:
                    raise ValueError('LP coefficient matrix must be integral')
                row.append(int(value))
            row += [-1 if k == i else 0 for k in range(n)]
            row += [1 if k == i else 0 for k in range(n)]
            T.append(row)
        self.T = T
        self.den = [1] * n
        self.bounds = [Fraction(value) for value in b]
        self.rhs = [1] + [0] * (n - 1)
        self.basis = [self.art_start + i for i in range(n)]
        self.cost = None
        self.w = None
        self.infeasible = False
        self.farkas = None
        # Phase one: artificial variables start basic; their cost is -1, so
        # the objective is -sum of artificial values, maximised at 0.
        self.w = self._reduced_cost(
            [ZERO] * self.art_start + [-ONE] * n)
        self._primal_simplex(budget=200000)
        for row in range(n):
            if self.basis[row] >= self.art_start:
                if self.rhs[row] != 0:
                    raise RuntimeError('phase one ended with a positive artificial')
                for j in range(self.art_start):
                    if self.T[row][j] != 0:
                        self._pivot(row, j)
                        break
                else:
                    raise RuntimeError('artificial variable cannot be removed')
        # Phase two: same right-hand side, real cost -b on the y block.
        self.cost = ([-value for value in self.bounds]
                     + [ZERO] * (self.width - r))
        self.w = self._reduced_cost(self.cost)
        try:
            self._primal_simplex()
        except _Unbounded as unbounded:
            self.infeasible = True
            self.farkas = [unbounded.ray[k] for k in range(self.r)]

    # -- tableau primitives -------------------------------------------------
    def _pivot(self, row, col):
        T = self.T
        den = self.den
        rhs = self.rhs
        pivot = T[row][col]
        if pivot < 0:
            T[row] = [-value for value in T[row]]
            rhs[row] = -rhs[row]
            pivot = -pivot
        pivot_row = T[row]
        # Row value = T[row] / den[row]; dividing by pivot / den[row] leaves
        # T[row] / pivot.
        den[row] = pivot
        w = self.w
        if w is not None:
            factor = w[col] / pivot
            if factor:
                self.w = [a - factor * b for a, b in zip(w, pivot_row)]
        self.basis[row] = col
        for i in range(self.n):
            if i == row:
                continue
            current = T[i]
            factor = current[col]
            if factor:
                updated = [a * pivot - factor * b
                           for a, b in zip(current, pivot_row)]
                new_rhs = rhs[i] * pivot - factor * rhs[row]
                denominator = den[i] * pivot
                if denominator < 0:
                    updated = [-value for value in updated]
                    new_rhs = -new_rhs
                    denominator = -denominator
                divisor = gcd(denominator, gcd(abs(new_rhs), _common(updated)))
                if divisor > 1:
                    updated = [value // divisor for value in updated]
                    new_rhs //= divisor
                    denominator //= divisor
                T[i] = updated
                rhs[i] = new_rhs
                den[i] = denominator
            # else the row is unchanged.

    def _reduced_cost(self, cost):
        w = list(cost)
        start = self.art_start
        for row in range(self.n):
            coefficient = cost[self.basis[row]]
            if coefficient:
                factor = coefficient / self.den[row]
                Trow = self.T[row]
                w = [a - factor * b for a, b in zip(w, Trow)]
        return w

    def _entering(self, bland, upper):
        w = self.w
        basic = set(self.basis)
        enter = -1
        best = ZERO
        for j in range(upper):
            if j in basic:
                continue
            value = w[j]
            if value > 0:
                if bland:
                    return j
                if value > best:
                    best = value
                    enter = j
        return enter

    def _primal_simplex(self, budget=400000):
        T = self.T
        den = self.den
        rhs = self.rhs
        basis = self.basis
        upper = self.art_start
        bland = False
        stalled = 0
        for _ in range(budget):
            enter = self._entering(bland, upper)
            if enter < 0:
                return
            leave = -1
            best_ratio = None
            best_basic = None
            for i in range(self.n):
                numerator = T[i][enter]
                if numerator > 0:
                    ratio = Fraction(rhs[i], numerator)
                    if (best_ratio is None or ratio < best_ratio
                            or (ratio == best_ratio
                                and basis[i] < best_basic)):
                        best_ratio = ratio
                        leave = i
                        best_basic = basis[i]
            if leave < 0:
                ray = [ZERO] * self.width
                ray[enter] = ONE
                for i in range(self.n):
                    ray[basis[i]] = -Fraction(T[i][enter], den[i])
                raise _Unbounded(ray)
            stalled = stalled + 1 if best_ratio == 0 else 0
            if stalled >= _BLAND_AFTER:
                bland = True
            self._pivot(leave, enter)
        raise _Limit('primal simplex iteration limit')

    def _load_rhs(self, d):
        # Rows are kept as T[i]/den[i]; the right-hand side shares the same
        # denominator, so B^-1 d is just the numerator sum below.
        T = self.T
        start = self.art_start
        rhs = []
        for i in range(self.n):
            row = T[i]
            total = 0
            for k in range(self.n):
                value = d[k]
                if value:
                    total += value * row[start + k]
            rhs.append(total)
        self.rhs = rhs

    def _dual_simplex(self, budget=400000):
        """Restore primal feasibility keeping reduced costs ``<= 0``.

        The leaving row is the most infeasible basic row.  The entering column
        minimises ``w_j / a_rj`` over nonbasic ``j`` with ``a_rj < 0``; because
        ``w_j <= 0`` the minimum is nonnegative and keeps every reduced cost
        nonpositive.  Basic columns are excluded: a basic variable is already
        in the basis, and letting it "re-enter" would make the basis singular.
        Bland's rule on both choices breaks ties and prevents cycling.
        """
        T = self.T
        den = self.den
        rhs = self.rhs
        bland = False
        stalled = 0
        for _ in range(budget):
            w = self.w
            leave = -1
            for i in range(self.n):
                if rhs[i] < 0:
                    if bland:
                        if leave < 0 or self.basis[i] < self.basis[leave]:
                            leave = i
                    elif leave < 0 or rhs[i] * den[leave] < rhs[leave] * den[i]:
                        leave = i
            if leave < 0:
                return
            enter = -1
            best = None
            best_index = None
            scale = den[leave]
            basic = set(self.basis)
            for j in range(self.art_start):
                if j in basic:
                    continue
                numerator = T[leave][j]
                if numerator < 0:
                    ratio = w[j] * scale / numerator
                    if (best is None or ratio < best
                            or (ratio == best and j < best_index)):
                        best = ratio
                        best_index = j
                        enter = j
            if enter < 0:
                raise RuntimeError('dual simplex lost dual feasibility')
            stalled = stalled + 1 if best == 0 else 0
            if stalled >= _BLAND_AFTER:
                bland = True
            self._pivot(leave, enter)
        raise _Limit('dual simplex iteration limit')

    def bound(self, index, sign):
        """Return ``(x, value, y)`` for ``max sign * x[index]``.

        ``x`` is a feasible primal witness of length ``n``; ``value`` is
        ``sign * x[index]``; ``y`` has length ``R`` with ``A^T y >= c`` and
        ``b.y == value``.  ``None`` means the polyhedron is empty.
        """
        if self.infeasible:
            return None
        d = [0] * self.n
        d[index] = sign
        self._load_rhs(d)
        self._dual_simplex()
        den = self.den
        y = [ZERO] * self.r
        for row, column in enumerate(self.basis):
            if column < self.r:
                y[column] = Fraction(self.rhs[row], den[row])
        objective = ZERO
        for row, column in enumerate(self.basis):
            coefficient = self.cost[column]
            if coefficient:
                objective += coefficient * Fraction(self.rhs[row], den[row])
        pi = []
        start = self.art_start
        for k in range(self.n):
            column = start + k
            total = ZERO
            for row in range(self.n):
                coefficient = self.cost[self.basis[row]]
                if coefficient:
                    total += coefficient * Fraction(
                        self.T[row][column], den[row])
            pi.append(total)
        x = [-value for value in pi]
        return x, -objective, y
