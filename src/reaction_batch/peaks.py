"""Shared-drift chromatographic peak identity reconciliation.

One batch holds a reference catalogue (known formula + zero-drift retention
window) and several runs (observed peaks + a shared instrument drift domain).
Inside one run every peak may be explained by at most one reference and every
reference explains at most one peak.  Peak ``p`` with time ``t`` may pair with
reference ``r`` whose window is ``[lo, hi]`` exactly when their formula element
counts are equal and ``lo + d <= t <= hi + d`` for the shared drift ``d``; the
pair is therefore feasible precisely on the closed drift interval
``[t - hi, t - lo]`` intersected with the run drift domain.

For a drift ``d`` an *explanation* is a matching of the bipartite graph of
feasible pairs; ``K`` is the largest explanation size over all drifts.  The
report carries ``K``, every drift that admits a ``K``-explanation (merged
closed intervals), and for every peak the union of reference identities it
takes in some ``K``-explanation plus whether some ``K``-explanation leaves it
unmatched.  A canonical witness explanation is also reported.

Algorithm (no enumeration of matchings or permutations)
-------------------------------------------------------
The edge set ``E(d)`` changes only at interval endpoints, so the drift domain
is swept left to right over *regions*: each endpoint point ``x`` and each open
segment between consecutive endpoints.  ``E(x)`` at a point is a superset of
the adjacent segment edge sets, so the matching size at a point dominates its
neighbours and the optimal-drift set is always a union of closed regions whose
closures merge into the reported closed intervals (single points survive).

A maximum matching is maintained incrementally: moving into a point only adds
edges, moving into a segment only removes edges, and each change is repaired
with breadth-first augmenting paths, so the sweep never recomputes from
scratch and never enumerates matchings.  For each region whose size equals the
current best, the *allowed* edges (those appearing in some maximum matching of
that region) and the *freeable* peaks are derived in linear time from one
maximum matching ``M``: orient unmatched edges left->right and matched edges
right->left, then an unmatched edge ``(u, v)`` is allowed iff ``u`` is
reachable from an unmatched left vertex, or ``v`` reaches an unmatched right
vertex, or both ends share a strongly connected component; a peak can be
unmatched iff it is reachable from an unmatched left vertex (itself included).
These characterisations are the standard alternating-path ones, so supports
come from whole regions at once rather than from enumerating explanations.

The witness uses the smallest optimal drift ``d*`` and greedily fixes the
lexicographically smallest assignment vector: peaks are processed in
``peak_id`` order, each receives the smallest ``reference_id`` that still
allows a ``K``-matching of the remaining graph (null only when no reference
does), which is exactly the documented minimum because null sorts after every
``reference_id``.
"""

from collections import deque

from .chemistry import formula, quantity
from .common import ident

_REFERENCE_LIMIT = 40
_RUN_LIMIT = 8
_PEAK_LIMIT = 40


def _object(value, fields):
    if not isinstance(value, dict) or set(value) != set(fields):
        raise ValueError('expected object fields: ' + ' '.join(fields))
    return value


def _array(value, minimum, maximum):
    if not isinstance(value, list) or not minimum <= len(value) <= maximum:
        raise ValueError('array length outside allowed bounds')
    return value


def _drift(value):
    """Drift bound: quantity syntax with one optional leading minus."""
    if not isinstance(value, str):
        raise ValueError('drift bound must be a string')
    if value.startswith('-'):
        return -quantity(value[1:])
    return quantity(value)


def _interval(value, parse, label):
    _object(value, ('lower', 'upper'))
    lower = parse(value['lower'])
    upper = parse(value['upper'])
    if lower > upper:
        raise ValueError(label + ' lower must not exceed upper')
    return lower, upper


class _Reference:
    __slots__ = ('reference_id', 'counts', 'lower', 'upper')

    def __init__(self, reference_id, counts, lower, upper):
        self.reference_id = reference_id
        self.counts = counts
        self.lower = lower
        self.upper = upper


class _Peak:
    __slots__ = ('peak_id', 'counts', 'time')

    def __init__(self, peak_id, counts, time):
        self.peak_id = peak_id
        self.counts = counts
        self.time = time


class _Run:
    __slots__ = ('run_id', 'drift_lower', 'drift_upper', 'peaks')

    def __init__(self, run_id, drift_lower, drift_upper, peaks):
        self.run_id = run_id
        self.drift_lower = drift_lower
        self.drift_upper = drift_upper
        self.peaks = peaks


def _parse(data):
    """Validate the whole batch before any computation or publication."""
    _object(data, ('references', 'runs'))
    references = []
    seen_references = set()
    for item in _array(data['references'], 1, _REFERENCE_LIMIT):
        _object(item, ('reference_id', 'formula', 'window'))
        reference_id = ident(item['reference_id'])
        if reference_id in seen_references:
            raise ValueError('duplicate reference_id: ' + reference_id)
        seen_references.add(reference_id)
        counts = formula(item['formula'])
        lower, upper = _interval(item['window'], quantity, 'window')
        references.append(_Reference(reference_id, counts, lower, upper))
    runs = []
    seen_runs = set()
    for item in _array(data['runs'], 1, _RUN_LIMIT):
        _object(item, ('run_id', 'drift', 'peaks'))
        run_id = ident(item['run_id'])
        if run_id in seen_runs:
            raise ValueError('duplicate run_id: ' + run_id)
        seen_runs.add(run_id)
        drift_lower, drift_upper = _interval(item['drift'], _drift, 'drift')
        peaks = []
        seen_peaks = set()
        for raw in _array(item['peaks'], 0, _PEAK_LIMIT):
            _object(raw, ('peak_id', 'formula', 'time'))
            peak_id = ident(raw['peak_id'])
            if peak_id in seen_peaks:
                raise ValueError('duplicate peak_id: ' + peak_id)
            seen_peaks.add(peak_id)
            peaks.append(
                _Peak(peak_id, formula(raw['formula']), quantity(raw['time'])))
        runs.append(_Run(run_id, drift_lower, drift_upper, peaks))
    return references, runs


class _Matching:
    """Bipartite graph with toggleable edges and an incremental maximum matching."""

    def __init__(self, left, right, edges):
        self.left = left
        self.right = right
        self.edges = edges
        self.adjacent = [[] for _ in range(left)]
        for index, (u, _v) in enumerate(edges):
            self.adjacent[u].append(index)
        self.active = [False] * len(edges)
        self.match_left = [None] * left
        self.match_right = [None] * right

    def activate(self, index):
        self.active[index] = True

    def deactivate(self, index):
        """Remove an edge; report whether the current matching lost it."""
        self.active[index] = False
        u, v = self.edges[index]
        if self.match_left[u] == v:
            self.match_left[u] = None
            self.match_right[v] = None
            return True
        return False

    def activate_all(self):
        for index in range(len(self.edges)):
            self.active[index] = True

    def size(self):
        return sum(1 for v in self.match_left if v is not None)

    def augment(self):
        """Restore maximality after edge changes via augmenting paths."""
        while self._augment_once():
            pass

    def _augment_once(self):
        """Find one augmenting path by multi-source BFS; apply and report it."""
        left, right = self.left, self.right
        parent = [-2] * (left + right)
        queue = deque()
        for u in range(left):
            if self.match_left[u] is None:
                parent[u] = -1
                queue.append(u)
        while queue:
            x = queue.popleft()
            if x < left:
                for index in self.adjacent[x]:
                    if not self.active[index]:
                        continue
                    _u, v = self.edges[index]
                    if self.match_left[x] == v:
                        continue
                    y = left + v
                    if parent[y] != -2:
                        continue
                    parent[y] = x
                    if self.match_right[v] is None:
                        self._splice(parent, x, v)
                        return True
                    queue.append(y)
            else:
                u = self.match_right[x - left]
                if u is not None and parent[u] == -2:
                    parent[u] = x
                    queue.append(u)
        return False

    def _splice(self, parent, u, v):
        """Flip the alternating path ending at free right vertex ``v``."""
        while True:
            self.match_left[u] = v
            self.match_right[v] = u
            if parent[u] == -1:
                return
            v = parent[u] - self.left
            u = parent[self.left + v]


def _reach(arcs, sources):
    seen = [False] * len(arcs)
    queue = deque()
    for start in sources:
        if not seen[start]:
            seen[start] = True
            queue.append(start)
    while queue:
        x = queue.popleft()
        for y in arcs[x]:
            if not seen[y]:
                seen[y] = True
                queue.append(y)
    return seen


def _components(arcs, reverse):
    """Iterative Kosaraju strongly connected components."""
    size = len(arcs)
    visited = [False] * size
    order = []
    for start in range(size):
        if visited[start]:
            continue
        visited[start] = True
        stack = [(start, 0)]
        while stack:
            node, step = stack[-1]
            if step < len(arcs[node]):
                stack[-1] = (node, step + 1)
                nxt = arcs[node][step]
                if not visited[nxt]:
                    visited[nxt] = True
                    stack.append((nxt, 0))
            else:
                order.append(node)
                stack.pop()
    comp = [-1] * size
    label = 0
    for start in reversed(order):
        if comp[start] != -1:
            continue
        comp[start] = label
        stack = [start]
        while stack:
            node = stack.pop()
            for nxt in reverse[node]:
                if comp[nxt] == -1:
                    comp[nxt] = label
                    stack.append(nxt)
        label += 1
    return comp


def _region_support(matching):
    """Allowed references per peak and freeable peaks of one region.

    Uses one maximum matching of the region: an unmatched edge is allowed iff
    its left end is reachable from a free left vertex, its right end reaches a
    free right vertex, or both ends share a strongly connected component of
    the alternating orientation; matched edges are always allowed.  A peak is
    freeable iff it is reachable from a free left vertex (itself included).
    """
    left, right = matching.left, matching.right
    arcs = [[] for _ in range(left + right)]
    reverse = [[] for _ in range(left + right)]
    for index, (u, v) in enumerate(matching.edges):
        if not matching.active[index]:
            continue
        if matching.match_left[u] == v:
            tail, head = left + v, u
        else:
            tail, head = u, left + v
        arcs[tail].append(head)
        reverse[head].append(tail)
    reach_left = _reach(
        arcs, [u for u in range(left) if matching.match_left[u] is None])
    to_free_right = _reach(
        reverse,
        [left + v for v in range(right) if matching.match_right[v] is None])
    comp = _components(arcs, reverse)
    allowed = [set() for _ in range(left)]
    for index, (u, v) in enumerate(matching.edges):
        if not matching.active[index]:
            continue
        if (matching.match_left[u] == v or reach_left[u]
                or to_free_right[left + v] or comp[u] == comp[left + v]):
            allowed[u].add(v)
    freeable = [matching.match_left[u] is None or reach_left[u]
                for u in range(left)]
    return allowed, freeable


def _build_edges(references, peaks, drift_lower, drift_upper):
    """Feasible pairs with their closed drift intervals clipped to the domain."""
    edges = []
    for u, peak in enumerate(peaks):
        for v, reference in enumerate(references):
            if peak.counts != reference.counts:
                continue
            lower = peak.time - reference.upper
            upper = peak.time - reference.lower
            if lower < drift_lower:
                lower = drift_lower
            if upper > drift_upper:
                upper = drift_upper
            if lower <= upper:
                edges.append((u, v, lower, upper))
    return edges


def _merge_regions(regions):
    """Merge ordered region closures; touching closed intervals coalesce."""
    merged = []
    for lower, upper in regions:
        if merged and lower <= merged[-1][1]:
            if upper > merged[-1][1]:
                merged[-1][1] = upper
        else:
            merged.append([lower, upper])
    return merged


def _witness(references, peaks, edges, drift, target):
    """Lexicographically smallest assignment vector at one optimal drift.

    Peaks are fixed in ``peak_id`` order; each takes the smallest
    ``reference_id`` that still leaves a ``target``-sized matching of the
    remaining graph, or null when every remaining maximum matching frees it.
    """
    left = len(peaks)
    right = len(references)
    base_edges = [(u, v) for (u, v, lower, upper) in edges
                  if lower <= drift <= upper]
    alive_left = [True] * left
    alive_right = [True] * right
    assignment = [None] * left
    remaining = target
    for u in sorted(range(left), key=lambda i: peaks[i].peak_id):
        current = [(a, b) for (a, b) in base_edges
                   if alive_left[a] and alive_right[b]]
        matching = _Matching(left, right, current)
        matching.activate_all()
        matching.augment()
        if matching.size() != remaining:
            raise RuntimeError('witness matching lost the optimal size')
        allowed, _freeable = _region_support(matching)
        alive_left[u] = False
        if allowed[u]:
            choice = min(allowed[u],
                         key=lambda v: references[v].reference_id)
            assignment[u] = choice
            alive_right[choice] = False
            remaining -= 1
    if remaining != 0:
        raise RuntimeError('witness did not reach the optimal size')
    return assignment


def _solve_run(references, run):
    peaks = run.peaks
    left = len(peaks)
    drift_lower, drift_upper = run.drift_lower, run.drift_upper
    edges = _build_edges(references, peaks, drift_lower, drift_upper)
    points = sorted({drift_lower, drift_upper}
                    | {edge[2] for edge in edges}
                    | {edge[3] for edge in edges})
    add_at = {point: [] for point in points}
    drop_at = {point: [] for point in points}
    for index, (_u, _v, lower, upper) in enumerate(edges):
        add_at[lower].append(index)
        drop_at[upper].append(index)
    matching = _Matching(left, len(references),
                         [(u, v) for (u, v, _a, _b) in edges])
    best = -1
    best_regions = []
    support = [set() for _ in range(left)]
    freeable = [False] * left

    def visit(region_lower, region_upper):
        nonlocal best, best_regions, support, freeable
        size = matching.size()
        if size < best:
            return
        allowed, can_free = _region_support(matching)
        if size > best:
            best = size
            best_regions = []
            support = [set() for _ in range(left)]
            freeable = [False] * left
        best_regions.append((region_lower, region_upper))
        for u in range(left):
            support[u] |= allowed[u]
            freeable[u] = freeable[u] or can_free[u]

    for position, point in enumerate(points):
        for index in add_at[point]:
            matching.activate(index)
        matching.augment()
        visit(point, point)
        if position + 1 < len(points):
            lost = False
            for index in drop_at[point]:
                lost = matching.deactivate(index) or lost
            if lost:
                matching.augment()
            visit(point, points[position + 1])

    best_drifts = _merge_regions(best_regions)
    assignment = _witness(references, peaks, edges, best_drifts[0][0], best)
    order = sorted(range(left), key=lambda i: peaks[i].peak_id)
    return {
        'run_id': run.run_id,
        'maximum_matches': best,
        'best_drifts': [{'lower': str(lower), 'upper': str(upper)}
                        for lower, upper in best_drifts],
        'peaks': [
            {'peak_id': peaks[u].peak_id,
             'reference_ids': sorted(references[v].reference_id
                                     for v in support[u]),
             'can_be_unmatched': freeable[u]}
            for u in order
        ],
        'witness': {
            'drift': str(best_drifts[0][0]),
            'assignments': [
                {'peak_id': peaks[u].peak_id,
                 'reference_id': (None if assignment[u] is None
                                  else references[assignment[u]].reference_id)}
                for u in order
            ],
        },
    }


def analyze(data):
    references, runs = _parse(data)
    return {'runs': [_solve_run(references, run) for run in runs]}
