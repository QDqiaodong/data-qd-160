"""Leakage-safe train/test planning for chemistry classification samples.

Two records are linked when parsing their ``formula`` yields exactly the same
element counts (so ``H2O`` and ``OH2`` are linked while ``H4O2`` is not), or
when they share an ``experiment_id``.  The undirected transitive closure of
those two relations forms indivisible groups; labels never cut a group apart.

Within that fixed grouping the planner honours every per-record ``pin``,
then chooses a whole-group assignment with both test-label counts inside
their closed quota intervals.  Among all such assignments it minimises
``|P-positive.target| + |N-negative.target|`` and, on a tie, the 0/1
assignment vector over record ids in Unicode order (train=0, test=1).

The optimiser keeps, for every suffix of free groups, the set of reachable
``(positive, negative)`` record counts as one integer per positive count,
whose bits are the reachable negative counts.  That is a bounded state space
of 301 x 301 counts regardless of how many subsets exist.  The canonical
assignment is rebuilt left to right over groups ordered by their minimum
record id, preferring ``train`` whenever the suffix reachability set still
admits an optimal completion; no subset enumeration, random seed or greedy
count choice is involved.
"""

from .chemistry import formula
from .common import ident

_RECORD_LIMIT = 600
_LABEL_LIMIT = 300
_QUOTA_LIMIT = 300
_COUNT_MAX = 300

_LABELS = ('positive', 'negative')
_PINS = ('free', 'train', 'test')


def _object(value, fields):
    if not isinstance(value, dict) or set(value) != set(fields):
        raise ValueError('expected object fields: ' + ' '.join(fields))
    return value


def _integer(value, minimum, maximum, label):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(label + ' must be an integer within its bounds')
    return value


def _quota(value, label):
    _object(value, ('lower', 'upper', 'target'))
    lower = _integer(value['lower'], 0, _QUOTA_LIMIT, label + '.lower')
    upper = _integer(value['upper'], 0, _QUOTA_LIMIT, label + '.upper')
    target = _integer(value['target'], 0, _QUOTA_LIMIT, label + '.target')
    if lower > upper:
        raise ValueError(label + ' lower exceeds upper')
    return {'lower': lower, 'upper': upper, 'target': target}


def _parse(data):
    """Validate the entire batch before any grouping or planning."""
    _object(data, ('records', 'test_quotas'))
    raw_records = data['records']
    if not isinstance(raw_records, list) or len(raw_records) > _RECORD_LIMIT:
        raise ValueError('records must be an array of at most 600 entries')
    quotas = _object(data['test_quotas'], ('positive', 'negative'))
    quota = {'positive': _quota(quotas['positive'], 'positive'),
             'negative': _quota(quotas['negative'], 'negative')}

    records = []
    seen = set()
    label_totals = {label: 0 for label in _LABELS}
    for raw in raw_records:
        _object(raw, ('record_id', 'formula', 'experiment_id', 'label', 'pin'))
        record_id = ident(raw['record_id'])
        if record_id in seen:
            raise ValueError('duplicate record_id: ' + record_id)
        seen.add(record_id)
        counts = formula(raw['formula'])
        experiment_id = ident(raw['experiment_id'])
        label = raw['label']
        if label not in _LABELS:
            raise ValueError('label must be positive or negative')
        pin = raw['pin']
        if pin not in _PINS:
            raise ValueError('pin must be free, train or test')
        label_totals[label] += 1
        if label_totals[label] > _LABEL_LIMIT:
            raise ValueError('too many records for one label')
        records.append({'record_id': record_id,
                        'counts': tuple(sorted(counts.items())),
                        'experiment_id': experiment_id,
                        'label': label,
                        'pin': pin})
    return records, quota


class _DisjointSet:
    def __init__(self, size):
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, item):
        parent = self.parent
        root = item
        while parent[root] != root:
            root = parent[root]
        while parent[item] != root:
            parent[item], item = root, parent[item]
        return root

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left == right:
            return
        if self.rank[left] < self.rank[right]:
            left, right = right, left
        self.parent[right] = left
        if self.rank[left] == self.rank[right]:
            self.rank[left] += 1


def _components(records):
    """Transitive closure over equal element counts and equal experiment id."""
    dsu = _DisjointSet(len(records))
    first_by_formula = {}
    first_by_experiment = {}
    for index, record in enumerate(records):
        counts = record['counts']
        if counts in first_by_formula:
            dsu.union(index, first_by_formula[counts])
        else:
            first_by_formula[counts] = index
        experiment = record['experiment_id']
        if experiment in first_by_experiment:
            dsu.union(index, first_by_experiment[experiment])
        else:
            first_by_experiment[experiment] = index

    members = {}
    for index in range(len(records)):
        members.setdefault(dsu.find(index), []).append(index)
    components = []
    for indices in members.values():
        record_ids = sorted(records[index]['record_id'] for index in indices)
        positive = sum(1 for index in indices
                       if records[index]['label'] == 'positive')
        negative = len(indices) - positive
        pins = {records[index]['pin'] for index in indices}
        if 'train' in pins and 'test' in pins:
            forced = 'conflict'
        elif 'train' in pins:
            forced = 'train'
        elif 'test' in pins:
            forced = 'test'
        else:
            forced = 'free'
        components.append({'group_id': record_ids[0],
                           'record_ids': record_ids,
                           'indices': indices,
                           'positive': positive,
                           'negative': negative,
                           'forced': forced})
    components.sort(key=lambda item: item['group_id'])
    return components


def _group_reports(components):
    return [{'group_id': component['group_id'],
             'record_ids': list(component['record_ids']),
             'counts': {'positive': component['positive'],
                        'negative': component['negative']},
             'forced': component['forced']}
            for component in components]


def _suffix_reachability(free):
    """``suffix[i]``: reachable count pairs using free groups ``i..``.

    One row per positive count holds a bit per negative count.  Rows are
    immutable integers, so snapshots share unchanged rows and each snapshot
    stays at most ``301 * 301`` reachable bits, independent of group count.
    """
    count = len(free)
    suffix = [None] * (count + 1)
    empty = [0] * (_COUNT_MAX + 1)
    empty[0] = 1  # the always-reachable ``(0, 0)`` suffix
    suffix[count] = tuple(empty)
    keep = (1 << (_COUNT_MAX + 1)) - 1
    for index in range(count - 1, -1, -1):
        positive, negative = free[index]['positive'], free[index]['negative']
        previous = suffix[index + 1]
        rows = list(previous)
        last = _COUNT_MAX - positive
        for base in range(last + 1):
            rows[base + positive] |= (previous[base] << negative) & keep
        suffix[index] = tuple(rows)
    return suffix


def _acceptable_pairs(quota, deviation):
    """All ``(P, N)`` inside the quotas whose deviation equals ``deviation``."""
    positive = quota['positive']
    negative = quota['negative']
    pairs = []
    for p_count in range(positive['lower'], positive['upper'] + 1):
        first = abs(p_count - positive['target'])
        if first > deviation:
            continue
        second = deviation - first
        candidates = {negative['target'] - second,
                      negative['target'] + second}
        for n_count in sorted(candidates):
            if (negative['lower'] <= n_count <= negative['upper']
                    and abs(n_count - negative['target']) == second):
                pairs.append((p_count, n_count))
    return pairs


def _completable(reachable, base_p, base_n, pairs):
    """Whether some suffix pair completes an optimal, in-quota assignment."""
    for p_count, n_count in pairs:
        p_extra = p_count - base_p
        if not 0 <= p_extra <= _COUNT_MAX:
            continue
        n_extra = n_count - base_n
        if not 0 <= n_extra <= _COUNT_MAX:
            continue
        if (reachable[p_extra] >> n_extra) & 1:
            return True
    return False


def _best_deviation(reachable, base_p, base_n, quota):
    best = None
    positive = quota['positive']
    negative = quota['negative']
    for p_extra in range(_COUNT_MAX + 1):
        p_count = base_p + p_extra
        if p_count > positive['upper']:
            break
        if p_count < positive['lower']:
            continue
        row = reachable[p_extra]
        if not row:
            continue
        for n_extra in range(_COUNT_MAX + 1):
            n_count = base_n + n_extra
            if n_count > negative['upper']:
                break
            if n_count < negative['lower']:
                continue
            if not ((row >> n_extra) & 1):
                continue
            deviation = (abs(p_count - positive['target'])
                         + abs(n_count - negative['target']))
            if best is None or deviation < best:
                best = deviation
    return best


def analyze(data):
    records, quota = _parse(data)
    components = _components(records)
    groups = _group_reports(components)

    conflicts = sorted(component['group_id'] for component in components
                       if component['forced'] == 'conflict')
    if conflicts:
        return {'groups': groups, 'status': 'infeasible',
                'reason': 'pin_conflict', 'conflict_groups': conflicts}

    base_positive = sum(component['positive'] for component in components
                        if component['forced'] == 'test')
    base_negative = sum(component['negative'] for component in components
                        if component['forced'] == 'test')
    free = [component for component in components
            if component['forced'] == 'free']
    suffix = _suffix_reachability(free)
    reachable = suffix[0]

    deviation = _best_deviation(reachable, base_positive, base_negative, quota)
    if deviation is None:
        return {'groups': groups, 'status': 'infeasible', 'reason': 'quotas'}

    pairs = _acceptable_pairs(quota, deviation)
    chosen_positive = chosen_negative = 0
    rank = 0
    side = {}
    for component in components:
        if component['forced'] == 'train':
            selected = False
        elif component['forced'] == 'test':
            selected = True
        else:
            selected = False
            if not _completable(suffix[rank + 1],
                                base_positive + chosen_positive,
                                base_negative + chosen_negative, pairs):
                selected = True
                chosen_positive += component['positive']
                chosen_negative += component['negative']
            rank += 1
        for index in component['indices']:
            side[index] = selected

    train = []
    test = []
    test_positive = test_negative = 0
    for index in sorted(range(len(records)),
                        key=lambda item: records[item]['record_id']):
        record_id = records[index]['record_id']
        if side[index]:
            test.append(record_id)
            if records[index]['label'] == 'positive':
                test_positive += 1
            else:
                test_negative += 1
        else:
            train.append(record_id)

    return {'groups': groups, 'status': 'feasible',
            'split': {'train': train, 'test': test,
                      'test_counts': {'positive': test_positive,
                                      'negative': test_negative},
                      'deviation': (abs(test_positive - quota['positive']['target'])
                                    + abs(test_negative - quota['negative']['target']))}}
