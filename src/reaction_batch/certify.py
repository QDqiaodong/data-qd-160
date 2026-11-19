"""Traceable piecewise instrument calibration, corrected observations and
protected certificates for offline teaching concentration records.

The input is an ordered stream of batches, each carrying simultaneous
calibration / observation / certificate puts and deletes, plus historical
queries addressed by batch prefix.  Every batch is applied to a candidate view
on top of the last accepted state; the candidate is committed only when every
active certificate still passes, so a retroactive calibration correction that
would silently invalidate a protected certificate is rejected as a whole and
the previous state survives untouched.

All arithmetic is exact :class:`~fractions.Fraction` arithmetic.  Calibration
curves are piecewise linear over strictly increasing ``raw`` nodes on the
half-open validity interval ``[start, end)``; a time covered by no interval is
``no_calibration`` and a raw value outside the selected curve's node domain is
``outside_curve``.  Only the prefix views actually named by a query are
retained as snapshots; records are immutable, so a shallow dictionary copy
freezes a view forever and later retroactive corrections can never rewrite an
earlier prefix's values, sources or certificate membership.
"""

from bisect import bisect_right
from fractions import Fraction

from .chemistry import formula, quantity
from .common import ident

_BATCH_LIMIT = 300
_CHANGE_LIMIT = 6000
_QUERY_LIMIT = 60
_CALIBRATION_NODE_LIMIT = 20000
_CERTIFICATE_MEMBER_LIMIT = 2000
_TIME_MAX = 10**9
_POINT_MIN = 2
_POINT_MAX = 64
_MEMBER_MIN = 1
_MEMBER_MAX = 64

_ENTITIES = ('calibration', 'observation', 'certificate')
_REASON_PRIORITY = (
    'missing_delete',
    'overlap',
    'missing_member',
    'formula_mismatch',
    'unavailable_member',
    'mean_out_of_range',
    'range_exceeded',
)


def _object(value, fields):
    if not isinstance(value, dict) or set(value) != set(fields):
        raise ValueError('expected object fields: ' + ' '.join(fields))
    return value


def _integer(value, minimum, maximum, label):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(label + ' must be an integer within its bounds')
    return value


def _amount(value, label):
    amount = quantity(value)
    if amount < 0:
        raise ValueError(label + ' must be nonnegative')
    return amount


def _calibration(record):
    _object(record, ('id', 'instrument', 'start', 'end', 'points'))
    cal_id = ident(record['id'])
    instrument = ident(record['instrument'])
    start = _integer(record['start'], 0, _TIME_MAX, 'calibration.start')
    end = _integer(record['end'], 0, _TIME_MAX, 'calibration.end')
    if start >= end:
        raise ValueError('calibration interval must satisfy start < end')
    points_raw = record['points']
    if not isinstance(points_raw, list) or not _POINT_MIN <= len(points_raw) <= _POINT_MAX:
        raise ValueError('calibration points must be an array of 2..64 nodes')
    points = []
    previous = None
    for node in points_raw:
        _object(node, ('raw', 'value'))
        raw = _amount(node['raw'], 'calibration raw')
        value = _amount(node['value'], 'calibration value')
        if previous is not None and raw <= previous:
            raise ValueError('calibration raw nodes must be strictly increasing')
        previous = raw
        points.append({'raw': raw, 'value': value})
    return {'id': cal_id, 'instrument': instrument, 'start': start,
            'end': end, 'points': points}


def _observation(record):
    _object(record, ('id', 'instrument', 'time', 'formula', 'raw'))
    obs_id = ident(record['id'])
    instrument = ident(record['instrument'])
    time = _integer(record['time'], 0, _TIME_MAX, 'observation.time')
    formula_text = record['formula']
    counts = tuple(sorted(formula(formula_text).items()))
    raw = _amount(record['raw'], 'observation raw')
    return {'id': obs_id, 'instrument': instrument, 'time': time,
            'formula': formula_text, 'counts': counts, 'raw': raw}


def _certificate(record):
    _object(record, ('id', 'formula', 'members', 'lower', 'upper', 'max_range'))
    cert_id = ident(record['id'])
    formula_text = record['formula']
    counts = tuple(sorted(formula(formula_text).items()))
    members_raw = record['members']
    if not isinstance(members_raw, list) or not _MEMBER_MIN <= len(members_raw) <= _MEMBER_MAX:
        raise ValueError('certificate members must be an array of 1..64 entries')
    members = []
    seen = set()
    for member in members_raw:
        _object(member, ('observation_id', 'weight'))
        observation_id = ident(member['observation_id'])
        if observation_id in seen:
            raise ValueError('duplicate certificate member: ' + observation_id)
        seen.add(observation_id)
        weight = _amount(member['weight'], 'certificate weight')
        if weight <= 0:
            raise ValueError('certificate weight must be greater than zero')
        members.append({'observation_id': observation_id, 'weight': weight})
    lower = _amount(record['lower'], 'certificate lower')
    upper = _amount(record['upper'], 'certificate upper')
    if lower > upper:
        raise ValueError('certificate lower exceeds upper')
    max_range = _amount(record['max_range'], 'certificate max_range')
    return {'id': cert_id, 'formula': formula_text, 'counts': counts,
            'members': members, 'lower': lower, 'upper': upper,
            'max_range': max_range}


_RECORD_PARSERS = {
    'calibration': _calibration,
    'observation': _observation,
    'certificate': _certificate,
}


def _change(raw):
    if not isinstance(raw, dict):
        raise ValueError('change must be an object')
    action = raw.get('action')
    entity = raw.get('entity')
    if entity not in _ENTITIES:
        raise ValueError('change entity must be calibration, observation or certificate')
    if action == 'put':
        _object(raw, ('entity', 'action', 'record'))
        record = _RECORD_PARSERS[entity](raw['record'])
        return {'entity': entity, 'action': 'put', 'id': record['id'],
                'record': record}
    if action == 'delete':
        _object(raw, ('entity', 'action', 'id'))
        return {'entity': entity, 'action': 'delete',
                'id': ident(raw['id'])}
    raise ValueError('change action must be put or delete')


def _parse(data):
    """Validate the whole file, including trailing records, before accounting."""
    _object(data, ('batches', 'queries'))
    raw_batches = data['batches']
    raw_queries = data['queries']
    if not isinstance(raw_batches, list) or len(raw_batches) > _BATCH_LIMIT:
        raise ValueError('batches must be an array of at most 300 entries')
    if not isinstance(raw_queries, list) or len(raw_queries) > _QUERY_LIMIT:
        raise ValueError('queries must be an array of at most 60 entries')

    batches = []
    batch_ids = set()
    change_total = 0
    calibration_nodes = 0
    certificate_members = 0
    for raw_batch in raw_batches:
        _object(raw_batch, ('batch_id', 'changes'))
        batch_id = ident(raw_batch['batch_id'])
        if batch_id in batch_ids:
            raise ValueError('duplicate batch_id: ' + batch_id)
        batch_ids.add(batch_id)
        raw_changes = raw_batch['changes']
        if not isinstance(raw_changes, list):
            raise ValueError('changes must be an array')
        change_total += len(raw_changes)
        if change_total > _CHANGE_LIMIT:
            raise ValueError('too many changes across the file')
        changes = []
        touched = set()
        for raw_change in raw_changes:
            change = _change(raw_change)
            key = (change['entity'], change['id'])
            if key in touched:
                raise ValueError('duplicate change key in one batch')
            touched.add(key)
            if change['action'] == 'put':
                if change['entity'] == 'calibration':
                    calibration_nodes += len(change['record']['points'])
                    if calibration_nodes > _CALIBRATION_NODE_LIMIT:
                        raise ValueError('too many calibration nodes across the file')
                elif change['entity'] == 'certificate':
                    certificate_members += len(change['record']['members'])
                    if certificate_members > _CERTIFICATE_MEMBER_LIMIT:
                        raise ValueError('too many certificate members across the file')
            changes.append(change)
        batches.append({'batch_id': batch_id, 'changes': changes})

    queries = []
    query_ids = set()
    for raw_query in raw_queries:
        _object(raw_query, ('query_id', 'as_of'))
        query_id = ident(raw_query['query_id'])
        if query_id in query_ids:
            raise ValueError('duplicate query_id: ' + query_id)
        query_ids.add(query_id)
        as_of = _integer(raw_query['as_of'], 0, len(batches), 'query.as_of')
        queries.append({'query_id': query_id, 'as_of': as_of})
    return batches, queries


def _empty_state():
    return {'calibration': {}, 'observation': {}, 'certificate': {}}


def _snapshot(state):
    """Freeze a prefix view.

    Records are never mutated after creation, so the snapshot only keeps
    references inside tuples; retaining up to 61 queried prefixes costs O(61)
    pointers per surviving record rather than 61 whole hash tables.
    """
    return {entity: tuple(table.items()) for entity, table in state.items()}


def _index_snapshot(snapshot):
    calibrations = dict(snapshot['calibration'])
    by_instrument = _index_calibrations(calibrations)
    return dict(snapshot['observation']), dict(snapshot['certificate']), by_instrument


def _index_calibrations(calibrations):
    """Group active calibrations per instrument, time-sorted for overlap and
    coverage checks.  Overlapping intervals never survive a commit, so at most
    one indexed interval covers any given time."""
    by_instrument = {}
    for cal in calibrations.values():
        by_instrument.setdefault(cal['instrument'], []).append(cal)
    for entries in by_instrument.values():
        entries.sort(key=lambda item: (item['start'], item['end'], item['id']))
    return by_instrument


def _interpolate(cal, raw):
    points = cal['points']
    raws = [node['raw'] for node in points]
    if raw < raws[0] or raw > raws[-1]:
        return None
    index = bisect_right(raws, raw) - 1
    if index == len(points) - 1:
        return points[-1]['value']
    x0, y0 = points[index]['raw'], points[index]['value']
    x1, y1 = points[index + 1]['raw'], points[index + 1]['value']
    return y0 + (raw - x0) * (y1 - y0) / (x1 - x0)


def evaluate(by_instrument, observation):
    """Resolve one observation against an already-indexed calibration view."""
    for cal in by_instrument.get(observation['instrument'], ()):
        if cal['start'] <= observation['time'] < cal['end']:
            value = _interpolate(cal, observation['raw'])
            if value is None:
                return {'status': 'unavailable', 'calibration_id': cal['id'],
                        'value': None, 'reason': 'outside_curve'}
            return {'status': 'available', 'calibration_id': cal['id'],
                    'value': value, 'reason': None}
    return {'status': 'unavailable', 'calibration_id': None,
            'value': None, 'reason': 'no_calibration'}


def _overlap_reason(calibrations):
    for entries in _index_calibrations(calibrations).values():
        previous_end = None
        for cal in entries:
            if previous_end is not None and cal['start'] < previous_end:
                return 'overlap'
            previous_end = cal['end']
    return None


def _certificate_reason(state, by_instrument):
    """First failing reason across all active certificates, in the fixed
    global priority.  Per-certificate checks follow the same order so a
    certificate failing several rules reports its highest-priority one."""
    worst = len(_REASON_PRIORITY)
    for cert in state['certificate'].values():
        evaluations = {}
        reason = None
        for member in cert['members']:
            observation = state['observation'].get(member['observation_id'])
            if observation is None:
                reason = 'missing_member'
                break
        if reason is None:
            for member in cert['members']:
                observation = state['observation'][member['observation_id']]
                if observation['counts'] != cert['counts']:
                    reason = 'formula_mismatch'
                    break
        if reason is None:
            for member in cert['members']:
                observation = state['observation'][member['observation_id']]
                result = evaluate(by_instrument, observation)
                evaluations[member['observation_id']] = result
                if result['status'] != 'available':
                    reason = 'unavailable_member'
                    break
        if reason is None:
            values = [evaluations[m['observation_id']]['value']
                      for m in cert['members']]
            total_weight = sum(member['weight'] for member in cert['members'])
            mean = sum(member['weight']
                       * evaluations[member['observation_id']]['value']
                       for member in cert['members']) / total_weight
            if not cert['lower'] <= mean <= cert['upper']:
                reason = 'mean_out_of_range'
            elif max(values) - min(values) > cert['max_range']:
                reason = 'range_exceeded'
        if reason is not None:
            worst = min(worst, _REASON_PRIORITY.index(reason))
    return None if worst == len(_REASON_PRIORITY) else _REASON_PRIORITY[worst]


def _apply(state, changes):
    """Apply all changes to a private candidate and validate the final view."""
    candidate = {entity: dict(table) for entity, table in state.items()}
    for change in changes:
        table = candidate[change['entity']]
        if change['action'] == 'delete':
            if change['id'] not in table:
                return 'missing_delete', None
            del table[change['id']]
        else:
            table[change['id']] = change['record']
    reason = _overlap_reason(candidate['calibration'])
    if reason is not None:
        return reason, None
    by_instrument = _index_calibrations(candidate['calibration'])
    reason = _certificate_reason(candidate, by_instrument)
    if reason is not None:
        return reason, None
    return None, candidate


def _observation_report(by_instrument, observation):
    result = evaluate(by_instrument, observation)
    return {
        'id': observation['id'],
        'instrument': observation['instrument'],
        'time': observation['time'],
        'formula': observation['formula'],
        'raw': str(observation['raw']),
        'status': result['status'],
        'calibration_id': result['calibration_id'],
        'value': None if result['value'] is None else str(result['value']),
        'reason': result['reason'],
    }


def _certificate_report(by_instrument, observations, cert):
    members = []
    values = []
    total_weight = Fraction(0)
    weighted = Fraction(0)
    for member in sorted(cert['members'], key=lambda item: item['observation_id']):
        observation = observations[member['observation_id']]
        result = evaluate(by_instrument, observation)
        value = result['value']
        values.append(value)
        total_weight += member['weight']
        weighted += member['weight'] * value
        members.append({
            'observation_id': member['observation_id'],
            'weight': str(member['weight']),
            'value': str(value),
            'calibration_id': result['calibration_id'],
        })
    mean = weighted / total_weight
    return {
        'id': cert['id'],
        'formula': cert['formula'],
        'lower': str(cert['lower']),
        'upper': str(cert['upper']),
        'max_range': str(cert['max_range']),
        'mean': str(mean),
        'range': str(max(values) - min(values)),
        'members': members,
    }


def _query_report(snapshot, query):
    observations, certificates, by_instrument = _index_snapshot(snapshot)
    observation_reports = [
        _observation_report(by_instrument, observations[obs_id])
        for obs_id in sorted(observations)
    ]
    certificate_reports = [
        _certificate_report(by_instrument, observations, certificates[cert_id])
        for cert_id in sorted(certificates)
    ]
    return {'query_id': query['query_id'], 'as_of': query['as_of'],
            'observations': observation_reports,
            'certificates': certificate_reports}


def analyze(data):
    batches, queries = _parse(data)
    state = _empty_state()
    wanted = {query['as_of'] for query in queries}
    snapshots = {}
    if 0 in wanted:
        snapshots[0] = _snapshot(state)

    batch_reports = []
    for index, batch in enumerate(batches, start=1):
        reason, candidate = _apply(state, batch['changes'])
        if reason is None:
            state = candidate
            batch_reports.append({'batch_id': batch['batch_id'],
                                  'status': 'accepted', 'reason': None})
        else:
            batch_reports.append({'batch_id': batch['batch_id'],
                                  'status': 'rejected', 'reason': reason})
        if index in wanted:
            snapshots[index] = _snapshot(state)

    query_reports = [
        _query_report(snapshots[query['as_of']], query)
        for query in queries]
    return {'batches': batch_reports, 'queries': query_reports}
