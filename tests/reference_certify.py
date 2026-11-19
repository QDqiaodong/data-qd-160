"""Independent, deliberately small reference implementation for reaction-certify.

This oracle is kept out of the production package and never imports
``reaction_batch.certify`` or any production calibration/certificate
evaluation function.  It re-derives every prefix view directly from the empty
state by replaying batches one by one on plain dictionaries, re-checks every
active certificate on every candidate, and renders observation/certificate
reports with its own interpolation and statistics code.  Only the baseline
grammar helpers (``chemistry.formula`` / ``chemistry.quantity``) are shared,
because those define the input literals rather than the certification rules.

It is intentionally O(batches^2) and suitable only for the small generated
sequences used by ``tests/test_certify.py``; the production implementation
keeps prefix snapshots instead.
"""

from bisect import bisect_right
from fractions import Fraction

from reaction_batch.chemistry import formula, quantity

_PRIORITY = ('missing_delete', 'overlap', 'missing_member',
             'formula_mismatch', 'unavailable_member',
             'mean_out_of_range', 'range_exceeded')


def _fraction(text):
    return quantity(text)


def _counts(text):
    return tuple(sorted(formula(text).items()))


def _record(change):
    record = change['record']
    kind = change['entity']
    if kind == 'calibration':
        return {'id': record['id'], 'instrument': record['instrument'],
                'start': record['start'], 'end': record['end'],
                'points': [{'raw': _fraction(p['raw']),
                            'value': _fraction(p['value'])}
                           for p in record['points']]}
    if kind == 'observation':
        return {'id': record['id'], 'instrument': record['instrument'],
                'time': record['time'], 'formula': record['formula'],
                'counts': _counts(record['formula']),
                'raw': _fraction(record['raw'])}
    return {'id': record['id'], 'formula': record['formula'],
            'counts': _counts(record['formula']),
            'members': [{'observation_id': m['observation_id'],
                         'weight': _fraction(m['weight'])}
                        for m in record['members']],
            'lower': _fraction(record['lower']),
            'upper': _fraction(record['upper']),
            'max_range': _fraction(record['max_range'])}


def _interpolate(calibration, raw):
    raws = [point['raw'] for point in calibration['points']]
    if raw < raws[0] or raw > raws[-1]:
        return None
    position = bisect_right(raws, raw) - 1
    if position == len(calibration['points']) - 1:
        return calibration['points'][-1]['value']
    left = calibration['points'][position]
    right = calibration['points'][position + 1]
    return (left['value']
            + (raw - left['raw']) * (right['value'] - left['value'])
            / (right['raw'] - left['raw']))


def _catalog(calibrations):
    grouped = {}
    for cal in calibrations.values():
        grouped.setdefault(cal['instrument'], []).append(cal)
    for entries in grouped.values():
        entries.sort(key=lambda item: (item['start'], item['end'], item['id']))
    return grouped


def _evaluate(grouped, observation):
    for cal in grouped.get(observation['instrument'], ()):
        if cal['start'] <= observation['time'] < cal['end']:
            value = _interpolate(cal, observation['raw'])
            if value is None:
                return ('unavailable', cal['id'], None, 'outside_curve')
            return ('available', cal['id'], value, None)
    return ('unavailable', None, None, 'no_calibration')


def _overlap(calibrations):
    for entries in _catalog(calibrations).values():
        cursor = None
        for cal in entries:
            if cursor is not None and cal['start'] < cursor:
                return True
            cursor = cal['end']
    return False


def _first_certificate_failure(state, grouped):
    rank = len(_PRIORITY)
    for cert in state['certificate'].values():
        reason = None
        for member in cert['members']:
            if member['observation_id'] not in state['observation']:
                reason = 'missing_member'
                break
        if reason is None:
            for member in cert['members']:
                obs = state['observation'][member['observation_id']]
                if obs['counts'] != cert['counts']:
                    reason = 'formula_mismatch'
                    break
        evaluated = {}
        if reason is None:
            for member in cert['members']:
                obs = state['observation'][member['observation_id']]
                evaluated[member['observation_id']] = _evaluate(grouped, obs)
                if evaluated[member['observation_id']][0] != 'available':
                    reason = 'unavailable_member'
                    break
        if reason is None:
            values = [evaluated[m['observation_id']][2]
                      for m in cert['members']]
            weight_total = sum(m['weight'] for m in cert['members'])
            mean = sum(m['weight'] * evaluated[m['observation_id']][2]
                       for m in cert['members']) / weight_total
            if not cert['lower'] <= mean <= cert['upper']:
                reason = 'mean_out_of_range'
            elif max(values) - min(values) > cert['max_range']:
                reason = 'range_exceeded'
        if reason is not None:
            rank = min(rank, _PRIORITY.index(reason))
    return None if rank == len(_PRIORITY) else _PRIORITY[rank]


def _candidate(state, changes):
    trial = {kind: dict(table) for kind, table in state.items()}
    for change in changes:
        table = trial[change['entity']]
        if change['action'] == 'delete':
            if change['id'] not in table:
                return 'missing_delete', None
            del table[change['id']]
        else:
            table[change['record']['id']] = _record(change)
    if _overlap(trial['calibration']):
        return 'overlap', None
    grouped = _catalog(trial['calibration'])
    reason = _first_certificate_failure(trial, grouped)
    if reason is not None:
        return reason, None
    return None, trial


def _observation_item(grouped, obs):
    status, cal_id, value, reason = _evaluate(grouped, obs)
    return {'id': obs['id'], 'instrument': obs['instrument'],
            'time': obs['time'], 'formula': obs['formula'],
            'raw': str(obs['raw']), 'status': status,
            'calibration_id': cal_id,
            'value': None if value is None else str(value),
            'reason': reason}


def _certificate_item(grouped, state, cert):
    members = []
    values = []
    weight_total = Fraction(0)
    weighted = Fraction(0)
    for member in sorted(cert['members'],
                         key=lambda item: item['observation_id']):
        obs = state['observation'][member['observation_id']]
        _, cal_id, value, _ = _evaluate(grouped, obs)
        values.append(value)
        weight_total += member['weight']
        weighted += member['weight'] * value
        members.append({'observation_id': member['observation_id'],
                        'weight': str(member['weight']),
                        'value': str(value), 'calibration_id': cal_id})
    mean = weighted / weight_total
    return {'id': cert['id'], 'formula': cert['formula'],
            'lower': str(cert['lower']), 'upper': str(cert['upper']),
            'max_range': str(cert['max_range']), 'mean': str(mean),
            'range': str(max(values) - min(values)), 'members': members}


def replay(data):
    """Return per-batch statuses and every prefix view's rendered queries."""
    state = {'calibration': {}, 'observation': {}, 'certificate': {}}
    statuses = []
    prefix_views = [_render(state)]
    for batch in data['batches']:
        reason, candidate = _candidate(state, batch['changes'])
        if reason is None:
            state = candidate
            statuses.append({'batch_id': batch['batch_id'],
                             'status': 'accepted', 'reason': None})
        else:
            statuses.append({'batch_id': batch['batch_id'],
                             'status': 'rejected', 'reason': reason})
        prefix_views.append(_render(state))
    queries = [
        {'query_id': query['query_id'], 'as_of': query['as_of'],
         'observations': prefix_views[query['as_of']]['observations'],
         'certificates': prefix_views[query['as_of']]['certificates']}
        for query in data['queries']
    ]
    return statuses, queries


def _render(state):
    grouped = _catalog(state['calibration'])
    return {
        'observations': [_observation_item(grouped, state['observation'][key])
                         for key in sorted(state['observation'])],
        'certificates': [_certificate_item(grouped, state,
                                          state['certificate'][key])
                         for key in sorted(state['certificate'])],
    }
