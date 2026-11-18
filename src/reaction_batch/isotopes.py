"""Exact molecular isotope joint distribution and selective channel accounting.

Every input atom is an independent draw from its element's user-supplied
isotope table.  A draw fixes both the nominal mass contribution and whether
the atom is tagged, so the joint law of ``(mass, tag count)`` is never the
product of the two marginal laws: assignments that differ only in which atom
position carries which isotope are distinct outcomes and their weights are
summed.

Only offsets ``S = M - B`` inside ``0..cutoff`` are kept as an explicit joint
state ``(S, T)``.  ``B`` is the base mass ``sum_e n_e * b_e`` where ``b_e`` is
the smallest mass number among *all* of element ``e``'s entries, including
zero-abundance ones.  Since every isotope's offset from ``b_e`` is
nonnegative, a partial offset above ``cutoff`` can never return; such states
are dropped from the DP and their probability and tag moments are recovered
exactly as ``total - kept`` (see ``_solve_query``).

Weights are accumulated as integers over the common denominator
``10 ** (9 * atoms)``: each abundance is a decimal string with at most nine
fraction digits, so its numerator over ``10**9`` is exact.  ``Fraction`` is
used only for the final reduced output values, keeping every arithmetic step
exact without repeated rational normalisation.
"""

import re
from fractions import Fraction

from .chemistry import formula, quantity
from .common import ident

_ELEMENT_LIMIT = 20
_ISOTOPE_LIMIT = 6
_QUERY_LIMIT = 4
_CHANNEL_LIMIT = 16
_ATOM_LIMIT = 80
_MASS_NUMBER_MAX = 300
_MASS_BOUND = 10 ** 9
_TAG_BOUND = 80
_CUTOFF_MAX = 128
_ABUNDANCE_SCALE = 10 ** 9
_ELEMENT_TOKEN = re.compile(r'[A-Z][a-z]?\Z')


def _object(value, fields):
    if not isinstance(value, dict) or set(value) != set(fields):
        raise ValueError('expected object fields: ' + ' '.join(fields))
    return value


def _array(value, minimum, maximum):
    if not isinstance(value, list) or not minimum <= len(value) <= maximum:
        raise ValueError('array length outside allowed bounds')
    return value


def _integer(value, minimum, maximum, label):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(label + ' must be an integer within its bounds')
    return value


def _fraction_text(value):
    return str(value)


def _parse(data):
    """Validate the whole batch; raise before any computation or publication."""
    _object(data, ('elements', 'queries'))
    elements = {}
    for raw in _array(data['elements'], 1, _ELEMENT_LIMIT):
        _object(raw, ('element', 'isotopes'))
        token = raw['element']
        if not isinstance(token, str) or _ELEMENT_TOKEN.fullmatch(token) is None:
            raise ValueError('invalid element token')
        if token in elements:
            raise ValueError('duplicate element: ' + token)
        isotopes = []
        seen_masses = set()
        abundance_total = Fraction(0)
        for item in _array(raw['isotopes'], 1, _ISOTOPE_LIMIT):
            _object(item, ('mass_number', 'abundance', 'tagged'))
            mass_number = _integer(item['mass_number'], 1, _MASS_NUMBER_MAX,
                                   'mass_number')
            if mass_number in seen_masses:
                raise ValueError('duplicate mass_number for ' + token)
            seen_masses.add(mass_number)
            abundance = quantity(item['abundance'])
            if abundance > 1:
                raise ValueError('isotope abundance exceeds one')
            tagged = item['tagged']
            if type(tagged) is not bool:
                raise ValueError('tagged must be a JSON boolean')
            abundance_total += abundance
            isotopes.append((mass_number, abundance, tagged))
        if abundance_total != 1:
            raise ValueError('isotope abundances must sum to one')
        elements[token] = isotopes

    queries = []
    seen_queries = set()
    for raw in _array(data['queries'], 1, _QUERY_LIMIT):
        _object(raw, ('query_id', 'formula', 'cutoff', 'channels'))
        query_id = ident(raw['query_id'])
        if query_id in seen_queries:
            raise ValueError('duplicate query_id: ' + query_id)
        seen_queries.add(query_id)
        counts = formula(raw['formula'])
        if sum(counts.values()) > _ATOM_LIMIT:
            raise ValueError('formula exceeds the atom limit')
        for token in counts:
            if token not in elements:
                raise ValueError('missing isotope table for ' + token)
        base_mass = sum(n * min(m for m, _a, _t in elements[token])
                        for token, n in counts.items())
        cutoff = _integer(raw['cutoff'], 0, _CUTOFF_MAX, 'cutoff')
        channels = []
        seen_channels = set()
        for item in _array(raw['channels'], 0, _CHANNEL_LIMIT):
            _object(item, ('channel_id', 'mass', 'tag_count'))
            channel_id = ident(item['channel_id'])
            if channel_id in seen_channels:
                raise ValueError('duplicate channel_id: ' + channel_id)
            seen_channels.add(channel_id)
            _object(item['mass'], ('lower', 'upper'))
            mass_lower = _integer(item['mass']['lower'], 0, _MASS_BOUND,
                                  'mass lower')
            mass_upper = _integer(item['mass']['upper'], 0, _MASS_BOUND,
                                  'mass upper')
            if mass_lower > mass_upper:
                raise ValueError('mass lower must not exceed upper')
            if mass_upper > base_mass + cutoff:
                raise ValueError('mass upper exceeds base mass plus cutoff')
            _object(item['tag_count'], ('lower', 'upper'))
            tag_lower = _integer(item['tag_count']['lower'], 0, _TAG_BOUND,
                                 'tag_count lower')
            tag_upper = _integer(item['tag_count']['upper'], 0, _TAG_BOUND,
                                 'tag_count upper')
            if tag_lower > tag_upper:
                raise ValueError('tag_count lower must not exceed upper')
            channels.append((channel_id, mass_lower, mass_upper,
                             tag_lower, tag_upper))
        queries.append({'query_id': query_id, 'counts': counts,
                        'cutoff': cutoff, 'channels': channels})
    return elements, queries


def _scaled_abundance(value):
    """Exact numerator of ``value`` over ``10**9`` (denominator divides it)."""
    return value.numerator * (_ABUNDANCE_SCALE // value.denominator)


def _joint(counts, elements, cutoff):
    """Joint weights for every retained ``(offset, tag count)`` state.

    Returns ``(states, denominator, base_mass)`` where ``states`` maps
    ``(s, t)`` to a nonnegative integer weight and the exact probability of a
    state is ``weight / denominator``.  Atom positions are distinct sampling
    events; equal ``(s, t)`` states accumulate all contributing weights.
    """
    base_mass = 0
    sequence = []
    for token, number in counts.items():
        rows = elements[token]
        base = min(mass for mass, _a, _t in rows)
        base_mass += number * base
        options = []
        for mass, abundance, tagged in rows:
            scaled = _scaled_abundance(abundance)
            if scaled:
                options.append((mass - base, 1 if tagged else 0, scaled))
        for _ in range(number):
            sequence.append(options)

    # suffix[i] is the minimum offset still obtainable from atoms i and later;
    # a state that cannot stay inside the cutoff even in the best case is
    # dropped early.  This is exact because offsets never decrease.
    suffix = [0] * (len(sequence) + 1)
    for index in range(len(sequence) - 1, -1, -1):
        suffix[index] = suffix[index + 1] + min(o[0] for o in sequence[index])

    states = {(0, 0): 1}
    for index, options in enumerate(sequence):
        remaining = suffix[index + 1]
        updated = {}
        for (offset, tags), weight in states.items():
            for delta, tag_delta, scaled in options:
                new_offset = offset + delta
                if new_offset > cutoff or new_offset + remaining > cutoff:
                    continue
                key = (new_offset, tags + tag_delta)
                updated[key] = updated.get(key, 0) + weight * scaled
        states = updated

    return states, _ABUNDANCE_SCALE ** len(sequence), base_mass


def _stats(probability, first_moment, second_factorial_moment):
    if probability == 0:
        return {'probability': '0', 'first_moment': '0',
                'second_factorial_moment': '0',
                'mean_tags': None, 'variance_tags': None}
    mean = first_moment / probability
    variance = (second_factorial_moment + first_moment) / probability - mean * mean
    return {'probability': _fraction_text(probability),
            'first_moment': _fraction_text(first_moment),
            'second_factorial_moment': _fraction_text(second_factorial_moment),
            'mean_tags': _fraction_text(mean),
            'variance_tags': _fraction_text(variance)}


def _solve_query(elements, query):
    counts = query['counts']
    cutoff = query['cutoff']
    channels = query['channels']

    states, denominator, base_mass = _joint(counts, elements, cutoff)

    # Untruncated total tag moments: each atom independently contributes 1 with
    # probability q_e (the tagged share of its element).
    total_first = Fraction(0)
    total_squares = Fraction(0)
    for token, number in counts.items():
        tagged_share = sum((abundance for _m, abundance, tagged
                            in elements[token] if tagged), Fraction(0))
        total_first += number * tagged_share
        total_squares += number * tagged_share * tagged_share
    total_second = total_first * total_first - total_squares

    envelope_weights = {}
    kept_count = kept_first = kept_second = 0
    channel_counts = [0] * len(channels)
    channel_first = [0] * len(channels)
    channel_second = [0] * len(channels)

    for (offset, tags), weight in states.items():
        envelope_weights[offset] = envelope_weights.get(offset, 0) + weight
        kept_count += weight
        kept_first += tags * weight
        kept_second += tags * (tags - 1) * weight
        mass = base_mass + offset
        for index, (_cid, mass_lower, mass_upper, tag_lower, tag_upper) \
                in enumerate(channels):
            if mass_lower <= mass <= mass_upper and \
                    tag_lower <= tags <= tag_upper:
                channel_counts[index] += weight
                channel_first[index] += tags * weight
                channel_second[index] += tags * (tags - 1) * weight
                break

    envelope = [
        {'mass': base_mass + offset,
         'probability': _fraction_text(Fraction(weight, denominator))}
        for offset, weight in sorted(envelope_weights.items()) if weight]

    channel_reports = []
    for index, channel in enumerate(channels):
        channel_reports.append({
            'channel_id': channel[0],
            'stats': _stats(Fraction(channel_counts[index], denominator),
                            Fraction(channel_first[index], denominator),
                            Fraction(channel_second[index], denominator)),
        })

    kept_probability = Fraction(kept_count, denominator)
    overflow = _stats(
        Fraction(1) - kept_probability,
        total_first - Fraction(kept_first, denominator),
        total_second - Fraction(kept_second, denominator))

    assigned_count = sum(channel_counts)
    assigned_first = sum(channel_first)
    assigned_second = sum(channel_second)
    missed = _stats(
        Fraction(denominator - assigned_count, denominator),
        total_first - Fraction(assigned_first, denominator),
        total_second - Fraction(assigned_second, denominator))

    total = _stats(Fraction(1), total_first, total_second)

    return {'query_id': query['query_id'],
            'base_mass': base_mass,
            'envelope': envelope,
            'overflow': overflow,
            'channels': channel_reports,
            'missed': missed,
            'total': total}


def analyze(data):
    elements, queries = _parse(data)
    return {'queries': [_solve_query(elements, query) for query in queries]}
