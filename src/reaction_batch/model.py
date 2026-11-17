from fractions import Fraction
from .common import require,rows
from .chemistry import balance,quantity

def analyze(data):
    require(data,'reactants products available_moles');left=rows(data['reactants'],20);right=rows(data['products'],20)
    if not left or not right or any(not isinstance(f,str) for f in left+right) or len(set(left+right))!=len(left+right):raise ValueError('unique nonempty chemical species required')
    coeff,comps=balance(left,right);available=data['available_moles']
    if not isinstance(available,dict) or set(available)!=set(left):raise ValueError('available moles required for exactly reactants')
    amounts=[quantity(available[f]) for f in left];extent=min(q/c for q,c in zip(amounts,coeff))
    consumed=[extent*c for c in coeff[:len(left)]];produced=[extent*c for c in coeff[len(left):]]
    return {'coefficients':dict(zip(left+right,coeff)),'extent_moles':str(extent),
            'limiting_reactants':[f for f,q,c in zip(left,amounts,coeff) if q/c==extent],
            'consumed_moles':dict(zip(left,map(str,consumed))),'remaining_moles':{f:str(q-c) for f,q,c in zip(left,amounts,consumed)},
            'produced_moles':dict(zip(right,map(str,produced))),
            'atom_totals_per_extent':{e:sum(c.get(e,0)*n for c,n in zip(comps[:len(left)],coeff)) for e in sorted(set().union(*(c.keys() for c in comps)))}}
