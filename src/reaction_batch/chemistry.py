"""Bounded formula grammar and exact one-dimensional nullspace balance."""
import re,math,functools
from fractions import Fraction
from collections import Counter

def formula(text):
    if not isinstance(text,str) or not 1<=len(text)<=100:raise ValueError('invalid formula length')
    i=0
    def count():
        nonlocal i
        start=i
        while i<len(text) and text[i].isdigit():i+=1
        if start==i:return 1
        token=text[start:i]
        if token[0]=='0' or len(token)>4 or int(token)>1000:raise ValueError('invalid formula multiplier')
        return int(token)
    def group(depth=0):
        nonlocal i
        if depth>10:raise ValueError('formula nesting limit')
        result=Counter();terms=0
        while i<len(text) and text[i]!=')':
            if text[i]=='(':
                i+=1;part=group(depth+1)
                if i>=len(text) or text[i]!=')':raise ValueError('unclosed formula')
                i+=1
            else:
                m=re.match(r'[A-Z][a-z]?',text[i:])
                if not m:raise ValueError('expected element token')
                token=m.group();i+=len(token);part=Counter({token:1})
            n=count();terms+=1
            for e,c in part.items():
                result[e]+=c*n
                if result[e]>10**12:raise ValueError('atom count limit')
        if terms==0:raise ValueError('empty formula group')
        return result
    result=group()
    if i!=len(text):raise ValueError('unexpected closing bracket')
    return dict(result)

def balance(left,right):
    comps=[formula(f) for f in left+right];elements=sorted(set().union(*(c.keys() for c in comps)))
    cols=len(comps);a=[[Fraction(c.get(e,0)*(1 if j<len(left) else -1)) for j,c in enumerate(comps)] for e in elements]
    pivots=[];row=0
    for col in range(cols):
        hit=next((k for k in range(row,len(a)) if a[k][col]),None)
        if hit is None:continue
        a[row],a[hit]=a[hit],a[row];v=a[row][col];a[row]=[x/v for x in a[row]]
        for k in range(len(a)):
            if k!=row:
                v=a[k][col];a[k]=[x-v*y for x,y in zip(a[k],a[row])]
        pivots.append(col);row+=1
        if row==len(a):break
    free=[j for j in range(cols) if j not in pivots]
    if len(free)!=1:raise ValueError('balance requires one-dimensional nullspace')
    vector=[Fraction(0)]*cols;vector[free[0]]=1
    for k,col in enumerate(pivots):vector[col]=-a[k][free[0]]
    if all(v<0 for v in vector):vector=[-v for v in vector]
    if any(v<=0 for v in vector):raise ValueError('no strictly positive balance')
    denominator=math.lcm(*(v.denominator for v in vector));ints=[int(v*denominator) for v in vector];g=functools.reduce(math.gcd,ints)
    return [v//g for v in ints],comps

def quantity(value):
    if not isinstance(value,str) or not re.fullmatch(r'(?:0|[1-9][0-9]{0,8})(?:\.[0-9]{1,9})?',value):raise ValueError('moles must be bounded nonnegative decimal string')
    return Fraction(value)
