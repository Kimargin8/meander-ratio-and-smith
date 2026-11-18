"""Integral polynomial row certificates by inserting rational Jones--Wenzl projectors.

The single-variable entries encode both rings. The independent verifier
checks when the homogeneous lift belongs to Z[X,Y]. Requires SymPy 1.14.
"""
import argparse
import json
from collections import defaultdict
from functools import cache
from pathlib import Path

from sympy.polys import ring, ZZ
from check_bezout import beraha_factors, certificate

R, x = ring("x", ZZ)

@cache
def multiply(p, q):
    """Stack two diagrams whose lower rows appear in reverse order on the boundary."""
    m = len(p)//2
    parent = list(range(3*m))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for h, offset in ((p,0),(q,m)):
        def node(i):
            return offset + (i if i<m else 3*m-1-i)
        for i,j in enumerate(h):
            if i<j:
                parent[find(node(i))]=find(node(j))
    outer=[find(i) for i in list(range(m))+list(range(3*m-1,2*m-1,-1))]
    indices=defaultdict(list)
    for i,root in enumerate(outer):
        indices[root].append(i)
    out=[None]*(2*m)
    for i,j in indices.values():
        out[i],out[j]=j,i
    loops=len({find(i) for i in range(m,2*m)}-set(outer))
    return tuple(out),loops

def projectors(n):
    quantum=[R.zero,R.one]
    for j in range(2,n+2):
        quantum.append(x*quantum[-1]-quantum[-2])
    nums={(1,0):R.one}; denominator=R.one
    all_projectors={1:(nums,denominator)}
    for m in range(2,n+1):
        embedded={}
        for p,v in nums.items():
            a=[None]*(2*m)
            shift=lambda i:i if i<m-1 else i+2
            for i,j in enumerate(p):
                a[shift(i)]=shift(j)
            a[m-1],a[m]=m,m-1
            embedded[tuple(a)]=v
        e=list(reversed(range(2*m)))
        e[m-2],e[m-1],e[m],e[m+1]=m-1,m-2,m+1,m
        tmp=defaultdict(lambda:R.zero)
        for p,v in embedded.items():
            diagram,loops=multiply(p,tuple(e))
            tmp[diagram]+=v*x**loops
        nums=defaultdict(lambda:R.zero,{p:v*quantum[m]*denominator for p,v in embedded.items()})
        for p,v in tmp.items():
            for q,w in embedded.items():
                diagram,loops=multiply(p,q)
                nums[diagram]-=quantum[m-1]*v*w*x**loops
        denominator=quantum[m]*denominator**2
        gcd=denominator
        for v in nums.values():
            gcd=gcd.gcd(v)
        denominator=denominator.exquo(gcd)
        nums={p:v.exquo(gcd) for p,v in nums.items() if v}
        assert nums[tuple(reversed(range(2*m)))]==denominator
        all_projectors[m]=(nums,denominator)
    return all_projectors

@cache
def diagrams(n):
    if n==0:
        return [()]
    return [(2*i+1,)+tuple(1+j for j in p)+(0,)+tuple(2*i+2+j for j in q)
            for i in range(n) for p in diagrams(i) for q in diagrams(n-i-1)]

def shape(p):
    heights=[];height=0
    for i,j in enumerate(p):
        height+=1 if j>i else -1
        heights.append(height)
    return sum(heights),max(heights),1+heights.index(max(heights))

def as_x_squared(coeff):
    return R.from_dict({(2*i,):ZZ(v) for i,v in enumerate(coeff) if v})

def encode(n):
    project=projectors(n)
    pats=sorted(diagrams(n),key=lambda p:(shape(p)[0],p))
    ids={p:i for i,p in enumerate(pats)}
    factors=beraha_factors(n+1)
    forms={}
    for k in range(1,n+1):
        weights=certificate(k+1,factors)["weights"]
        forms[k]=[]
        for sr,weight in weights.items():
            r=int(sr)
            if not weight:
                continue
            prefactor=as_x_squared(weight)
            for d in range(2,k+2):
                if r%d:
                    prefactor*=as_x_squared(factors[d])
            nums,den=project[r-1]
            forms[k].append((r-1,prefactor.exquo(den),nums))
    rows=[]
    for i,p in enumerate(pats):
        area,h,cut=shape(p)
        across=[j for j in range(cut) if p[j]>=cut]
        row=defaultdict(lambda:R.zero)
        for size,prefactor,nums in forms[h]:
            left=across[:size]
            ends=left+sorted(p[j] for j in left)
            for diagram,coefficient in nums.items():
                target=list(p)
                for a,b in enumerate(diagram):
                    target[ends[a]]=ends[b]
                j=ids[tuple(target)]
                assert i==j or shape(pats[j])[0]<area
                row[j]+=prefactor*coefficient
        assert row[i]==R.one
        rows.append([[j,[[int(d[0]),int(a)] for d,a in sorted(v.items())]] for j,v in sorted(row.items()) if v])
    return {"order":n,"matching_order":pats,"left_rows":rows}

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("order",type=int,choices=range(1,9))
    ap.add_argument("output",type=Path)
    args=ap.parse_args()
    data=encode(args.order)
    # The archived row fingerprints include a final CRLF on every platform.
    encoded = json.dumps(data, separators=(",", ":")).encode("ascii")
    args.output.write_bytes(encoded + b"\r\n")
    print(args.order,len(data["matching_order"]),args.output.stat().st_size,flush=True)
