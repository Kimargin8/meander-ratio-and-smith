"""Exact two-orientation occupation tables and an independent small-order oracle."""
import argparse
from collections import Counter
from functools import cache
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path

def require(ok,msg):
    if not ok:
        raise ValueError(msg)

def compositions(k,r):
    if r==0:
        return int(k==0)
    return comb(k-1,r-1) if k>=r else 0

def from_series(n,counts):
    # Mixed binomial moments are coefficients of the marked-run substitution.
    moments={}
    for i in range(n):
        for j in range(n-i):
            m=n-i-j
            ways=sum(comb(m,a)*comb(m-a,b)*compositions(i,a)*compositions(j,b)
                     for a in range(min(m,i)+1) for b in range(min(m-a,j)+1))
            moments[i,j]=ways*counts[m]
    # Two triangular inversions allow arbitrary integer precision without interpolation.
    undo_plus={(i,j):sum((-1)**(k-i)*comb(k,i)*moments[k,j] for k in range(i,n-j))
               for i,j in moments}
    exact={(i,j):sum((-1)**(k-j)*comb(k,j)*undo_plus[i,k] for k in range(j,n-i))
           for i,j in moments}
    require(all(v>=0 and v==exact[j,i] for (i,j),v in exact.items()),"invalid signs or asymmetry")
    require(sum(exact.values())==counts[n],"total meanders")
    if n>1:
        require(sum(i*v for (i,j),v in exact.items())==(n-1)*counts[n-1],"one-site marginal")
    return exact

def matchings(n):
    for opening in combinations(range(2*n),n):
        start=set(opening);stack=[];pair=[0]*(2*n)
        for i in range(2*n):
            if i in start:
                stack.append(i)
            elif stack:
                j=stack.pop();pair[i],pair[j]=j,i
            else:
                break
        else:
            yield pair

def oracle(n):
    catalogue=[]
    for p in matchings(n):
        first=sum(1<<i for i in range(n-1) if p[2*i+1]==2*i)
        last=sum(1<<i for i in range(n-1) if p[2*i+1]==2*i+2)
        catalogue.append((p,first,last))
    ans=Counter()
    for p,a,b in catalogue:
        for q,c,d in catalogue:
            j=q[p[0]];length=1
            while j:
                length+=1;j=q[p[j]]
            if length==n:
                ans[(a&d).bit_count(),(b&c).bit_count()]+=1
    return ans

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input",type=Path);ap.add_argument("output",type=Path)
    args=ap.parse_args()
    data=args.input.read_bytes()
    counts=dict(tuple(map(int,line.split())) for line in data.decode("utf-8-sig").splitlines()
                if line.strip() and not line.startswith("#"))
    rows={}
    for n in range(1,max(counts)+1):
        h=from_series(n,counts)
        if n<=8:
            require(oracle(n)==Counter({k:v for k,v in h.items() if v}),"direct enumeration")
        rows[n]=[[h[i,j] for j in range(n-i)] for i in range(n)]
    result={"indices":"rows i = positive count, columns j = negative count, i+j <= n-1",
            "direct_pairs_through":8,"input_sha256":sha256(data).hexdigest(),"orders":rows}
    args.output.write_text(json.dumps(result,separators=(",",":"))+"\n",encoding="ascii")
    print("Exact joint tables to order",max(counts),"direct ordered-pair enumeration to order 8.")

if __name__=="__main__":
    main()
