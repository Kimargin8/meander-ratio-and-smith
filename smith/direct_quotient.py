"""Check the direct, division-free formula for V against the archived multiplier L.

For each pair of matchings, trace the physical loops once and look for the
first repeated loop among the chosen crossing arcs. Compare exact coefficient
lists, using homogeneity to check the two-variable identity.
"""
import argparse
from collections import Counter
from functools import cache
from hashlib import sha256
import json
from pathlib import Path
import sys


def product(a,b):
    coeff=[0]*(len(a)+len(b)-1) if a and b else []
    for i,v in enumerate(a):
        if v:
            for j,w in enumerate(b):
                coeff[i+j]+=v*w
    return tuple(coeff)


def selected_arcs(p):
    h=peak=cut=0
    for i,v in enumerate(p):
        h+=1 if v>i else -1
        if h>peak: peak,cut=h,i+1
    return tuple(i for i in range(cut) if p[i]>=cut)


def loop_labels(p,q):
    labels=[-1]*len(p)
    count=0
    for i in range(len(p)):
        if labels[i]>=0: continue
        j=i
        while labels[j]<0:
            v=p[j]
            labels[j]=labels[v]=count
            j=q[v]
        count+=1
    return count,labels


def cutoff(arcs,labels):
    used=0
    for j,endpoint in enumerate(arcs):
        bit=1<<labels[endpoint]
        if used&bit: return j
        used|=bit
    return len(arcs)


def run(path,weight_file):
    # Shared code supplies only the independent Chebyshev diagonal and join
    # oracle. The topology formula below uses neither polynomial division nor
    # the projector generator to obtain V.
    from verify_rows import diagonal,parts,joins
    raw=json.loads(path.read_bytes())
    patterns=raw["matching_order"]; n=raw["order"]; size=len(patterns)
    weights={v["N"]:{int(r):a for r,a in v["weights"].items()}
             for v in json.loads(weight_file.read_bytes())["certificates"]}
    arcs=[selected_arcs(p) for p in patterns]
    partitions=[];complements=[]
    for p in patterns:
        perm=[p[2*i+1]//2 for i in range(n)]
        partitions.append(parts(perm))
        complements.append(parts([perm.index((i+1)%n) for i in range(n)]))
    colour=[len(blocks) for blocks in partitions]
    diag=diagonal(n)
    gram=[bytearray(size) for _ in patterns]; cap=[bytearray(size) for _ in patterns]
    for i,p in enumerate(patterns):
        for j in range(i+1):
            c,labels=loop_labels(p,patterns[j])
            gram[i][j]=gram[j][i]=c
            cap[i][j]=cutoff(arcs[i],labels);cap[j][i]=cutoff(arcs[j],labels)
            a=joins(partitions[i],partitions[j]);b=joins(complements[i],complements[j])
            if a+b!=c+1 or a-b!=colour[i]+colour[j]-n-1:
                raise ValueError(("join counts",i,j))

    @cache
    def direct(k,d,a,b):
        ans=Counter()
        for r, coeff in weights[k+1].items():
            if r>d+1 or not coeff: continue
            x=a-r//2; y=b-r//2
            if min(x,y)<0: raise ValueError("negative monomial exponent")
            for power,c in enumerate(coeff):
                if c: ans[x+power,y+power]+=c
        return tuple((i,j,v) for (i,j),v in sorted(ans.items()) if v)

    @cache
    def expanded(k,terms):
        if not terms: return ()
        # All terms have fixed degree difference, hence no two distinct
        # bivariate monomials coalesce at X=Y=x.
        if len({a-b for a,b,c in terms})!=1: raise ValueError("mixed grading")
        coefficients=[0]*(1+max(a+b for a,b,c in terms))
        for a,b,c in terms: coefficients[a+b]=c
        return product(coefficients,diag[k])

    actual_types=set();one_types=set();nonzero=0;forced=0
    quotient_hash=sha256()
    @cache
    def wire(terms):
        return (json.dumps(terms,separators=(",",":"))+"\n").encode("ascii")
    for i,left in enumerate(raw["left_rows"]):
        k=len(arcs[i]); acc=[[] for _ in range(size)]
        for j,terms in left:
            delta=colour[i]-colour[j]
            if any(l<abs(delta) or (l-delta)%2 for l,c in terms):
                raise ValueError("left multiplier has no homogeneous lift")
            for col,loops in enumerate(gram[j]):
                v=acc[col]
                for power,coeff in terms:
                    deg=power+loops
                    if len(v)<=deg: v.extend([0]*(deg+1-len(v)))
                    v[deg]+=coeff
        for j,coeff in enumerate(acc):
            while coeff and coeff[-1]==0: coeff.pop()
            delta=colour[i]+colour[j]-n-1
            a=(gram[i][j]+1+delta)//2; b=gram[i][j]+1-a
            terms=direct(k,cap[i][j],a,b)
            if tuple(coeff)!=expanded(k,terms): raise ValueError(("LG != DV",n,i,j))
            quotient_hash.update(wire(terms))
            if terms:
                nonzero+=1;actual_types.add(terms);one_types.add(tuple((x+y,v) for x,y,v in terms))
            if cap[i][j]+1<max(2,(k+1)//2+1): forced+=1
    return {"order":n,"dimension":size,"entries_per_ring":size**2,"nonzero":nonzero,
            "zeros_forced_by_cutoff":forced,"distinct_nonzero_bivariate_entries":len(actual_types),
            "distinct_nonzero_univariate_entries":len(one_types),"bivariate_row_major_sha256":quotient_hash.hexdigest(),
            "left_file_sha256":sha256(path.read_bytes()).hexdigest()}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("weights",type=Path)
    parser.add_argument("rows",type=Path,nargs="+")
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    sys.path.insert(0,str(args.weights.resolve().parent))
    results=[]
    for file in args.rows:
        result=run(file,args.weights);results.append(result)
        print(json.dumps(result),flush=True)
    if args.output: args.output.write_text(json.dumps(results,indent=2)+"\n",encoding="ascii")
