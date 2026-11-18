"""Verify both integral polynomial identities from sparse homogeneous rows.

Uses the Python standard library. No Jones--Wenzl coefficients, factor
weights, external packages, numeric evaluations or random primes are used.
"""
import argparse
import hashlib
import json
from math import comb
from pathlib import Path

def require(condition, message):
    if not condition:
        raise ValueError(message)

def trim(a):
    while a and a[-1]==0:
        a.pop()
    return a

def divmonic(a,b):
    require(b[-1]==1,"nonmonic divisor")
    a=trim(a.copy())
    q=[0]*max(0,1+len(a)-len(b))
    while len(a)>=len(b):
        shift=len(a)-len(b); v=a[-1]; q[shift]=v
        for j,w in enumerate(b):
            a[shift+j]-=v*w
        trim(a)
    require(not a,"nonzero polynomial remainder")
    return q

def diagonal(n):
    prev=[]; current=[1]; prim={}; target=[[1]]
    for d in range(2,n+2):
        val=[0]+current
        for i,v in enumerate(prev):
            val[i]-=v
        prev,current=current,val
        factor=val.copy()
        for j in range(2,d):
            if d%j==0:
                factor=divmonic(factor,prim[j])
        prim[d]=factor
        v=[0]*(len(factor)+len(target[-1])-1)
        for i,a in enumerate(target[-1]):
            for j,b in enumerate(factor):
                v[i+j]+=a*b
        target.append(v)
    return target

def parts(perm):
    rem=(1<<len(perm))-1; blocks=[]
    while rem:
        mask=0; i=(rem&-rem).bit_length()-1
        while not mask>>i&1:
            mask|=1<<i; i=perm[i]
        rem^=mask;blocks.append(mask)
    return blocks

def joins(blocks_a, blocks_b):
    result=blocks_a.copy()
    for block in blocks_b:
        keep=[]
        for mask in result:
            if mask&block:
                block|=mask
            else:
                keep.append(mask)
        keep.append(block)
        result=keep
    return len(result)

def loops(p,q):
    seen=0; cycles=0
    for i in range(len(p)):
        if (seen>>i)&1:
            continue
        cycles+=1; j=i
        while not (seen>>j)&1:
            v=p[j]
            seen|=(1<<j)|(1<<v)
            j=q[v]
    return cycles

def check(file):
    src=file.read_bytes(); obj=json.loads(src); n=obj["order"]
    patterns=[tuple(p) for p in obj["matching_order"]]
    size=comb(2*n,n)//(n+1)
    require(len(patterns)==len(set(patterns))==size,"catalogue size")
    heights=[];areas=[];black=[];white=[]
    for p in patterns:
        require(len(p)==2*n and sorted(p)==list(range(2*n)),"invalid matching")
        stack=[]; h=area=0
        for i,j in enumerate(p):
            require(i!=j and p[j]==i,"not a perfect matching")
            if i<j:
                stack.append(j)
            else:
                require(stack.pop()==i,"crossing")
            h=max(h,len(stack));area+=len(stack)
        heights.append(h);areas.append(area)
        perm=[p[2*i+1]//2 for i in range(n)]
        black.append(parts(perm));white.append(parts([perm.index((i+1)%n) for i in range(n)]))
    diag=diagonal(n)
    require(sum(len(diag[h])-1 for h in heights)==n*size,"degree budget")
    colours=[len(b) for b in black]
    G=[bytearray(size) for _ in patterns]
    for i,p in enumerate(patterns):
        for j in range(i+1):
            c=loops(p,patterns[j]);a=joins(black[i],black[j]);b=joins(white[i],white[j])
            require(a+b==c+1 and a-b==colours[i]+colours[j]-n-1,"join/loop identification")
            require((c==n)==(i==j),"unexpected maximum loop count")
            G[i][j]=G[j][i]=c
    require(len(obj["left_rows"])==size,"row count")
    sparse_count=0; nonzero_v=0;max_deg=0
    for i,sparse_row in enumerate(obj["left_rows"]):
        seen=set(); acc=[[] for _ in range(size)]
        for j,terms in sparse_row:
            require(0<=j<=i and j not in seen and terms,"row index");seen.add(j)
            require((terms==[[0,1]]) if i==j else areas[j]<areas[i],"unitriangularity")
            sparse_count+=1
            delta=colours[i]-colours[j]
            require(len({k for k,v in terms})==len(terms),"duplicate exponent")
            for k,v in terms:
                require(type(k)==int and type(v)==int and v!=0,"coefficient domain")
                require(k>=abs(delta) and (k-delta)%2==0,"left multiplier has no integral two-colour lift")
            for col,c in enumerate(G[j]):
                a=acc[col]
                for k,v in terms:
                    d=k+c
                    if len(a)<=d:
                        a.extend([0]*(d+1-len(a)))
                    a[d]+=v
        require(i in seen,"missing unit pivot")
        for col,a in enumerate(acc):
            quot=divmonic(a,diag[heights[i]])
            nonzero_v+=bool(quot)
            delta=colours[i]+colours[col]-n-1
            for degree,v in enumerate(quot):
                if v:
                    require(degree>=abs(delta) and (degree-delta)%2==0,"quotient has no integral two-colour lift")
                    max_deg=max(max_deg,degree)
    return {"order":n,"dimension":size,"each_ring_entries":size**2,"nonzero_left":sparse_count,
            "nonzero_quotient":nonzero_v,"max_quotient_degree":max_deg,"det_left":1,"det_quotient":1,
            "rings":["Z[x]","Z[X,Y]"],"input_sha256":hashlib.sha256(src).hexdigest()}

if __name__=="__main__":
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument("files",type=Path,nargs="+")
    for file in arg.parse_args().files:
        print(json.dumps(check(file)),flush=True)
