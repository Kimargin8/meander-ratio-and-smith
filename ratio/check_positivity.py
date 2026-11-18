"""Verify the sharp shift using exact finite differences and two-sign avoidance."""
from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path


def compose(outer,inner,n):
    # Horner's rule modulo z^(n+1); both series have constant term zero.
    v=[0]*(n+1)
    for coefficient in reversed(outer):
        v=[sum(v[i-j]*inner[j] for j in range(1,i+1)) for i in range(n+1)]
        v[0]+=coefficient
    return v


def run(mfile,histfile):
    mbytes=mfile.read_bytes()
    m={int(a):int(b) for a,b in (r.split() for r in mbytes.decode("utf-8-sig").splitlines()
         if r.strip() and not r.startswith("#"))}
    hist=json.loads(histfile.read_bytes())["orders"]
    n=max(m); seq=[m[i] for i in range(1,n+1)]
    # First member of each successive row of (E-2) finite differences.
    differences=[];c=[0]
    while seq:
        differences.append(seq);c.append(seq[0])
        if min(seq)<0: raise ValueError("negative shifted difference")
        seq=[v-2*u for u,v in zip(seq,seq[1:])]

    # The positive tree recursion for h = z + q.
    q=[0]*(n+1);squared=[0]*(n+1)
    for i in range(1,n+1):
        squared[i]=sum(q[j]*q[i-j] for j in range(1,i))
        q[i]=2*(i==3)+q[i-1]+(4*q[i-2] if i>=2 else 0)+squared[i]+2*squared[i-1]
    h=q.copy();h[1]+=1
    avoid=[0]+[hist[str(i)][0][0] for i in range(1,n+1)]
    if compose(avoid,h,n)!=c: raise ValueError("decoration disagrees with the signed difference")

    def shifted(i,t):
        return sum(comb(i-1,j)*(-t)**j*m[i-j] for j in range(i))
    for i in range(1,n+1):
        if c[i]!=shifted(i,2): raise ValueError("difference indexing")
        B=shifted(i,1); expected=sum(hist[str(i)][0])
        if B!=expected or B!=sum(comb(i-1,j)*c[j+1] for j in range(i)):
            raise ValueError("double smoothing")
    for k,row in enumerate(differences):
        for i,value in enumerate(row,1):
            if value!=sum(comb(i-1,j)*2**(i-1-j)*c[j+k+1] for j in range(i)):
                raise ValueError("shifted positivity representation")
    if shifted(2,Fraction(201,100))!=Fraction(-1,100):
        raise ValueError("the sharpness control failed")
    swapped=avoid.copy();swapped[3]+=1
    if compose(swapped,h,n)==c: raise ValueError("mutated double-avoidance count went undetected")
    return {"M_source_sha256":sha256(mbytes).hexdigest(),"orders":n,
            "checked_differences":sum(map(len,differences)),
            "direct_joint_counts_through":8,"first_order":1,"C":c[1:],
            "negative_control_t_201_over_100_at_order_2":"-1/100",
            "corrupted_double_avoider_at_order_3_rejected":True}


if __name__=="__main__":
    ap=ArgumentParser(description=__doc__);ap.add_argument("input",type=Path)
    ap.add_argument("joint",type=Path);ap.add_argument("--output",type=Path)
    a=ap.parse_args();res=run(a.input,a.joint)
    if a.output: a.output.write_text(json.dumps(res,indent=2)+"\n",encoding="utf-8")
    print("orders",res["orders"],"differences",res["checked_differences"],"C38",res["C"][-1])
