"""Longest positive run conditioned on a uniformly chosen binary site set."""
import argparse
from math import comb
import json
from pathlib import Path

def shorter(N,x,r):
    if r==0:
        return 0
    g=N-x+1
    return sum((-1)**j*comb(g,j)*comb(N-r*j,x-r*j)
               for j in range(1+min(g,x//r)))

def make(rows):
    results={}
    for sn,hist in rows.items():
        N=int(sn)-1
        # Each prescribed set of the same size has the same count.
        each=[]
        for x,row in enumerate(hist):
            q,rem=divmod(sum(row),comb(N,x))
            if rem: raise ArithmeticError("nonintegral occupation weight")
            each.append(q)
        cumul=[sum(v*shorter(N,x,r) for x,v in enumerate(each)) for r in range(N+2)]
        out=[b-a for a,b in zip(cumul,cumul[1:])]
        if min(out)<0 or sum(out)!=sum(map(sum,hist)):
            raise ArithmeticError("invalid distribution")
        if N<=12:
            oracle=[0]*(N+1)
            for mask in range(1<<N):
                longest=max(map(len,format(mask,"b").split("0")))
                oracle[longest]+=each[mask.bit_count()]
            if oracle!=out: raise ArithmeticError("binary word oracle")
        results[sn]=out
    return {"indices":"entry k counts diagrams with maximum positive run exactly k",
            "independent_binary_words_through_number_of_sites":12,"orders":results}

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("joint",type=Path);parser.add_argument("output",type=Path)
    a=parser.parse_args()
    result=make(json.loads(a.joint.read_bytes())["orders"])
    a.output.write_text(json.dumps(result,separators=(",",":"))+"\n",encoding="ascii")
