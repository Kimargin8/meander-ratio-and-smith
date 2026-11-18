"""Jensen matching-catalogue oracle for universal two-step-block returns."""
from collections import Counter
from functools import cache
from itertools import product
import argparse
import json
from math import comb
from pathlib import Path


@cache
def matchings(m):
    if m == 0:
        return ((),)
    return tuple((2*j+1,)+tuple(x+1 for x in p)+(0,)+
                 tuple(x+2*j+2 for x in q)
                 for j in range(m) for p in matchings(j) for q in matchings(m-1-j))


def advance(state,move):
    """Top ports run far to near, bottom ports near to far."""
    if state is None:
        return None
    a,p=state
    b=len(p)-a
    if move=="R":
        return (a+1,p) if b else None
    if move=="L":
        return (a-1,p) if a else None
    if move=="U":
        lifted=tuple(v if v<a else v+2 for v in p)
        return a+1,lifted[:a]+(a+1,a)+lifted[a:]
    if move!="D":
        raise ValueError(move)
    if not a or not b or p[a-1]==a:
        return None
    old=list(p)
    x,y=p[a-1],p[a]
    old[x],old[y]=y,x
    return a-1,tuple(v if v<a-1 else v-2 for v in old[:a-1]+old[a+1:])


def walk(state,word):
    for move in word:
        state=advance(state,move)
    return state


def blocks(k,allow_identity):
    for word in product("UDLR",repeat=2*k):
        w="".join(word)
        if allow_identity or all(w[i:i+2] not in ("RL","LR") for i in range(0,2*k,2)):
            yield w


def transfer_row(k):
    e=(1,(1,0))
    plain=Counter({e:1});removed=Counter({e:1});single=Counter({e:1})
    out=[]
    for step in range(k+1):
        out.append((plain[e],removed[e],single[e]))
        if step==k:
            break
        # Each remaining step can erase at most four ports before the final e.
        cap=2+4*(k-step-1)
        for cur,names in [(plain,list(blocks(1,True))),(removed,list(blocks(1,False))),
                          (single,[w for w in blocks(1,True) if w!="LR"])]:
            nxt=Counter()
            for state,count in cur.items():
                for w in names:
                    dest=walk(state,w)
                    if dest is not None and len(dest[1])<=cap:
                        nxt[dest]+=count
            cur.clear()
            cur.update(nxt)
    return out


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write",type=Path,default=Path(__file__).with_name("return_words.json"))
    args=cli.parse_args()
    e=(1,(1,0))
    returns={k:sorted(w for w in blocks(k,False) if walk(e,w)==e) for k in (2,3,4)}
    if [len(ws) for ws in returns.values()]!=[4,10,70]:
        raise ArithmeticError("base return census")
    coverage=[]
    for m in range(1,9):
        states=0
        for p in matchings(m):
            for a in range(1,2*m,2):
                s=(a,p);states+=1
                if walk(s,"LR")!=s or walk(s,"RL")!=s:
                    raise ArithmeticError((m,s,"identity"))
                for ws in returns.values():
                    for w in ws:
                        if walk(s,w)!=s:
                            raise ArithmeticError((m,s,w))
        coverage.append({"port_pairs":m,"all_odd_split_matchings":states})
    values=transfer_row(7)
    ordinary=[1,2,8,42,262,1828,13820,110954]
    if [row[0] for row in values]!=ordinary: raise ArithmeticError("M1..M8")
    for n,(_,d,b) in enumerate(values):
        if d!=sum((-2)**(n-j)*comb(n,j)*ordinary[j] for j in range(n+1)):
            raise ArithmeticError("double identity inverse")
        if b!=sum((-1)**(n-j)*comb(n,j)*ordinary[j] for j in range(n+1)):
            raise ArithmeticError("oriented avoidance inverse")
    result={"coverage":coverage,"states":sum(r["all_odd_split_matchings"] for r in coverage),
            "Q_returns_at_e":returns,"words_returning_every_test_state":{k:len(ws) for k,ws in returns.items()},
            "entries_M_Q_B_at_e":values}
    args.write.write_bytes((json.dumps(result,indent=2)+"\n").encode("ascii"))
    print(json.dumps(result))


if __name__=="__main__":
    main()
