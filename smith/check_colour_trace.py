"""Rebuild shaded Jones--Wenzl projectors from turnback equations over Z[X,Y].

Use Python 3.10+ and SymPy 1.14.0. No project modules or input tables.
Ports run clockwise from the upper left; the first m lie on the top.
"""
import argparse
from functools import cache
from hashlib import sha256
import json
from pathlib import Path

from sympy.polys.domains import ZZ
from sympy.polys.matrices import DomainMatrix


R = ZZ.poly_ring("X", "Y")
X, Y = R.gens
WEIGHTS = (X, Y)


def require(condition, name):
    if not condition:
        raise ArithmeticError(name)


@cache
def matchings(m):
    if m == 0:
        return ((),)
    result = []
    for left_pairs in range(m):
        j = 2 * left_pairs + 1
        for u in matchings(left_pairs):
            for v in matchings(m - 1 - left_pairs):
                result.append((j,) + tuple(k+1 for k in u) + (0,) +
                              tuple(k+j+1 for k in v))
    return tuple(result)


def reconnect(p, i, j):
    """Apply an adjacent cap/cup generator; return an unweighted diagram."""
    if p[i] == j:
        return p
    out = list(p)
    a, b = p[i], p[j]
    out[a], out[b] = b, a
    out[i], out[j] = j, i
    return tuple(out)


def turnback_system(m, shade, top):
    states = matchings(m)
    equations = {}
    for i in range(m-1):
        a, b = (i, i+1) if top else (2*m-1-i, 2*m-2-i)
        # Both upper and lower discs enclose the gap after i+1 strands
        # from the left. Its colour alternates with i.
        loop = WEIGHTS[shade ^ ((i+1) % 2)]
        for j, pattern in enumerate(states):
            out = reconnect(pattern, a, b)
            eq = equations.setdefault((i,out), {})
            eq[j] = loop if pattern[a] == b else R.one
    return equations


def projector(m, shade):
    """Solve e_i f = 0 and f[identity] = 1, retaining no trace input."""
    states = matchings(m)
    identity = tuple(reversed(range(2*m)))
    pivot = states.index(identity)
    top = turnback_system(m,shade,True)
    if m <= 1:
        return (R.one,), R.one, 0
    chosen = []
    for state in states:
        if state == identity:
            continue
        i = next(k for k in range(m-1) if state[k] == k+1)
        chosen.append(top[i,state])
    chosen.append({pivot:R.one})
    dense = {(i,j):v for i,row in enumerate(chosen) for j,v in row.items()}
    size = len(states)
    matrix = DomainMatrix.from_dok(dense,(size,size),R)
    rhs = DomainMatrix.from_dok({(size-1,0):R.one},(size,1),R)
    solution, denominator = matrix.solve_den(rhs)
    nums = [row[0] for row in solution.to_list()]
    common = denominator
    for val in nums:
        common = common.gcd(val)
    denominator = denominator.exquo(common)
    nums = [val.exquo(common) for val in nums]
    if denominator.LC < 0:
        denominator, nums = -denominator, [-v for v in nums]
    require(nums[pivot] == denominator,"identity coefficient")
    bottom = turnback_system(m,shade,False)
    for eq in (*top.values(), *bottom.values()):
        require(not sum((v*nums[j] for j,v in eq.items()),R.zero),"turnback residual")
    return tuple(nums), denominator, len(top)+len(bottom)


def close_one(p, side):
    """Delete the joined boundary ports and bridge partners, recording a free loop."""
    m = len(p)//2
    i,j = (m-1,m) if side=="right" else (0,2*m-1)
    closed = p[i]==j
    q = reconnect(p,i,j)
    keep = [k for k in range(2*m) if k not in (i,j)]
    return tuple(keep.index(q[k]) for k in keep), closed


def quantum(upto, shade):
    """Positive-weight form of (-1)^(j-1)[j] for the indicated left shade."""
    q = [R.zero,R.one]
    for j in range(1,upto):
        q.append(WEIGHTS[shade^(j%2)]*q[-1] - q[-2])
    return q


@cache
def expected_den(m,shade):
    q = quantum(m,shade)
    primitive = {}
    target = R.one
    for j in range(2,m+1):
        val=q[j]
        for d in range(2,j):
            if not j%d:
                val=val.exquo(primitive[d])
        primitive[j]=val
        if (m+1)%j:
            target *= val
    return target


def sphere_face_counts(p, q, shade):
    """Glue two matching disks and count checkerboard faces by boundary intervals."""
    N=len(p)
    parent=list(range(N))
    def find(i):
        while i!=parent[i]:
            parent[i]=parent[parent[i]]
            i=parent[i]
        return i
    for disk in (p,q):
        for i,j in enumerate(disk):
            if i<j:
                parent[find((i-1)%N)]=find(j)
                parent[find(i)]=find((j-1)%N)
    colours=[set(),set()]
    for i in range(N):
        colours[shade^((i+1)%2)].add(find(i))
    require(colours[0].isdisjoint(colours[1]),"sphere colouring")
    return len(colours[0]),len(colours[1])


def check(upto):
    data={}
    report=[]
    mutations={}
    for m in range(upto+1):
        states=matchings(m)
        for shade in range(2):
            nums,den,ann=projector(m,shade)
            data[m,shade]=nums,den
            require(den==expected_den(m,shade),"common reduced denominator")
            if m==0:
                continue
            q=quantum(m+1,shade)
            checks=0
            for side in ("right","left"):
                remaining_shade=shade if side=="right" else shade^1
                fixed_shade=shade if side=="right" else shade^(m%2)
                disc_shade=shade^(m%2) if side=="right" else shade
                partial={p:R.zero for p in matchings(m-1)}
                for p,v in zip(states,nums):
                    reduced,disk=close_one(p,side)
                    partial[reduced]+=v*(WEIGHTS[disc_shade] if disk else R.one)
                target, lower_den = data[m-1,remaining_shade]
                qs=quantum(m+1,fixed_shade)
                for k,p in enumerate(matchings(m-1)):
                    require(partial[p]*lower_den*qs[m] == den*target[k]*qs[m+1],
                            side+" partial trace")
                    checks+=1
            identity=tuple(reversed(range(2*m)))
            predicted=WEIGHTS[shade]*q[m+1]
            observed=None
            for exterior in states:
                total=R.zero
                for p,num in zip(states,nums):
                    a,b=sphere_face_counts(p,exterior,shade)
                    total+=num*X**a*Y**b
                require(total==(den*predicted if exterior==identity else R.zero),
                        "all exterior matching traces")
                if exterior==identity:
                    observed=total.exquo(den)
            if m==1 and shade==0:
                require(observed!=q[2],"outer-face control")
                mutations["omit_spherical_face"]={"order":1,"left":"X","correct":"X*Y","mutated":"Y"}
            if m==2 and shade==0:
                wrong=quantum(m,1)[2]
                require(den!=wrong,"wrong-shade denominator control")
                require(observed!=WEIGHTS[shade^1]*q[m+1],
                        "swapped exterior control")
                mutations["swap_denominator_shade"]={"order":2,"left":"X","correct":str(den),"mutated":str(wrong)}
                mutations["swap_outermost_face"]={"order":2,"left":"X","correct":str(predicted),
                                                  "mutated":str(WEIGHTS[shade^1]*q[m+1])}
            report.append({"strands":m,"left_shade":str(WEIGHTS[shade]),"basis":len(states),
                           "left_and_right_turnback_equations":ann,
                           "left_and_right_partial_trace_coefficients":checks,
                           "spherical_matchings_checked":len(states),"least_denominator":str(den),
                           "unique_nonzero_spherical_closure":str(predicted)})
    return {"ring":"Z[X,Y], with projector coefficients in Q(X,Y)",
            "max_strands":upto,"results":report,"incorrect_normalizations_detected":mutations}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-strands",type=int,choices=range(1,7),default=6)
    parser.add_argument("--write",type=Path)
    args=parser.parse_args()
    results=check(args.max_strands)
    text=json.dumps(results,indent=2)+"\n"
    if args.write:
        args.write.write_bytes(text.encode("utf-8"))
    print(text if not args.write else "PASS "+sha256(text.encode()).hexdigest())
