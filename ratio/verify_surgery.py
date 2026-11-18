"""Independent finite checks for adjacent-arch surgery. Python 3.10+, standard library."""

import json
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from math import comb, prod
from pathlib import Path
from platform import python_version


def matchings(n):
    """Stack-match all balanced bracket words with n opening symbols."""
    for opens in combinations(range(2*n), n):
        markers = set(opens)
        stack = []
        mate = [-1]*(2*n)
        for v in range(2*n):
            if v in markers:
                stack.append(v)
            elif stack:
                u = stack.pop()
                mate[u], mate[v] = v, u
            else:
                break
        else:
            yield tuple(mate)


def one_component(top, bottom):
    """Contract each upper edge, then use graph reachability on the quotient."""
    visited = {0, top[0]}
    pending = [0, top[0]]
    while pending:
        v = bottom[pending.pop()]
        if v not in visited:
            u = top[v]
            visited.update((u, v))
            pending.extend((u, v))
    return len(visited) == len(top)


def compatible(top, bottom):
    plus = minus = 0
    for center in range(1, len(top)-1, 2):
        flag = 1 << (center//2)
        if (top[center], bottom[center]) == (center-1, center+1):
            plus |= flag
        elif (top[center], bottom[center]) == (center+1, center-1):
            minus |= flag
    return plus, minus


def intervals(n):
    return [(a, b) for a in range(2*n) for b in range(a+2, 2*n, 2)]


def exposed(matching, blocks):
    return sum(1 << k for k,(a,b) in enumerate(blocks)
               if sum(not (a <= matching[u] <= b) for u in range(a,b+1)) == 1)


def disjoint_families(blocks, pos=0, mask=0, orders=()):
    yield mask, orders
    for i in range(pos,len(blocks)):
        a,b = blocks[i]
        after = i+1
        while after<len(blocks) and blocks[after][0]<=b:
            after += 1
        yield from disjoint_families(blocks, after, mask | (1 << i), orders+((b-a)//2,))


@lru_cache(None)
def vertex_groups(n, chosen):
    group = [0]
    for v in range(1,2*n):
        anchor = (v-1)//2
        collapsed = bool(chosen & (1 << anchor))
        group.append(group[-1]+int(not collapsed))
    return group


def collapse(matching, groups):
    reduced = [-1]*(1+groups[-1])
    ports = [-1]*len(reduced)
    for u,v in enumerate(matching):
        i,j = groups[u],groups[v]
        if i!=j:
            assert ports[i] == -1
            reduced[i],ports[i] = j,u
    assert -1 not in ports
    return tuple(reduced),ports


def expand(small, groups, signs_positive, colour):
    n = len(groups)
    large = [-1]*n
    # The first/second edge at each site alternates with the prescribed colour.
    for a in range(0,n-2,2):
        if groups[a] == groups[a+2]:
            first_colour = 0 if signs_positive & (1 << (a//2)) else 1
            u = a + int(colour != first_colour)
            large[u],large[u+1] = u+1,u
    ports = [-1]*len(small)
    for u,partner in enumerate(large):
        if partner == -1:
            i = groups[u]
            assert ports[i] == -1
            ports[i] = u
    assert -1 not in ports
    for i,j in enumerate(small):
        large[ports[i]] = ports[j]
    return tuple(large)


def check_surgery(top, bottom, positive, negative, catalogs):
    active = positive | negative
    chosen = active
    total = 0
    while chosen:
        group = vertex_groups(len(top)//2,chosen)
        t, _ = collapse(top,group)
        b, _ = collapse(bottom,group)
        assert t in catalogs[len(t)//2] and b in catalogs[len(t)//2] and one_component(t,b)
        assert (expand(t,group,positive & chosen,0),expand(b,group,positive & chosen,1)) == (top,bottom)
        if chosen == positive:
            assert compatible(t,b)[0] == 0
        total += 1
        chosen = (chosen-1) & active
    return total


def join(left, right):
    t,b = left
    v,w = right
    d = len(t)
    top = list(t)+[e+d for e in v]
    i,j = top[d-1],top[d]
    top[d-1],top[d] = d,d-1
    top[i],top[j] = j,i
    return tuple(top),b+tuple(e+d for e in w)


def split_at(joined, d):
    top,bottom = joined
    bridges = [(u,v) for u,v in enumerate(top) if u<d<=v]
    assert (d-1,d) in bridges and len(bridges)==2
    i,j = next(pair for pair in bridges if pair != (d-1,d))
    top = list(top)
    top[i],top[d-1],top[d],top[j] = d-1,i,j,d
    assert all(v<d for v in bottom[:d])
    return ((tuple(top[:d]),bottom[:d]),
            (tuple(x-d for x in top[d:]),tuple(x-d for x in bottom[d:])))


def run(max_order=8, surgery_through=7):
    catalogs = {n:set(matchings(n)) for n in range(1,max_order+1)}
    counts, avoiding, report = {}, {}, []
    stored_diagrams = {}
    for n in range(1,max_order+1):
        patterns, block_patterns, positive_sets = Counter(),Counter(),Counter()
        blocks = intervals(n)
        exposed_mask = {t:exposed(t,blocks) for t in catalogs[n]}
        stored_diagrams[n] = []
        surgeries = 0
        for t,b in product(catalogs[n], repeat=2):
            if not one_component(t,b):
                continue
            plus,minus = compatible(t,b)
            patterns[plus,minus] += 1
            positive_sets[plus] += 1
            block_patterns[exposed_mask[t] & exposed_mask[b]] += 1
            if n<max_order:
                stored_diagrams[n].append(((t,b),plus.bit_count(),minus.bit_count()))
            if n<=surgery_through:
                surgeries += check_surgery(t,b,plus,minus,catalogs)
        total = counts[n] = sum(patterns.values())
        avoiding[n] = positive_sets[0]
        for mask in range(1 << (n-1)):
            assert positive_sets[mask] == avoiding[n-mask.bit_count()]
        signed_tests = 0
        for signs in product((0,1,-1),repeat=n-1):
            pos = sum(1<<i for i,s in enumerate(signs) if s==1)
            neg = sum(1<<i for i,s in enumerate(signs) if s==-1)
            adjacent = (pos << 1) & neg or (neg << 1) & pos
            actual = sum(v for (p,m),v in patterns.items() if p&pos==pos and m&neg==neg)
            assert actual == (0 if adjacent else counts[n-(pos|neg).bit_count()])
            signed_tests += 1
        families = 0
        for mask,orders in disjoint_families(blocks):
            if mask:
                actual = sum(v for p,v in block_patterns.items() if p&mask==mask)
                assert actual == counts[n-sum(orders)]*prod(counts[j+1] for j in orders)
                families += 1
        for chosen in range(1 << (n-1)):
            exterior = ((chosen << 1 | chosen >> 1) & ((1 << (n-1))-1)) & ~chosen
            d, k = exterior.bit_count(), chosen.bit_count()
            actual = sum(v for (p,m),v in patterns.items() if not p and m&chosen==chosen)
            expected = sum(comb(d,a)*avoiding[n-k-a] for a in range(d+1))
            assert actual == expected
        glued = all_seams = both_safe = 0
        sign_changes = Counter()
        for r in range(1,n):
            for (a,ap,am),(b,bp,bm) in product(stored_diagrams[r],stored_diagrams[n-r]):
                t,h = join(a,b)
                assert t in catalogs[n] and h in catalogs[n] and one_component(t,h)
                assert split_at((t,h),2*r)==(a,b)
                pos,neg = compatible(t,h)
                dp,dm = pos.bit_count()-ap-bp,neg.bit_count()-am-bm
                assert dp in {-1,0} and dm in {0,1}
                sign_changes[dp,dm] += 1
                all_seams += 1
                if ap==0 and bp==0:
                    assert not pos
                    glued += 1
                    if not am and not bm and r>=2:
                        if a[1][-1] != 2*r-2:
                            result = (t,h)
                        else:
                            assert a[0][-1] != 2*r-2
                            h,t = join(a[::-1],b[::-1])
                            result = (t,h)
                        assert compatible(*result)==(0,0) and one_component(*result)
                        top_cross = any(v>=2*r for v in result[0][:2*r])
                        if top_cross:
                            assert split_at(result,2*r)==(a,b)
                        else:
                            x,y=split_at(result[::-1],2*r)
                            assert (x[::-1],y[::-1])==(a,b)
                        both_safe += 1
        report.append(dict(n=n, ordered_pairs=len(catalogs[n])**2, M=total, B=avoiding[n],
                           signed_prescriptions=signed_tests, positive_subsets=len(positive_sets),
                           nonempty_disjoint_block_families=families,
                           contraction_expansion_roundtrips=surgeries,
                           avoiding_seam_roundtrips=glued, arbitrary_seam_roundtrips=all_seams,
                           two_sign_avoidance=patterns[0,0], avoiding_both_seam_roundtrips=both_safe,
                           conditioned_subset_identities=1 << (n-1),
                           seam_sign_changes={str(key):val for key,val in sorted(sign_changes.items())}))
    return report


if __name__=="__main__":
    print(json.dumps(dict(parameters=dict(max_order=8, surgery_through=7, rooted=True,
                                         swap_colours=False, wraparound=False),
                          python=python_version(),
                          script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                          input="Generated from bracket words; no sequence data is read.",
                          orders=run()),indent=2))
