"""Exact small-order checks for the two weights on alternating meander triples."""
import argparse
from collections import Counter, defaultdict
from functools import cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_return_words import advance, matchings


E = (1, (1, 0))
MOVES = "UDLR"
FLIP = str.maketrans("LR", "RL")
SLABS = tuple(product(MOVES, repeat=2))


def exponents(previous, center, following):
    # Positive = adjacent upper arch on the left, adjacent lower on the right.
    plus = previous in "UR" and center == "L" and following in "DR"
    minus = previous in "UL" and center == "R" and following in "DL"
    return int(plus), int(minus)


def weighted_walk(state, previous, word):
    plus = minus = 0
    for i in range(0, len(word), 2):
        center, following = word[i:i+2]
        p, q = exponents(previous, center, following)
        plus += p
        minus += q
        state = advance(advance(state, center), following)
        previous = following
    return state, previous, plus, minus


def catalogue_tests(max_pairs):
    words = {"R": ("ULDR", "ULLRDR"), "D": ("LUURDD", "UULLDRRD"),
             "U_lower": ("URDDLU", "ULDRRDLU")}
    words["L"] = tuple(w.translate(FLIP) for w in words["R"])
    words["U_upper"] = tuple(w.translate(FLIP) for w in words["U_lower"])
    tests = Counter()
    for m in range(1, max_pairs+1):
        for pairing in matchings(m):
            for upper in range(1, 2*m, 2):
                state = upper, pairing
                for case, candidates in words.items():
                    if case.startswith("U_"):
                        if pairing[upper-1] != upper:
                            continue
                        side = upper if case == "U_upper" else len(pairing)-upper
                        if side < 3:
                            continue
                    previous = case[0]
                    for w in candidates:
                        if weighted_walk(state, previous, w) != (state, previous, 0, 0):
                            raise ValueError(("invalid unmarked return", state, case, w))
                        tests[case] += 1
    if any(weighted_walk(E, "U", w)[0] == E for case in ("U_upper", "U_lower")
           for w in words[case]):
        raise ValueError("start exception")
    return {"max_port_pairs": max_pairs, "words": words, "word_state_tests": dict(tests)}


def score_arch_pairs(n):
    """Compare signs on individual Dyck pairs, including disconnected systems."""
    catalogue = []
    for p in matchings(n):
        directions = tuple(p[i] > i for i in range(2*n))
        left = sum(1 << i for i in range(n-1) if p[2*i+1] == 2*i)
        right = sum(1 << i for i in range(n-1) if p[2*i+1] == 2*i+2)
        catalogue.append((directions, left, right))
    labels = {(True, True): "U", (False, False): "D",
              (True, False): "R", (False, True): "L"}
    for top, top_l, top_r in catalogue:
        for bottom, bottom_l, bottom_r in catalogue:
            word = [labels[pair] for pair in zip(top,bottom)]
            for i in range(n-1):
                arch = (int(bool((top_l & bottom_r) & (1 << i))),
                        int(bool((bottom_l & top_r) & (1 << i))))
                if exponents(*word[2*i:2*i+3]) != arch:
                    raise ValueError(("individual signed site", n, word, i, arch))
    return (n-1)*len(catalogue)**2


@cache
def options(state, previous):
    result = []
    for center, following in SLABS:
        dest = advance(advance(state, center), following)
        if dest is not None:
            result.append((dest, following, *exponents(previous, center, following)))
    return result


def table(n):
    frontier = {(E,"U"): Counter({(0,0): 1})}
    for blocks_used in range(1,n):
        future = 2*(n-1-blocks_used)
        nxt = defaultdict(Counter)
        for (state, previous), poly in frontier.items():
            for dest, following, p, q in options(state,previous):
                upper, pairing = dest
                if max(upper, len(pairing)-upper) > 1+future:
                    continue
                target = nxt[dest,following]
                for (a,b), value in poly.items():
                    target[a+p,b+q] += value
        frontier = nxt
    total = Counter()
    for (state,_), poly in frontier.items():
        if state == E:
            total.update(poly)
    return [[total[i,j] for j in range(n-i)] for i in range(n)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--joint",type=Path,
                        default=Path(__file__).with_name("joint_counts.json"))
    parser.add_argument("--write",type=Path,default=Path(__file__).with_name("weighted_checks.json"))
    args = parser.parse_args()
    source = args.joint.read_bytes()
    existing = json.loads(source)["orders"]
    certificate = {"joint_sha256": sha256(source).hexdigest(),
                   "convention": "R moves the nearest lower end to the upper side. The positive activity counts L-centered sites (upper adjacent arch then lower adjacent arch). Reflect the side names for the paper's R from left to right with left called upper.",
                   "returns": catalogue_tests(8),
                   "individual_signed_sites_by_order": {n:score_arch_pairs(n) for n in range(1,8)},
                   "orientation_control": {}, "transfer_orders": {}}
    for word,expected in [("ULRD",(1,0)),("URLD",(0,1))]:
        found = weighted_walk(E,word[0],word[1:-1])
        if found[0] != E or found[2:] != expected:
            raise ValueError("individual sign control")
        certificate["orientation_control"][word] = found[2:]
    for n in range(1,12):
        computed = table(n)
        if computed != existing[str(n)]:
            raise ValueError(("joint discrepancy",n))
        certificate["transfer_orders"][n] = computed
        options.cache_clear()
    args.write.write_bytes((json.dumps(certificate,indent=2)+"\n").encode("ascii"))
    print("Every coefficient agrees through order 11. Return tests:",
          certificate["returns"]["word_state_tests"])


if __name__ == "__main__":
    main()
