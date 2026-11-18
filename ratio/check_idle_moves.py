"""Check odd-cut return moves and count return-word avoiders in matching pairs."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def matchings(size):
    """Generate noncrossing partners via balanced parentheses."""
    partners = [-1] * (2 * size)

    def visit(position, opened, stack):
        if position == 2 * size:
            yield tuple(partners)
            return
        if opened < size:
            yield from visit(position + 1, opened + 1, stack + (position,))
        if stack:
            other = stack[-1]
            partners[position] = other
            partners[other] = position
            yield from visit(position + 1, opened, stack[:-1])

    yield from visit(0, 0, ())


def pair_oracle(size):
    masks = []
    for p in matchings(size):
        up_down = 0
        down_up = 0
        for i, a in enumerate(range(1, 2 * size - 1, 2)):
            if p[a] > a and p[a+1] < a+1:
                up_down |= 1 << i
            if p[a] < a and p[a+1] > a+1:
                down_up |= 1 << i
        masks.append((p, up_down, down_up))
    connected = 0
    avoiders = 0
    for p, a, b in masks:
        for q, c, d in masks:
            vertex = q[p[0]]
            visited = 1
            while vertex:
                vertex = q[p[vertex]]
                visited += 1
            if visited == size:
                connected += 1
                if not ((a & d) | (b & c)):
                    avoiders += 1
    return connected, avoiders


def canonical(left, right):
    relabel = {}
    ends = tuple(relabel.setdefault(x, len(relabel)) for x in left + right)
    return ends[:len(left)], ends[len(left):]


def step(state, letter):
    left, right = state
    if letter == "U":
        fresh = len(left) + len(right)
        return canonical((fresh,) + left, (fresh,) + right)
    if letter == "R":
        return canonical(left[1:], left[:1] + right) if left else None
    if letter == "L":
        return canonical(right[:1] + left, right[1:]) if right else None
    if not left or not right or left[0] == right[0]:
        return None
    source, target = left[0], right[0]
    merged = tuple(target if x == source else x for x in left[1:] + right[1:])
    return canonical(merged[:len(left)-1], merged[len(left)-1:])


def walk(state, word):
    for letter in word:
        state = step(state, letter)
        if state is None:
            break
    return state


def transfer_check(levels, exact):
    words = ["".join(p) for p in product("UDLR", repeat=2)]
    idle = ("LR", "RL", "URDL", "UURDDL")
    base = ((0,), (0,))
    full = Counter({base: 1})
    restricted = Counter({base: 1})
    verified = set()
    widths = []
    for m in range(levels + 1):
        require(full[base] == exact[m][0], f"full transfer order {m+1}")
        require(restricted[base] == exact[m][1], f"restricted transfer order {m+1}")
        widths.append(len(full))
        new_full, new_restricted = Counter(), Counter()
        for state, count in full.items():
            if state not in verified:
                require(all(len(side) % 2 for side in state), "odd side parity")
                for pattern in idle:
                    require(walk(state, pattern) == state, f"invalid return {pattern}")
                verified.add(state)
            if m == levels:
                continue
            for word in words:
                result = walk(state, word)
                if result is not None:
                    new_full[result] += count
                    if word not in idle[:2] and state in restricted:
                        new_restricted[result] += restricted[state]
        full, restricted = new_full, new_restricted
    # UD tries to close its fresh arc. A pair of identical shifts can run out of ends.
    require(walk(base, "UD") is None and walk(base, "RR") is None, "control words")
    return {"words": list(idle), "distinct_odd_states_checked": len(verified),
            "reachable_widths": widths, "control_words_rejected_at_base": ["UD", "RR"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("counts", type=Path, help="two columns n and M_n")
    ap.add_argument("--output", type=Path)
    ap.add_argument("--through", type=int, default=8)
    ap.add_argument("--transfer-through", type=int, default=6)
    args = ap.parse_args()
    require(1 <= args.transfer_through <= args.through, "invalid audit ranges")
    raw = args.counts.read_bytes()
    known = dict(tuple(map(int, line.split())) for line in raw.decode("utf-8-sig").splitlines()
                 if line.strip() and not line.startswith("#"))
    exact = []
    for n in range(1, args.through + 1):
        counts = pair_oracle(n)
        c = sum(comb(n-1, j)*(-2)**j*known[n-j] for j in range(n))
        require(counts == (known[n], c), f"matching enumeration order {n}")
        exact.append(counts)
    data = {"input_sha256": sha256(raw).hexdigest(), "first_order": 1,
            "matching_orders": args.through, "M_C_by_order": exact,
            "transfer_orders": args.transfer_through,
            "transfer": transfer_check(args.transfer_through - 1, exact)}
    output = json.dumps(data, indent=2) + "\n"
    if args.output:
        args.output.write_text(output, encoding="ascii")
    print(output, end="")


if __name__ == "__main__":
    main()
