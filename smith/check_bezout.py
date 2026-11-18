"""Generate and verify the univariate weights in the integral row construction.

Usage: python check_bezout.py [N] [--write file.json] [--verify file.json]
N is k+1; lists encode polynomials by ascending powers, with [] denoting zero.
Python standard library only.
"""

import argparse
import json
from fractions import Fraction
from math import prod
from pathlib import Path


def trim(a):
    a = list(a)
    while a and not a[-1]:
        a.pop()
    return tuple(a)


def add(a, b):
    size = max(len(a), len(b))
    return trim((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                for i in range(size))


def neg(a):
    return tuple(-v for v in a)


def mul(a, b):
    if not a or not b:
        return ()
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return trim(result)


def product(polynomials):
    result = (1,)
    for f in polynomials:
        result = mul(result, f)
    return result


def divide(a, b):
    if not b:
        raise ZeroDivisionError
    a = list(a)
    q = [0] * max(0, 1 + len(a) - len(b))
    while a and len(a) >= len(b):
        shift = len(a) - len(b)
        c = a[-1] if b[-1] == 1 else Fraction(a[-1], b[-1])
        q[shift] = c
        for j, y in enumerate(b):
            a[j + shift] -= c*y
        while a and not a[-1]:
            a.pop()
    return trim(q), trim(a)


def integral(a):
    if any(Fraction(x).denominator != 1 for x in a):
        raise ArithmeticError("nonintegral coefficient in a Bezout identity")
    return tuple(map(int, a))


def bezout(a, b):
    old_r, r = a, b
    old_s, s = (1,), ()
    old_t, t = (), (1,)
    while r:
        q, new_r = divide(old_r, r)
        old_r, r = r, new_r
        old_s, s = s, add(old_s, neg(mul(q, s)))
        old_t, t = t, add(old_t, neg(mul(q, t)))
    if len(old_r) != 1:
        raise ArithmeticError("the two quotients share a root")
    s = integral([Fraction(c, old_r[0]) for c in old_s])
    t = integral([Fraction(c, old_r[0]) for c in old_t])
    if add(mul(s, a), mul(t, b)) != (1,):
        raise ArithmeticError("Bezout identity failed")
    return s, t


def beraha_factors(N):
    """Extract primitive factors of the Chebyshev recurrence, with B_2(z)=z."""
    prev, current = (), (1,)
    primitive, result = {}, {2:(0, 1)}
    for n in range(2, N + 1):
        prev, current = current, add((0,) + current, neg(prev))
        known = product(primitive[d] for d in range(2, n) if n % d == 0)
        f, remainder = divide(current, known)
        if remainder:
            raise ArithmeticError("primitive polynomial division failed")
        primitive[n] = f
        if n > 2:
            if any(f[1::2]):
                raise ArithmeticError("a primitive factor has incorrect parity")
            result[n] = f[::2]
    return result


def certificate(N, factors):
    # If U,V are downward closed under divisibility, all pairs in
    # (U\V) x (V\U) are incomparable. Their factor products are comaximal.
    covered = set()
    coeffs = {}
    for r in range(max(2, N//2+1), N+1):
        divisors = {d for d in range(2, r+1) if r % d == 0}
        if not coeffs:
            covered = divisors
            coeffs[r] = (1,)
            continue
        left = product(factors[d] for d in sorted(divisors - covered))
        right = product(factors[d] for d in sorted(covered - divisors))
        u, v = bezout(left, right)
        coeffs = {i:mul(u, a) for i, a in coeffs.items()}
        coeffs[r] = v
        covered |= divisors
    if covered != set(range(2, N+1)):
        raise ArithmeticError("the divisibility cover is incomplete")
    return {"N":N, "weights":{str(i):list(a) for i,a in coeffs.items()}}


def evaluate(f, z):
    value = 0
    for c in reversed(f):
        value = value*z + c
    return value


def verify(rows):
    if not rows:
        raise ValueError("empty certificate")
    factors = beraha_factors(max(row["N"] for row in rows))
    for row in rows:
        N = row["N"]
        if N < 2 or set(row["weights"]) != {str(r) for r in range(max(2, N//2+1), N+1)}:
            raise ValueError("wrong set of indices for N="+str(N))
        total = ()
        for rs, a in row["weights"].items():
            if any(type(v) is not int for v in a):
                raise TypeError("nonintegral certificate")
            t = product(factors[d] for d in range(2, N+1) if int(rs) % d)
            total = add(total, mul(a, t))
        if total != (1,):
            raise ValueError("polynomial verification failed for N="+str(N))
        # An independent definition of the primitive factors at square inputs:
        # invert the divisor product directly on the integer recurrence.
        for x in (3, 4, 7):
            p, q = 0, 1
            for d in range(2, N+1):
                p, q = q, x*q-p
                value = prod(x if j == 2 else evaluate(factors[j], x*x)
                             for j in range(2, d+1) if d % j == 0)
                if value != q:
                    raise ArithmeticError("integer recurrence check failed")
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int, nargs="?", default=16)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path)
    group.add_argument("--verify", type=Path)
    args = parser.parse_args()
    if args.verify:
        data = json.loads(args.verify.read_text(encoding="utf-8"))
    else:
        if args.N < 2:
            parser.error("N must be at least 2")
        factors = beraha_factors(args.N)
        data = {"schema":1, "variable":"z", "coefficient_order":"ascending",
                "certificates":[certificate(N,factors) for N in range(2, args.N+1)]}
    if data["schema"] != 1 or data["variable"] != "z" or data["coefficient_order"] != "ascending":
        raise ValueError("unsupported certificate format")
    count = verify(data["certificates"])
    if args.write:
        args.write.write_text(json.dumps(data,separators=(",",":"))+"\n", encoding="utf-8")
    print(f"Verified {count} polynomial identities. Largest N = "
          f"{max(r['N'] for r in data['certificates'])}.")


if __name__ == "__main__":
    main()
