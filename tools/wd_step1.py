#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wd_step1.py — foundations for replicating Huang–Kosters–Yeo's Table (CRYPTO 2015,
Sec. 5.2): degree of regularity of Weil-descent systems from the bivariate
summation polynomial S_3 over F_{2^n}.

This file: F_{2^n} arithmetic (polynomial basis), ordinary curve
y^2 + xy = x^3 + a2 x^2 + a6, and S_3 in characteristic 2 -- DERIVED from
HKY's general b2/b4/b6/b8 formula (their Sec. 5.1), then VERIFIED on actual
point triples P0 + P1 + P2 = O before anything is built on it.

Derivation (a1=1, a3=0, a4=0 in the ordinary char-2 model):
  b2 = a1^2 + 4a2 = 1,  b4 = a1a3 + 2a4 = 0,  b6 = a3^2 + 4a6 = 0,
  b8 = a1^2 a6 - a1a3a4 + a2a3^2 + 4a2a6 - a4^2 = a6.
  In char 2 the "-2(...)" middle term of S_3 vanishes and signs are irrelevant:
    S_3 = X0^2X1^2 + X0^2X2^2 + X1^2X2^2 + X0X1X2 + a6.
  (a2 drops out entirely -- worth noticing, and worth checking numerically.)
"""

import random

# ---------------------------------------------------------------- F_{2^n}

class GF2n:
    """F_{2^n} = F_2[t]/(modulus); elements are ints (bit i = coeff of t^i)."""
    def __init__(self, n, modulus=None):
        self.n = n
        self.mod = modulus if modulus is not None else self.find_irreducible(n)
        self.mask = (1 << n) - 1

    @staticmethod
    def is_irreducible(poly, n):
        # Ben-Or test: x^(2^i) mod poly, gcd with poly for i <= n/2
        def mulmod(a, b):
            r = 0
            while b:
                if b & 1: r ^= a
                b >>= 1; a <<= 1
                if a >> n & 1: a ^= poly
            return r
        def gcd(a, b):
            while b:
                # polynomial gcd over F_2
                while a.bit_length() >= b.bit_length() and a:
                    a ^= b << (a.bit_length() - b.bit_length())
                a, b = b, a
            return a
        x = 2
        for i in range(1, n // 2 + 1):
            x = mulmod(x, x)
            if gcd(poly, x ^ 2) != 1:
                return False
        return True

    @classmethod
    def find_irreducible(cls, n):
        for poly in range((1 << n) | 1, 1 << (n + 1), 2):
            if cls.is_irreducible(poly, n):
                return poly
        raise ValueError

    def add(self, a, b): return a ^ b
    def mul(self, a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a >> self.n & 1: a ^= self.mod
        return r
    def sq(self, a): return self.mul(a, a)
    def pow(self, a, e):
        r = 1
        while e:
            if e & 1: r = self.mul(r, a)
            a = self.mul(a, a); e >>= 1
        return r
    def inv(self, a):
        assert a != 0
        return self.pow(a, (1 << self.n) - 2)
    def trace(self, a):
        t, x = 0, a
        for _ in range(self.n):
            t ^= x; x = self.sq(x)
        return t  # 0 or 1
    def rand(self, rng): return rng.randrange(1 << self.n)


# ------------------------------------------------- ordinary curve, char 2
# E: y^2 + xy = x^3 + a2 x^2 + a6   (a6 != 0), ordinary (a1 = 1).

INF = None

class BinCurve:
    def __init__(self, F, a2, a6):
        self.F, self.a2, self.a6 = F, a2, a6

    def on_curve(self, P):
        if P is INF: return True
        x, y = P; F = self.F
        lhs = F.add(F.sq(y), F.mul(x, y))
        rhs = F.add(F.add(F.mul(F.sq(x), x), F.mul(self.a2, F.sq(x))), self.a6)
        return lhs == rhs

    def neg(self, P):
        if P is INF: return INF
        x, y = P
        return (x, self.F.add(x, y))

    def add(self, P, Q):
        F = self.F
        if P is INF: return Q
        if Q is INF: return P
        x1, y1 = P; x2, y2 = Q
        if x1 == x2:
            if F.add(y1, y2) == x1:   # Q = -P
                return INF
            # doubling: lambda = x1 + y1/x1
            if x1 == 0: return INF
            lam = F.add(x1, F.mul(y1, F.inv(x1)))
            x3 = F.add(F.add(F.sq(lam), lam), self.a2)
            y3 = F.add(F.add(F.sq(x1), F.mul(lam, x3)), x3)
            return (x3, y3)
        lam = F.mul(F.add(y1, y2), F.inv(F.add(x1, x2)))
        x3 = F.add(F.add(F.add(F.add(F.sq(lam), lam), x1), x2), self.a2)
        y3 = F.add(F.add(F.mul(lam, F.add(x1, x3)), x3), y1)
        return (x3, y3)

    def random_point(self, rng):
        F = self.F
        while True:
            x = F.rand(rng)
            if x == 0: continue
            # y^2 + xy = c  ->  (y/x)^2 + (y/x) = c/x^2 ; solvable iff Tr(c/x^2)=0
            c = F.add(F.add(F.mul(F.sq(x), x), F.mul(self.a2, F.sq(x))), self.a6)
            d = F.mul(c, F.inv(F.sq(x)))
            if F.trace(d) != 0: continue
            # solve z^2 + z = d by half-trace (n odd) or brute force (small n)
            z = self.solve_quad(d)
            if z is None: continue
            return (x, F.mul(z, x))

    def solve_quad(self, d):
        """Solve z^2 + z = d. The map z -> z^2 + z is F_2-LINEAR, so this is a
        linear system over F_2 in the bit-coordinates of z: precompute its n x n
        matrix once, then Gaussian-eliminate. (A brute-force fallback over 2^n
        elements -- the first draft, for even n -- was the cause of a >300 s
        timeout at n=20: ~10^6 field ops per point, hundreds of points.)"""
        F = self.F
        if not hasattr(self, "_Lmat"):
            cols = [F.add(F.sq(1 << i), 1 << i) for i in range(F.n)]  # L(t^i)
            self._Lmat = cols
        cols = self._Lmat
        n = F.n
        # solve sum_i z_i * cols[i] = d over F_2 by elimination on augmented rows
        # build augmented matrix rows: for each bit position r, row = (coeffs of z_i at bit r | d_r)
        rows = []
        for r in range(n):
            coeff = 0
            for i in range(n):
                if cols[i] >> r & 1: coeff |= 1 << i
            rows.append((coeff, d >> r & 1))
        # gaussian elimination
        pivots = {}
        for coeff, rhs in rows:
            for pc, (pcoeff, prhs) in pivots.items():
                if coeff >> pc & 1:
                    coeff ^= pcoeff; rhs ^= prhs
            if coeff == 0:
                if rhs: return None  # inconsistent (trace != 0)
                continue
            pc = (coeff & -coeff).bit_length() - 1
            # eliminate pc from existing pivots
            for k in list(pivots):
                kc, kr = pivots[k]
                if kc >> pc & 1:
                    pivots[k] = (kc ^ coeff, kr ^ rhs)
            pivots[pc] = (coeff, rhs)
        z = 0
        for pc, (coeff, rhs) in pivots.items():
            # after full elimination each pivot row is coeff with only free vars besides pc; set free vars = 0
            if rhs: z |= 1 << pc
        if F.add(F.sq(z), z) == d: return z
        # free variable choices might matter: try flipping the kernel element 1
        z2 = z ^ 1
        return z2 if F.add(F.sq(z2), z2) == d else None


def S3_char2(F, a6, x0, x1, x2):
    """S_3 = X0^2X1^2 + X0^2X2^2 + X1^2X2^2 + X0X1X2 + a6  (char 2, a1=1)."""
    s0, s1, s2 = F.sq(x0), F.sq(x1), F.sq(x2)
    t = F.add(F.add(F.mul(s0, s1), F.mul(s0, s2)), F.mul(s1, s2))
    t = F.add(t, F.mul(F.mul(x0, x1), x2))
    return F.add(t, a6)


if __name__ == "__main__":
    rng = random.Random(1)
    for n in (12, 16, 20):
        F = GF2n(n)
        a2, a6 = F.rand(rng), F.rand(rng) or 1
        E = BinCurve(F, a2, a6)
        # sanity: group law closes on the curve
        P, Q = E.random_point(rng), E.random_point(rng)
        assert E.on_curve(P) and E.on_curve(Q) and E.on_curve(E.add(P, Q))
        # S_3 must vanish exactly when x-coords come from P0+P1+P2 = O
        ok, bad = 0, 0
        for _ in range(300):
            P0, P1 = E.random_point(rng), E.random_point(rng)
            P2 = E.neg(E.add(P0, P1))          # P0 + P1 + P2 = O
            ok += (S3_char2(F, a6, P0[0], P1[0], P2[0]) == 0)
            # negative control: random third x should (almost) never vanish
            bad += (S3_char2(F, a6, P0[0], P1[0], F.rand(rng)) == 0)
        print(f"n={n:>2}  modulus=0x{F.mod:x}  S_3 vanishes on true triples: {ok}/300"
              f"   on random triples: {bad}/300   (a2 does not appear in S_3 -- consistent with derivation)")
