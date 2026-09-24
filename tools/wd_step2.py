#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wd_step2.py — Weil descent of S_3(X0, X1, x(Q)) over F_{2^n} down to a quadratic
system over F_2, following Huang-Kosters-Yeo (CRYPTO 2015, Sec. 5.2).

Setup (theirs): E/F_{2^n} ordinary, Q a random nonzero point, X0 and X1 forced
into a random F_2-subspace V of dimension ceil(n/2). Their description of the
resulting system: "about n variables and about n quadratic equations together
with field equations Y_i^2 + Y_i".

Why the descended system is QUADRATIC (worth spelling out -- it is the whole
reason this approach is tractable at all in char 2): squaring is F_2-LINEAR, so
writing X0 = sum u_i v_i with u_i in F_2 and v_i a basis of V gives
    X0^2 = sum u_i^2 v_i^2 = sum u_i v_i^2       (u_i^2 = u_i for u_i in F_2)
i.e. X0^2 is LINEAR in the bit variables. With
    S_3 = X0^2 X1^2 + X0^2 xQ^2 + X1^2 xQ^2 + X0 X1 xQ + a6   (char 2, a1 = 1)
the terms X0^2 X1^2 and X0 X1 xQ are (linear)x(linear) = quadratic, and
X0^2 xQ^2, X1^2 xQ^2 are linear. So every descended equation has degree exactly 2.

The descent itself: S_3 evaluated at the bit variables is an element of F_{2^n}
of the form  sum_{i,j} u_i w_j c_ij + sum_i u_i d_i + sum_j w_j e_j + a6.
Taking each of the n polynomial-basis coordinates gives n equations over F_2.

Verification (run below): build a CONSISTENT instance -- pick P0, P1 with
x-coords actually in V, set xQ = x(P0+P1) -- and check that the bit vector of
(x(P0), x(P1)) satisfies all n descended equations. A system that cannot
recognise its own known solution is worthless for measuring anything.
"""

import random
import sys

sys.path.insert(0, ".")
from wd_step1 import GF2n, BinCurve, S3_char2, INF


def random_subspace(F, dim, rng):
    """Random F_2-subspace of F_{2^n} of given dimension: returns a basis."""
    while True:
        basis = [F.rand(rng) for _ in range(dim)]
        if independent(F, basis):
            return basis


def independent(F, vecs):
    """F_2-linear independence of elements of F_{2^n} (as bit vectors)."""
    pivots = []
    for v in vecs:
        for p in pivots:
            v = min(v, v ^ p)
        if v == 0:
            return False
        pivots.append(v)
        pivots.sort(reverse=True)
    return True


def subspace_elements(F, basis):
    """All 2^dim elements of the subspace (only used for small dim, to find
    curve points whose x lies in V)."""
    els = [0]
    for b in basis:
        els += [e ^ b for e in els]
    return els


def coords_in_subspace(F, basis, x):
    """Express x in the basis over F_2; return tuple of bits, or None."""
    # Gaussian elimination tracking combinations
    rows = [(b, 1 << i) for i, b in enumerate(basis)]
    cur = x
    combo = 0
    pivots = []
    for v, tag in rows:
        pivots.append([v, tag])
    pivots.sort(key=lambda r: -r[0].bit_length())
    for r in pivots:
        if cur == 0:
            break
        if cur.bit_length() == r[0].bit_length():
            cur ^= r[0]
            combo ^= r[1]
    if cur != 0:
        return None
    return tuple((combo >> i) & 1 for i in range(len(basis)))


class DescentSystem:
    """The n quadratic equations over F_2 obtained by Weil descent.

    Each equation is stored as (quad, lin, const) where quad is a dict
    {(i,j): 1} over variable index pairs, lin a dict {i: 1}, const in {0,1}.
    Variables: 0..nv-1 are the u's (for X0), nv..2*nv-1 the w's (for X1).
    """

    def __init__(self, F, a6, xQ, basis):
        self.F, self.n, self.nv = F, F.n, len(basis)
        self.nvars = 2 * self.nv
        c = {}
        for i, vi in enumerate(basis):
            for j, vj in enumerate(basis):
                # coefficient of u_i * w_j :  v_i^2 v_j^2 + v_i v_j xQ
                t1 = F.mul(F.sq(vi), F.sq(vj))
                t2 = F.mul(F.mul(vi, vj), xQ)
                c[(i, j)] = F.add(t1, t2)
        xQ2 = F.sq(xQ)
        d = [F.mul(F.sq(v), xQ2) for v in basis]   # coefficient of u_i
        e = [F.mul(F.sq(v), xQ2) for v in basis]   # coefficient of w_j
        self.eqs = []
        for r in range(self.n):
            quad = {}
            for (i, j), val in c.items():
                if (val >> r) & 1:
                    quad[(i, self.nv + j)] = 1
            lin = {}
            for i, val in enumerate(d):
                if (val >> r) & 1:
                    lin[i] = lin.get(i, 0) ^ 1
            for j, val in enumerate(e):
                if (val >> r) & 1:
                    k = self.nv + j
                    lin[k] = lin.get(k, 0) ^ 1
            lin = {k: v for k, v in lin.items() if v}
            const = (a6 >> r) & 1
            self.eqs.append((quad, lin, const))

    def evaluate(self, assignment):
        """assignment: tuple of bits, length nvars. Returns list of n bits."""
        out = []
        for quad, lin, const in self.eqs:
            v = const
            for (i, j) in quad:
                v ^= assignment[i] & assignment[j]
            for i in lin:
                v ^= assignment[i]
            out.append(v)
        return out

    def stats(self):
        nq = sum(len(q) for q, _, _ in self.eqs)
        nl = sum(len(l) for _, l, _ in self.eqs)
        return nq, nl


def build_consistent_instance(n, rng, seed_shift=0):
    """Return (F, E, a6, basis, xQ, known_solution_bits) with x(P0), x(P1) in V."""
    F = GF2n(n)
    a2, a6 = F.rand(rng), (F.rand(rng) or 1)
    E = BinCurve(F, a2, a6)
    dim = (n + 1) // 2
    for _ in range(200):
        basis = random_subspace(F, dim, rng)
        els = subspace_elements(F, basis)
        pts = []
        for x in els:
            if x == 0:
                continue
            c = F.add(F.add(F.mul(F.sq(x), x), F.mul(a2, F.sq(x))), a6)
            dd = F.mul(c, F.inv(F.sq(x)))
            if F.trace(dd) != 0:
                continue
            z = E.solve_quad(dd)
            if z is None:
                continue
            pts.append((x, F.mul(z, x)))
            if len(pts) >= 2:
                break
        if len(pts) >= 2:
            P0, P1 = pts[0], pts[1]
            R = E.add(P0, P1)
            if R is INF:
                continue
            xQ = R[0]
            u = coords_in_subspace(F, basis, P0[0])
            w = coords_in_subspace(F, basis, P1[0])
            if u is None or w is None:
                continue
            return F, E, a6, basis, xQ, tuple(u) + tuple(w), (P0, P1, R)
    raise RuntimeError("could not build a consistent instance")


if __name__ == "__main__":
    rng = random.Random(7)
    print("Weil descent of S_3(X0,X1,xQ) with X0,X1 in a random F_2-subspace V, dim V = ceil(n/2)")
    print(f"  {'n':>3} {'dim V':>6} {'#vars':>6} {'#eqs':>5} {'quad terms':>11} {'lin terms':>10}"
          f"  {'known soln satisfies all eqs':>30}  {'random assign zero?':>20}")
    for n in (8, 10, 12, 14, 16):
        F, E, a6, basis, xQ, sol, (P0, P1, R) = build_consistent_instance(n, rng)
        # independent check that the instance really is a relation
        assert S3_char2(F, a6, P0[0], P1[0], xQ) == 0, "instance is not a true relation"
        S = DescentSystem(F, a6, xQ, basis)
        res = S.evaluate(sol)
        ok = all(b == 0 for b in res)
        # negative control: a random assignment should not satisfy all n equations
        bad = 0
        for _ in range(200):
            a = tuple(rng.randrange(2) for _ in range(S.nvars))
            if all(b == 0 for b in S.evaluate(a)):
                bad += 1
        nq, nl = S.stats()
        print(f"  {n:>3} {len(basis):>6} {S.nvars:>6} {len(S.eqs):>5} {nq:>11} {nl:>10}"
              f"  {str(ok):>30}  {bad:>3}/200")
    print("\n  All equations are degree exactly 2 (see module docstring for why squaring being")
    print("  F_2-linear forces this). #vars ~ n and #eqs = n, matching HKY's description.")
