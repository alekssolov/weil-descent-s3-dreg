#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wd_instances.py -- instance generation for the HKY Sec. 5.2 replication, with
two fixes relative to wd_step2.build_consistent_instance:

  * P0, P1 are sampled UNIFORMLY from the factor base F_V = {P : x(P) in V}
    (build_consistent_instance takes the first two points of a fixed
    enumeration of V, which biases the solution to be very sparse);
  * a mode='randomQ' construction that follows HKY literally: Q is a random
    nonzero point of E(k); the system then has 0, 2, 4, ... solutions
    (a random Q of a curve over F_{2^n} decomposes as P0+P1 with P0,P1 in F_V
    only with probability ~ 1 - e^{-1/2} for even n).

count_solutions_bilinear exploits the bilinear structure (only u_i*w_j
products occur): fix u, solve a linear system in w.  Works for any n.
"""

import random
import sys

sys.path.insert(0, ".")
from wd_step1 import GF2n, BinCurve, INF
from wd_step2 import random_subspace, DescentSystem


def coords_in_subspace(F, basis, x):
    """Bits (u_0..u_{d-1}) with x = sum u_i basis[i], or None if x not in V.
    (wd_step2.coords_in_subspace assumes an echelonised basis and can return
    None for genuine members of V; this version does proper elimination.)"""
    piv = {}                      # leading bit -> (vector, tag)
    for i, b in enumerate(basis):
        v, tag = b, 1 << i
        while v:
            lb = v.bit_length() - 1
            if lb in piv:
                pv, pt = piv[lb]
                v ^= pv; tag ^= pt
            else:
                piv[lb] = (v, tag)
                break
    v, tag = x, 0
    while v:
        lb = v.bit_length() - 1
        if lb not in piv:
            return None
        pv, pt = piv[lb]
        v ^= pv; tag ^= pt
    return tuple((tag >> i) & 1 for i in range(len(basis)))


def random_factor_base_point(F, E, basis, rng, exclude_x=()):
    """Uniform random point of E with x-coordinate in V \\ {0}."""
    dim = len(basis)
    while True:
        u = rng.randrange(1, 1 << dim)
        x = 0
        for i in range(dim):
            if u >> i & 1:
                x ^= basis[i]
        if x in exclude_x:
            continue
        c = F.add(F.add(F.mul(F.sq(x), x), F.mul(E.a2, F.sq(x))), E.a6)
        d = F.mul(c, F.inv(F.sq(x)))
        if F.trace(d) != 0:
            continue
        z = E.solve_quad(d)
        if z is None:
            continue
        if rng.randrange(2):        # pick one of the two points +-P uniformly
            z ^= 1
        return (x, F.mul(z, x))


def make_instance(n, rng, mode="consistent", curve=None, basis=None):
    """Return dict with F, E, basis, xQ, eqs, nvars, solution (or None), meta.

    mode='consistent': Q = P0 + P1 with P0, P1 uniform in F_V (distinct x).
    mode='randomQ'   : Q uniform random nonzero point of E (HKY's setup).
    curve=(a2,a6) / basis: reuse fixed curve / subspace (for dependence studies).
    """
    F = GF2n(n)
    if curve is None:
        a2, a6 = F.rand(rng), (F.rand(rng) or 1)
    else:
        a2, a6 = curve
    E = BinCurve(F, a2, a6)
    dim = (n + 1) // 2
    if basis is None:
        basis = random_subspace(F, dim, rng)
    sol = None
    if mode == "consistent":
        while True:
            P0 = random_factor_base_point(F, E, basis, rng)
            P1 = random_factor_base_point(F, E, basis, rng, exclude_x=(P0[0],))
            Q = E.add(P0, P1)
            if Q is not INF:
                break
        u = coords_in_subspace(F, basis, P0[0])
        w = coords_in_subspace(F, basis, P1[0])
        sol = tuple(u) + tuple(w)
    elif mode == "randomQ":
        Q = E.random_point(rng)
    else:
        raise ValueError(mode)
    xQ = Q[0]
    S = DescentSystem(F, a6, xQ, basis)
    if sol is not None:
        assert all(v == 0 for v in S.evaluate(sol)), "known solution violates system"
    return {"F": F, "E": E, "basis": basis, "xQ": xQ, "eqs": S.eqs, "nvars": S.nvars,
            "nv": len(basis), "sol": sol, "a2": a2, "a6": a6, "mode": mode}


def count_solutions_bilinear(nvars, nv, eqs):
    """#solutions in F_2^nvars, using that quadratic terms are u_i*w_j only."""
    n_eq = len(eqs)
    # per equation: c[i] = bitmask over j of coefficients of u_i w_j; d = mask over i; e = mask over j
    C, D, EE, K = [], [], [], []
    for quad, lin, const in eqs:
        c = [0] * nv
        d = 0
        e = 0
        for (i, j) in quad:
            assert i < nv <= j, "equation is not bilinear in (u, w)"
            c[i] |= 1 << (j - nv)
        for i in lin:
            if i < nv:
                d |= 1 << i
            else:
                e |= 1 << (i - nv)
        C.append(c); D.append(d); EE.append(e); K.append(const)
    total = 0
    for u in range(1 << nv):
        # rows: (coeff mask over w, rhs)
        piv = {}
        consistent = True
        for r in range(n_eq):
            row = EE[r]
            cr = C[r]
            uu = u
            i = 0
            while uu:
                if uu & 1:
                    row ^= cr[i]
                uu >>= 1
                i += 1
            rhs = K[r] ^ (bin(D[r] & u).count("1") & 1)
            while row:
                b = row.bit_length() - 1
                if b in piv:
                    prow, prhs = piv[b]
                    row ^= prow
                    rhs ^= prhs
                else:
                    piv[b] = (row, rhs)
                    break
            if row == 0 and rhs:
                consistent = False
                break
        if consistent:
            total += 1 << (nv - len(piv))
    return total


if __name__ == "__main__":
    from wd_vspace import count_solutions
    rng = random.Random(5)
    for n in (8, 10, 12, 14):
        for mode in ("consistent", "randomQ"):
            I = make_instance(n, rng, mode)
            a = count_solutions_bilinear(I["nvars"], I["nv"], I["eqs"])
            b = count_solutions(I["nvars"], I["eqs"])
            print(f"n={n:>2} {mode:>10}: #sol bilinear={a:>3} brute={b:>3}  {'OK' if a == b else 'MISMATCH'}")
