#!/usr/bin/env python3
"""
wd_verify_gb.py -- independent check of wd_vspace on the seed-2025 n=16 instance.

1. Enumerate all solutions of the descended system (brute force).
2. Compute the reduced grevlex Groebner basis of the vanishing ideal of those
   points directly by the Buchberger-Moeller algorithm (pure linear algebra on
   evaluation vectors; no ideal operations at all).
3. Compute V_3 with wd_vspace, put it in reduced row echelon form, and extract
   the rows whose pivot is a minimal element of the pivot set.
4. The two sets must coincide -> V_3 really contains the reduced GB, i.e. the
   instance is solvable without ever exceeding degree 3.
"""
import random
import sys

sys.path.insert(0, ".")
from wd_step2 import build_consistent_instance, DescentSystem
from wd_vspace import MonomialTable, closure, count_solutions
import numpy as np


def all_solutions(nvars, eqs):
    N = nvars
    a = np.arange(1 << N, dtype=np.uint32)
    bits = [((a >> i) & 1).astype(np.uint8) for i in range(N)]
    alive = np.ones(1 << N, dtype=bool)
    for quad, lin, const in eqs:
        v = np.full(1 << N, const, dtype=np.uint8)
        for (i, j) in quad:
            v ^= bits[i] & bits[j]
        for i in lin:
            v ^= bits[i]
        alive &= (v == 0)
    return [int(x) for x in np.nonzero(alive)[0]]


def buchberger_moeller(points, table):
    """Reduced grevlex GB of the ideal of a finite set of points in F_2^N.
    Monomials are processed in increasing order; a monomial is standard iff its
    evaluation vector is independent of those of the standard monomials so far;
    otherwise (if not a multiple of an earlier leading monomial) it gives a GB
    element m + sum(standard monomials)."""
    npts = len(points)
    std_cols = []          # column indices of standard monomials
    # echelon on evaluation vectors (ints with npts bits), tracking combination
    piv = {}               # pivot bit -> (evalvec, combination as poly int)
    gb = []
    lead = []              # leading monomial masks of GB so far
    for k, m in enumerate(table.masks):
        if any((l & m) == l for l in lead):
            continue       # multiple of a leading monomial: not a minimal generator
        ev = 0
        for t, p in enumerate(points):
            if (p & m) == m:
                ev |= 1 << t
        comb = 1 << k
        v = ev
        while v:
            b = v.bit_length() - 1
            if b in piv:
                pv, pc = piv[b]
                v ^= pv
                comb ^= pc
            else:
                break
        if v:
            piv[v.bit_length() - 1] = (v, comb)
            std_cols.append(k)
        else:
            gb.append(comb)   # m + (combination of standard monomials) vanishes on all points
            lead.append(m)
    return gb, std_cols


def reduced_gb_from_space(E, table):
    piv = dict(E.pivots)
    cols = sorted(piv)
    # full reduction (RREF)
    for c in cols:
        r = piv[c]
        v = r ^ (1 << c)
        while v:
            b = v.bit_length() - 1
            if b in piv:
                r ^= piv[b]
                v ^= piv[b]
                v &= (1 << b) - 1  # bits above b already handled
            else:
                v ^= 1 << b
        piv[c] = r
    masks = {c: table.masks[c] for c in cols}
    gb = []
    for c in cols:
        m = masks[c]
        minimal = not any(c2 != c and (masks[c2] & m) == masks[c2] for c2 in cols)
        if minimal:
            gb.append(piv[c])
    return gb


def show(p, table, nvars):
    terms = []
    v = p
    while v:
        b = v.bit_length() - 1
        m = table.masks[b]
        terms.append("*".join(f"X{i}" for i in range(nvars) if m >> i & 1) or "1")
        v ^= 1 << b
    return " + ".join(terms)


rng = random.Random(2025)
for n in (12, 14, 16):
    F, E_, a6, basis, xQ, sol, _ = build_consistent_instance(n, rng)
    if n != 16:
        continue
    S = DescentSystem(F, a6, xQ, basis)
    N = S.nvars
    pts = all_solutions(N, S.eqs)
    print(f"n={n}: {len(pts)} solutions:", [format(p, f'0{N}b')[::-1] for p in pts])
    table = MonomialTable(N, 3)
    gens = [table.from_terms(q, l, c) for (q, l, c) in S.eqs]
    V3 = closure(table, gens, 3)
    gb_space = reduced_gb_from_space(V3, table)
    gb_bm, std = buchberger_moeller(pts, table)
    degs = lambda gb: sorted(table.degree_of_poly(p) for p in gb)
    print(f"  GB from V_3          : {len(gb_space)} elements, degrees {degs(gb_space)}")
    print(f"  GB from points (B-M) : {len(gb_bm)} elements, degrees {degs(gb_bm)}")
    print(f"  identical as sets    : {set(gb_space) == set(gb_bm)}")
    print("  standard monomials   :", [show(1 << k, table, N) for k in std])
    for p in sorted(gb_space, key=lambda p: -table.degree_of_poly(p))[:3]:
        print("  e.g.", show(p, table, N))
