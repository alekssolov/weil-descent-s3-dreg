#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wd_vspace.py -- intrinsic (monomial-order-free, tool-free) degree measures for
Weil-descent systems over F_2, following HKY (CRYPTO 2015), Section 2:

    V_i = smallest F_2-vector space containing {f in F : deg f <= i} and closed
          under g -> h*g whenever deg(h*g) <= i.                  (Def. 1)

We work in the Boolean ring B = F_2[x]/(x_i^2 + x_i); by the argument in the
notes below V_i(F cup {field eqs}) mod field eqs equals the closure of F under
multiplication by single variables applied to elements of degree <= i-1.

From the V_i we read off:
  * first fall degree  = min c >= 2 with  V_c cap B_{<=c-1} != V_{c-1}
  * last  fall degree  = max c        with  V_c cap B_{<=c-1} != V_{c-1}   (Sec. 2, after Remark 1)
  * D_solve            = min D with reduced Groebner basis of the ideal contained in V_D
                         (equivalently: the leading monomials of V_D generate LM(J);
                          tested by #standard monomials == #solutions).

D_solve is a LOWER BOUND for the max step degree of any F4-type algorithm
(Magma F4 step degree, M4GB lcmdegree): everything such an algorithm produces
at step degree <= D lies in V_D.  When it coincides with M4GB's number the
tool-reported D_reg is certified to be a property of the instance, not of the
implementation.

Representation: a Boolean polynomial is a Python int; bit k <-> the k-th
monomial in grevlex order (degree first, then grevlex).  Highest set bit =
leading monomial.  Monomials themselves are bitmasks over the N variables.
"""

import itertools
import sys
import time

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None


# ----------------------------------------------------------------- monomials

class MonomialTable:
    """Index all squarefree monomials of degree <= dmax in grevlex order."""

    def __init__(self, nvars, dmax):
        self.n, self.dmax = nvars, dmax
        masks = []
        for d in range(dmax + 1):
            layer = []
            for comb in itertools.combinations(range(nvars), d):
                m = 0
                for v in comb:
                    m |= 1 << v
                layer.append(m)
            # grevlex within a degree: m1 > m2  <=>  m1 < m2 as integers
            layer.sort(reverse=True)
            masks.extend(layer)
        self.masks = masks
        self.index = {m: k for k, m in enumerate(masks)}
        self.deg = [bin(m).count("1") for m in masks]
        # first column index of each degree layer (columns are degree-sorted)
        self.deg_start = [0] * (dmax + 2)
        for k, d in enumerate(self.deg):
            self.deg_start[d + 1] = k + 1
        for d in range(1, dmax + 2):
            self.deg_start[d] = max(self.deg_start[d], self.deg_start[d - 1])
        # multiplication tables: mul[a][k] = index of x_a * monomial_k (deg(m_k) < dmax)
        self.mul = []
        for a in range(nvars):
            bit = 1 << a
            row = [None] * len(masks)
            for k, m in enumerate(masks):
                if self.deg[k] < dmax or (m & bit):
                    row[k] = self.index[m | bit]
            self.mul.append(row)

    def degree_of_poly(self, p):
        if p == 0:
            return -1
        return self.deg[p.bit_length() - 1]

    def times_var(self, p, a):
        """x_a * p in the Boolean ring (requires deg p <= dmax-1)."""
        row = self.mul[a]
        out = 0
        while p:
            low = p & -p
            k = low.bit_length() - 1
            out ^= 1 << row[k]
            p ^= low
        return out

    def from_terms(self, quad, lin, const):
        """(quad {(i,j):1}, lin {i:1}, const) as in wd_step2 -> int polynomial."""
        p = 0
        for (i, j) in quad:
            m = (1 << i) | (1 << j)
            p ^= 1 << self.index[m]
        for i in lin:
            p ^= 1 << self.index[1 << i]
        if const:
            p ^= 1 << self.index[0]
        return p


# ------------------------------------------------------------ echelon space

class EchelonSpace:
    """Row-echelon basis of an F_2-subspace of Boolean polynomials.
    pivots: leading-column -> row (each row's highest bit is its pivot)."""

    def __init__(self):
        self.pivots = {}

    def reduce(self, v):
        piv = self.pivots
        while v:
            c = v.bit_length() - 1
            r = piv.get(c)
            if r is None:
                return v
            v ^= r
        return 0

    def insert(self, v):
        """Reduce v; if nonzero add as new pivot row. Returns the new row or 0."""
        v = self.reduce(v)
        if v:
            self.pivots[v.bit_length() - 1] = v
        return v

    def dim(self):
        return len(self.pivots)

    def dim_upto_degree(self, table, j):
        """dim(V cap B_{<=j}) -- valid because columns are degree-sorted."""
        lim = table.deg_start[j + 1]
        return sum(1 for c in self.pivots if c < lim)


def closure(table, gens, cap, seed_rows=None):
    """Compute V_cap as an EchelonSpace.

    gens: list of int polynomials (all of degree <= cap are used).
    seed_rows: optional rows of V_{cap-1} (already in V_cap) -- only their
    degree-(cap-1) members need to be multiplied further.
    """
    E = EchelonSpace()
    todo = []
    if seed_rows:
        for r in seed_rows:
            r2 = E.insert(r)
            if r2 and table.degree_of_poly(r2) == cap - 1:
                todo.append(r2)
    for g in gens:
        if table.degree_of_poly(g) <= cap:
            r = E.insert(g)
            if r and table.degree_of_poly(r) <= cap - 1:
                todo.append(r)
    n = table.n
    while todo:
        r = todo.pop()
        for a in range(n):
            v = E.insert(table.times_var(r, a))
            if v and table.degree_of_poly(v) <= cap - 1:
                todo.append(v)
    return E


# ---------------------------------------------------------- solution count

def count_solutions(nvars, eqs):
    """Brute-force number of F_2^N solutions of a list of (quad, lin, const)."""
    N = nvars
    total = 1 << N
    if np is None or N > 22:
        # slow generic fallback
        cnt = 0
        for a in range(total):
            bits = [(a >> i) & 1 for i in range(N)]
            ok = True
            for quad, lin, const in eqs:
                v = const
                for (i, j) in quad:
                    v ^= bits[i] & bits[j]
                for i in lin:
                    v ^= bits[i]
                if v:
                    ok = False
                    break
            cnt += ok
        return cnt
    a = np.arange(total, dtype=np.uint32)
    bits = [((a >> i) & 1).astype(np.uint8) for i in range(N)]
    alive = np.ones(total, dtype=bool)
    for quad, lin, const in eqs:
        v = np.full(total, const, dtype=np.uint8)
        for (i, j) in quad:
            v ^= bits[i] & bits[j]
        for i in lin:
            v ^= bits[i]
        alive &= (v == 0)
    return int(alive.sum())


def count_standard(nvars, pivot_masks):
    """#squarefree monomials NOT divisible by any mask in pivot_masks."""
    N = nvars
    if np is None:
        raise RuntimeError("numpy needed")
    nonstd = np.zeros(1 << N, dtype=bool)
    for m in pivot_masks:
        nonstd[m] = True
    # superset closure: m nonstandard if any subset obtained by removing a variable is
    for a in range(N):
        v = nonstd.reshape(1 << (N - 1 - a), 2, 1 << a)
        v[:, 1, :] |= v[:, 0, :]
    return int((1 << N) - nonstd.sum())


# --------------------------------------------------------------- analysis

def analyze(nvars, eqs, dmax=6, verbose=False, nsol=None):
    """Return dict with first/last fall degree, D_solve, dims, #solutions."""
    t0 = time.time()
    if nsol is None:
        nsol = count_solutions(nvars, eqs)
    table = MonomialTable(nvars, dmax)
    gens = [table.from_terms(q, l, c) for (q, l, c) in eqs]
    gens = [g for g in gens if g]
    dims = {}          # dims[c] = dim V_c
    dims_low = {}      # dims_low[c] = dim(V_c cap B_{<=c-1})
    std = {}
    prev_rows = None
    D_solve = None
    dims[1] = 0
    for c in range(2, dmax + 1):
        E = closure(table, gens, c, seed_rows=prev_rows)
        dims[c] = E.dim()
        dims_low[c] = E.dim_upto_degree(table, c - 1)
        piv = [table.masks[k] for k in E.pivots]
        std[c] = count_standard(nvars, piv)
        if verbose:
            print(f"   V_{c}: dim={dims[c]:6d}  dim(V_{c} cap B<={c-1})={dims_low[c]:6d}"
                  f"  #std(LM V_{c})={std[c]:6d}  (#sol={nsol})  [{time.time()-t0:.1f}s]",
                  flush=True)
        prev_rows = list(E.pivots.values())
        if std[c] == nsol:
            D_solve = c
            break
    falls = [c for c in dims_low if dims_low[c] > dims[c - 1]]
    ffd = min(falls) if falls else None
    lfd = max(falls) if falls else 0
    return {
        "nsol": nsol, "D_solve": D_solve, "first_fall": ffd, "last_fall": lfd,
        "dims": dims, "dims_low": dims_low, "std": std, "time": time.time() - t0,
    }


if __name__ == "__main__":
    # smoke test on a tiny system: x0 x1 + x2 = 0, x1 + x2 + 1 = 0, x0 x2 + x0 + x1 = 0
    eqs = [({(0, 1): 1}, {2: 1}, 0), ({}, {1: 1, 2: 1}, 1), ({(0, 2): 1}, {0: 1, 1: 1}, 0)]
    print(analyze(3, eqs, dmax=3, verbose=True))
