#!/usr/bin/env python3
"""
wd_m4gb_emul.py -- Buchberger over F_2 with EXPLICIT field equations, whose
pair management mirrors m4gb.hpp::update()/selection() (M4GB, 2018-01-31):

  * selection: all critical pairs of minimal lcm degree (cap MAXSELECTION=512)
  * update(p):
      B-criterion on old pairs (lines 1601-1612)
      new pairs (p,g): 'bad' if LMs coprime (product criterion, line 1627)
      M-criterion (lines 1634-1644):  pair A is dropped if ANY other new pair B
        has lcm(B) | lcm(A)  -- NON-STRICT divisibility, so two good pairs with
        EQUAL lcm eliminate each other (Gebauer-Moeller would keep exactly one)
      IMMEDIATE_BASIS_REDUCE (lines 1649-1686): old g with LM(p) | LM(g) is
        replaced by g reduced by p; pair (p,g) erased.
  * --fixed : Gebauer-Moeller F-criterion done right (keep one per equal lcm,
              none if a coprime pair shares that lcm)

Monomials live in the polynomial ring with exponents <= 2 only via the field
equations x_a^2+x_a, which are basis elements: lcm(g, x_a^2+x_a) = x_a*LM(g)
(contains x_a^2).  We represent an lcm as (squarefree mask, mask of squared
variables); polynomials themselves are kept Boolean-reduced (that is exactly
reduction by the field equations).
"""
import random
import sys
import time

sys.path.insert(0, ".")
from wd_step2 import build_consistent_instance, DescentSystem
from wd_vspace import MonomialTable, count_solutions, count_standard

MAXSELECTION = 512


def popcount(x):
    return bin(x).count("1")


class Lcm:
    __slots__ = ("m", "s")

    def __init__(self, m, s=0):
        self.m, self.s = m, s

    def degree(self):
        return popcount(self.m) + popcount(self.s)

    def divides(self, o):
        return (self.m & ~o.m) == 0 and (self.s & ~o.s) == 0

    def __eq__(self, o):
        return self.m == o.m and self.s == o.s

    def __hash__(self):
        return hash((self.m, self.s))


def key_of(x):
    """basis key: ('p', mask) for a polynomial, ('f', a) for x_a^2+x_a."""
    return x


def lm_lcm(k1, k2):
    """lcm of leading monomials of two basis keys (polynomial ring)."""
    if k1[0] == "p" and k2[0] == "p":
        return Lcm(k1[1] | k2[1], 0)
    if k1[0] == "f" and k2[0] == "f":
        return Lcm((1 << k1[1]) | (1 << k2[1]), (1 << k1[1]) | (1 << k2[1]))
    p, f = (k1, k2) if k1[0] == "p" else (k2, k1)
    return Lcm(p[1] | (1 << f[1]), 1 << f[1])


def deg_of(k):
    return popcount(k[1]) if k[0] == "p" else 2


class M4GBEmul:
    def __init__(self, table, gens, fixed=False):
        self.T = table
        self.N = table.n
        self.fixed = fixed
        self.basis = {}          # ('p', mask) -> poly ;  ('f', a) -> None
        self.shadow = {}         # removed elements (IMMEDIATE_BASIS_REDUCE) still referenced by pairs
        self.CP = {}             # (k1, k2) -> Lcm
        self.log = []
        for a in range(self.N):
            self.basis[("f", a)] = None
        for g in gens:
            self.update(g)

    # ---------------------------------------------------------- polynomials
    def lm_mask(self, p):
        return self.T.masks[p.bit_length() - 1]

    def mul_mono(self, p, t):
        if t == 0:
            return p
        out = 0
        T = self.T
        while p:
            low = p & -p
            k = low.bit_length() - 1
            out ^= 1 << T.index[T.masks[k] | t]
            p ^= low
        return out

    def divisor(self, m):
        s = m
        while True:
            if ("p", s) in self.basis:
                return s
            if s == 0:
                return None
            s = (s - 1) & m

    def reduce_full(self, p):
        result = 0
        while p:
            k = p.bit_length() - 1
            m = self.T.masks[k]
            d = self.divisor(m)
            if d is None:
                result ^= 1 << k
                p ^= 1 << k
            else:
                p ^= self.mul_mono(self.basis[("p", d)], m ^ d)
        return result

    # ------------------------------------------------------------- update()
    def update(self, p):
        p = self.reduce_full(p)
        if p == 0:
            return False
        m = self.lm_mask(p)
        pk = ("p", m)
        pdeg = popcount(m)
        pL = Lcm(m, 0)
        self.basis[pk] = p
        # B-criterion on old pairs
        for pair in list(self.CP):
            L = self.CP[pair]
            if pdeg < L.degree() and pL.divides(L) and lm_lcm(pk, pair[1]) != L and lm_lcm(pk, pair[0]) != L:
                del self.CP[pair]
        # new pairs
        new = []
        for gk in self.basis:
            if gk == pk:
                continue
            L = lm_lcm(pk, gk)
            coprime = (pdeg + deg_of(gk) == L.degree())
            new.append([gk, L, not coprime, coprime])      # [key, lcm, good, coprime]
        # order like std::map<crit_pair_t>: by lcm (degree, then monomial), then p2
        new.sort(key=lambda e: (e[1].degree(), e[1].m, e[1].s, str(e[0])))
        if not self.fixed:
            # --- M4GB as written: any other new pair whose lcm divides mine (non-strict) kills me
            for i, e in enumerate(new):
                if not e[2]:
                    continue
                for j, e2 in enumerate(new):
                    if e2[1].degree() > e[1].degree():
                        break
                    if i != j and e2[1].divides(e[1]):
                        e[2] = False
                        break
        else:
            # --- Gebauer-Moeller: proper divisor kills; equal lcm: keep the first good one,
            #     unless a coprime pair has that lcm (then drop all)
            for i, e in enumerate(new):
                if not e[2]:
                    continue
                for j, e2 in enumerate(new):
                    if e2[1].degree() > e[1].degree():
                        break
                    if i == j or not e2[1].divides(e[1]):
                        continue
                    if e2[1] != e[1]:
                        e[2] = False; break            # proper divisor
                    if e2[3]:
                        e[2] = False; break            # equal lcm with a coprime pair
                    if j < i and e2[2]:
                        e[2] = False; break            # an earlier kept pair with equal lcm
        for gk, L, good, _ in new:
            if good:
                self.CP[(pk, gk)] = L
        # IMMEDIATE_BASIS_REDUCE
        redo = []
        for gk in list(self.basis):
            if gk[0] != "p" or gk == pk:
                continue
            if pdeg < popcount(gk[1]) and (m & gk[1]) == m:
                self.CP.pop((pk, gk), None)      # M4GB erases only the pair (p, g); other pairs of g survive
                g = self.basis.pop(gk)
                self.shadow[gk] = g
                redo.append(g)
        for g in redo:
            self.update(g)                       # M4GB: goes to to_add, then update()
        return True

    # ------------------------------------------------------------ selection
    def poly_of(self, k):
        return self.basis[k] if k in self.basis else self.shadow[k]

    def spoly(self, k1, k2, L):
        if k1[0] == "f" and k2[0] == "f":
            return 0
        if k1[0] == "f" or k2[0] == "f":
            pk, fk = (k2, k1) if k1[0] == "f" else (k1, k2)
            return self.mul_mono(self.poly_of(pk), 1 << fk[1])      # x_a * g  (Boolean form of S(g, x_a^2+x_a))
        m1, m2 = k1[1], k2[1]
        lcm = m1 | m2
        return self.mul_mono(self.poly_of(k1), lcm ^ m1) ^ self.mul_mono(self.poly_of(k2), lcm ^ m2)

    def run(self, verbose=True):
        while self.CP:
            d = min(L.degree() for L in self.CP.values())
            batch = [pair for pair, L in self.CP.items() if L.degree() == d]
            batch.sort(key=lambda pr: (self.CP[pr].m, self.CP[pr].s, str(pr)))
            batch = batch[:MAXSELECTION]
            ncp_before = len(self.CP)
            rows = []
            for pair in batch:
                L = self.CP.pop(pair)
                k1, k2 = pair
                s = self.reduce_full(self.spoly(k1, k2, L))
                if s:
                    rows.append(s)
            # joint row reduction of the batch (M4GB: matrix.rowreduce), then update() one by one
            ech = {}
            reduced = []
            for r in rows:
                v = r
                while v:
                    c = v.bit_length() - 1
                    if c in ech:
                        v ^= ech[c]
                    else:
                        ech[c] = v
                        reduced.append(v)
                        break
            new = 0
            nb = len([k for k in self.basis if k[0] == "p"])
            for r in sorted(reduced, key=lambda x: x.bit_length()):
                if self.update(r):
                    new += 1
            nlin = sum(1 for k in self.basis if k[0] == "p" and popcount(k[1]) == 1)
            self.log.append((d, len(batch), len(reduced), new))
            if verbose:
                print(f"   lcmdegree={d}  selected={len(batch):4d} (of {ncp_before})  nonzero rows={len(reduced):4d}"
                      f"  basis {nb}->{len([k for k in self.basis if k[0]=='p'])}  linear={nlin}", flush=True)
        return self.log


if __name__ == "__main__":
    fixed = "--fixed" in sys.argv
    which = [int(a) for a in sys.argv[1:] if a.isdigit()] or [14, 16]
    rng = random.Random(2025)
    for n in (12, 14, 16, 18, 20):
        F, E, a6, basis, xQ, sol, _ = build_consistent_instance(n, rng)
        if n not in which:
            continue
        S = DescentSystem(F, a6, xQ, basis)
        N = S.nvars
        T = MonomialTable(N, 5)
        gens = [T.from_terms(q, l, c) for (q, l, c) in S.eqs]
        nsol = count_solutions(N, S.eqs)
        print(f"n={n}: #sol={nsol}   [{'Gebauer-Moeller fixed' if fixed else 'M4GB criteria as written'}]")
        t0 = time.time()
        B = M4GBEmul(T, gens, fixed=fixed)
        B.run()
        std = count_standard(N, [k[1] for k in B.basis if k[0] == "p"])
        maxd = max(d for d, sel, nz, new in B.log if new > 0)
        print(f"   done {time.time()-t0:.1f}s  GB complete: {std == nsol}   max lcmdegree with new polys: {maxd}")
