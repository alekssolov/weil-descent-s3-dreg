#!/usr/bin/env python3
"""
wd_buchberger.py -- plain Buchberger in the Boolean ring, normal selection
strategy (all pairs of minimal lcm degree per iteration), fully interreduced
basis, NO pair-elimination criteria (so nothing is skipped).

Pairs: (b,b') with lcm(LM b, LM b') of degree d -> S = (lcm/LM b) b + (lcm/LM b') b'
       (b, x_a^2+x_a) for x_a | LM(b)             -> S = x_a * b  (lcm degree deg b + 1)
Question: after ALL pairs of lcm degree <= 3 are exhausted, does the basis
contain the linear polynomials that V_3 contains?  Prints per-iteration
(lcmdegree, #new) like M4GB.
"""
import random
import sys
import time

sys.path.insert(0, ".")
from wd_step2 import build_consistent_instance, DescentSystem
from wd_vspace import MonomialTable, closure, count_solutions, count_standard


def popcount(x):
    return bin(x).count("1")


class Buchberger:
    def __init__(self, table, gens, maxdeg, strategy="lcm", product_criterion=False, interreduce=True):
        self.product_criterion = product_criterion
        self.interreduce = interreduce
        self.T = table
        self.N = table.n
        self.maxdeg = maxdeg
        self.strategy = strategy   # "lcm": key = lcm degree; "sugar": key = sugar degree
        self.basis = {}         # LM mask -> poly (fully reduced)
        self.sugar = {}         # LM mask -> sugar of that basis element
        self.origin = {}        # LM mask -> origin tag
        self._tag = "input"
        self.lin_parents = []
        self.log = []
        self.pairs = {}         # key -> list of (kind, m1, m2, lcmdeg)
        for g in gens:
            self.insert(g, popcount(self.T.masks[g.bit_length()-1]) if g else 0)

    # ---- polynomial helpers (int bitset representation, columns = table)
    def lm_mask(self, p):
        return self.T.masks[p.bit_length() - 1]

    def mul_mono(self, p, t):
        """p * monomial t (bitmask) in the Boolean ring."""
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
        """a basis LM dividing monomial mask m (enumerate submasks), or None."""
        s = m
        while True:
            if s in self.basis:
                return s
            if s == 0:
                return None
            s = (s - 1) & m

    def reduce_full(self, p, sug=0):
        """Full (top + tail) reduction modulo the basis; returns (remainder, sugar)."""
        result = 0
        while p:
            k = p.bit_length() - 1
            m = self.T.masks[k]
            d = self.divisor(m)
            if d is None:
                result ^= 1 << k
                p ^= 1 << k
            else:
                p ^= self.mul_mono(self.basis[d], m ^ d)   # m/d = m without d's bits
                if self.strategy != "label":
                    sug = max(sug, self.sugar[d] + popcount(m ^ d))
        return result, sug

    def add_pairs(self, m):
        deg = popcount(m)
        sg = self.sugar[m]
        for a in range(self.N):
            if m >> a & 1:
                key = deg + 1 if self.strategy == "lcm" else sg + 1
                self.pairs.setdefault(key, []).append(("field", m, a, deg + 1))
        for m2 in list(self.basis):
            if m2 == m:
                continue
            lcm = m | m2
            if self.product_criterion and (m & m2) == 0:
                continue
            ld = popcount(lcm)
            psug = max(sg + popcount(lcm ^ m), self.sugar[m2] + popcount(lcm ^ m2))
            key = ld if self.strategy == "lcm" else psug
            self.pairs.setdefault(key, []).append(("pair", m, m2, ld))

    def insert(self, p, sug):
        p, sug = self.reduce_full(p, sug)
        if p == 0:
            return False
        m = self.lm_mask(p)
        # elements whose LM is divisible by the new LM leave the basis and are re-inserted (re-reduced)
        redo = []
        if self.interreduce:
            redo = [(self.basis.pop(old), self.sugar.pop(old)) for old in list(self.basis) if (old & m) == m and old != m]
        self.basis[m] = p
        self.sugar[m] = sug
        self.origin[m] = self._tag
        if popcount(m) == 1 and self._tag != "input":
            self.lin_parents.append(self._tag)
        self.add_pairs(m)
        # tail-reduce the remaining elements by the new one (their LMs stay: not divisible by m)
        for k in (list(self.basis) if self.interreduce else []):
            if k == m:
                continue
            q = self.basis.pop(k)
            q2, s2 = self.reduce_full(q, self.sugar[k])
            assert q2 and self.lm_mask(q2) == k
            self.basis[k] = q2
            self.sugar[k] = s2
        for q, sq in redo:
            old_tag = self._tag
            self._tag = "redo(" + old_tag + ")" if not old_tag.startswith("redo") else old_tag
            self.insert(q, sq)
            self._tag = old_tag
        return True

    def run(self, verbose=True):
        while self.pairs:
            d = min(self.pairs)
            if d > self.maxdeg:
                break
            batch = self.pairs.pop(d)
            new = 0
            maxlcm = 0
            for kind, m1, m2, ld in batch:
                if m1 not in self.basis or (kind == "pair" and m2 not in self.basis):
                    continue        # element was replaced during interreduction
                b1 = self.basis[m1]
                if kind == "field":
                    s = self.mul_mono(b1, 1 << m2)
                    sug = self.sugar[m1] + 1
                    self._tag = f"field[{self.origin[m1]}:deg{popcount(m1)}]"
                else:
                    self._tag = f"pair[{self.origin[m1]}:deg{popcount(m1)} x {self.origin[m2]}:deg{popcount(m2)} lcm{ld}]"
                    lcm = m1 | m2
                    s = self.mul_mono(b1, lcm ^ m1) ^ self.mul_mono(self.basis[m2], lcm ^ m2)
                    sug = max(self.sugar[m1] + popcount(lcm ^ m1), self.sugar[m2] + popcount(lcm ^ m2))
                if self.insert(s, sug):
                    new += 1
                    maxlcm = max(maxlcm, ld)
            self.log.append((d, new))
            self._tag = "input"
            if verbose:
                nlin = sum(1 for k in self.basis if popcount(k) == 1)
                lab = "lcmdegree" if self.strategy == "lcm" else "sugar"
                print(f"   {lab}={d} #new={new:4d}  basis={len(self.basis):4d}  linear={nlin}"
                      + (f"  (max true lcm degree among productive pairs: {maxlcm})" if new else ""), flush=True)
        return self.log


if __name__ == "__main__":
    rng = random.Random(2025)
    strategy = "lcm"
    pc = False; ir = True
    args = sys.argv[1:]
    if args and args[0] in ("lcm", "sugar", "label"):
        strategy, args = args[0], args[1:]
    while args and args[0] in ("pc", "noir"):
        if args[0] == "pc": pc = True
        if args[0] == "noir": ir = False
        args = args[1:]
    which = [int(a) for a in args] or [14, 16]
    for n in (12, 14, 16, 18, 20):
        F, E, a6, basis, xQ, sol, _ = build_consistent_instance(n, rng)
        if n not in which:
            continue
        S = DescentSystem(F, a6, xQ, basis)
        N = S.nvars
        T = MonomialTable(N, 5)
        gens = [T.from_terms(q, l, c) for (q, l, c) in S.eqs]
        nsol = count_solutions(N, S.eqs)
        print(f"n={n}: #sol={nsol}")
        t0 = time.time()
        B = Buchberger(T, gens, maxdeg=5, strategy=strategy, product_criterion=pc, interreduce=ir)
        B.run()
        import collections
        print("   parents of linear results:")
        for tag, c in collections.Counter(B.lin_parents).most_common():
            print(f"      {c:3d} x {tag}")
        std = count_standard(N, list(B.basis))
        print(f"   done in {time.time()-t0:.1f}s; #std monomials={std} (#sol={nsol}) -> GB complete: {std == nsol}")
        maxd = max(d for d, new in B.log if new > 0)
        print(f"   max selection key with new polynomials: {maxd}   [strategy={strategy} product_criterion={pc} interreduce={ir}]")
