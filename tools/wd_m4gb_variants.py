#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wd_m4gb_variants.py -- M4GB inputs for a TOOL-dependence test at n=16.

The intrinsic degrees (V_i of HKY Sec. 2, last fall degree, D_solve) are
invariant under affine changes of the variables; F4-type step degrees are not.
We take the seed-2025 n=16 instance (M4GB said max lcmdegree 4, HKY/V_3 say 3)
and emit

   m4gb_variants/wd_n16_perm_<k>.in     random permutation of the variables
   m4gb_variants/wd_n16_aff_<k>.in      random invertible affine change  x = A y + b
   m4gb_variants/wd_n16_randQ_<k>.in    HKY-literal instances (random point Q), fresh curve/V

Each file also contains the known solution (if any) and the intrinsic degrees
computed here, as comments, so the M4GB result can be compared line by line.

Cygwin loop (after building bin/solver_m4gb_n16_gf2 with MAXVARS=16 FIELDSIZE=2):
  for f in m4gb_variants/*.in; do
    echo -n "$f: "; ./bin/solver_m4gb_n16_gf2 -s -i $f 2>&1 | grep -E "lcmdegree|#n:" | paste -sd' '
  done
Read D_reg as the largest lcmdegree of a step with #n > 0 (HKY's definition).
"""
import os
import random
import sys

sys.path.insert(0, ".")
from wd_step2 import build_consistent_instance, DescentSystem
from wd_instances import make_instance, count_solutions_bilinear
from wd_vspace import analyze

OUT = "m4gb_variants"
os.makedirs(OUT, exist_ok=True)


def write_m4gb(path, eqs, N, header, sol=None):
    lines = [f"# {h}" for h in header] + ["$fieldsize 2", f"$vars {N} X"]
    for quad, lin, const in eqs:
        terms = [f"X{i}*X{j}" for (i, j) in sorted(quad)] + [f"X{i}" for i in sorted(lin)] + (["1"] if const else [])
        if terms:
            lines.append("+".join(terms))
    lines += [f"X{i}*X{i}+X{i}" for i in range(N)]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    if sol is not None:
        with open(path[:-3] + ".ans", "w") as f:
            f.write(" ".join(map(str, sol)) + "\n")


def evaluate(eqs, bits):
    out = []
    for quad, lin, const in eqs:
        v = const
        for (i, j) in quad:
            v ^= bits[i] & bits[j]
        for i in lin:
            v ^= bits[i]
        out.append(v)
    return out


def permute(eqs, perm):
    """x_i -> x_{perm[i]}."""
    new = []
    for quad, lin, const in eqs:
        q = {}
        for (i, j) in quad:
            a, b = sorted((perm[i], perm[j]))
            q[(a, b)] = q.get((a, b), 0) ^ 1
        l = {}
        for i in lin:
            l[perm[i]] = l.get(perm[i], 0) ^ 1
        new.append(({k: 1 for k, v in q.items() if v}, {k: 1 for k, v in l.items() if v}, const))
    return new


def affine(eqs, N, A, b):
    """Substitute x_i = sum_j A[i][j] y_j + b[i]  (A rows as bitmasks over j)."""
    new = []
    for quad, lin, const in eqs:
        q = {}      # unordered pairs (a<b)
        l = {}      # linear
        c = const
        def add_lin(a):
            l[a] = l.get(a, 0) ^ 1
        def add_quad(a, bb):
            if a == bb:
                add_lin(a)          # y_a^2 = y_a
            else:
                k = (min(a, bb), max(a, bb))
                q[k] = q.get(k, 0) ^ 1
        for (i, j) in quad:
            Li, ci = A[i], b[i]
            Lj, cj = A[j], b[j]
            ii = [a for a in range(N) if Li >> a & 1]
            jj = [a for a in range(N) if Lj >> a & 1]
            for a in ii:
                for bb in jj:
                    add_quad(a, bb)
            if cj:
                for a in ii:
                    add_lin(a)
            if ci:
                for a in jj:
                    add_lin(a)
            c ^= ci & cj
        for i in lin:
            for a in range(N):
                if A[i] >> a & 1:
                    add_lin(a)
            c ^= b[i]
        new.append(({k: 1 for k, v in q.items() if v}, {k: 1 for k, v in l.items() if v}, c))
    return new


def random_gl(N, rng):
    while True:
        A = [rng.randrange(1, 1 << N) for _ in range(N)]
        # invertible?  rank via elimination
        piv = {}
        ok = True
        for row in A:
            v = row
            while v:
                t = v.bit_length() - 1
                if t in piv:
                    v ^= piv[t]
                else:
                    piv[t] = v
                    break
            if v == 0:
                ok = False
                break
        if ok:
            return A


def solve_affine_inverse(A, b, x, N):
    """Find y with A y + b = x (A invertible)."""
    target = x ^ b_to_mask(b)
    # solve sum_j A[i][j] y_j = target_i for all i: matrix rows A[i] over j
    rows = [(A[i], (target >> i) & 1) for i in range(N)]
    piv = {}
    for row, rhs in rows:
        v, r = row, rhs
        while v:
            t = v.bit_length() - 1
            if t in piv:
                pv, pr = piv[t]
                v ^= pv; r ^= pr
            else:
                piv[t] = (v, r)
                break
    # back substitution
    y = 0
    for t in sorted(piv):
        v, r = piv[t]
        rest = v & ~(1 << t)
        val = r
        vv = rest
        while vv:
            s = vv.bit_length() - 1
            val ^= (y >> s) & 1
            vv ^= 1 << s
        y |= val << t
    return y


def b_to_mask(b):
    m = 0
    for i, v in enumerate(b):
        if v:
            m |= 1 << i
    return m


if __name__ == "__main__":
    # ---- the seed-2025 n=16 instance (RNG consumed as in wd_m4gb_export: n=12, 14 first)
    rng = random.Random(2025)
    for n in (12, 14, 16):
        F, E, a6, basis, xQ, sol, _ = build_consistent_instance(n, rng)
    S = DescentSystem(F, a6, xQ, basis)
    N = S.nvars
    eqs0 = S.eqs
    base = analyze(N, eqs0, dmax=6)
    print(f"base instance n=16 seed-2025: #sol={base['nsol']} D_solve={base['D_solve']} LFD={base['last_fall']}")
    vrng = random.Random(4242)
    x_sol = list(sol)
    for k in range(10):
        perm = list(range(N)); vrng.shuffle(perm)
        eqs = permute(eqs0, perm)
        ysol = [0] * N
        for i in range(N):
            ysol[perm[i]] = x_sol[i]
        assert not any(evaluate(eqs, ysol))
        r = analyze(N, eqs, dmax=6)
        assert (r["D_solve"], r["last_fall"], r["dims"]) == (base["D_solve"], base["last_fall"], base["dims"])
        write_m4gb(f"{OUT}/wd_n16_perm_{k}.in", eqs, N,
                   [f"seed-2025 n=16 instance, variables permuted: x_i -> X{perm}",
                    f"intrinsic: #sol={r['nsol']} D_solve={r['D_solve']} first_fall={r['first_fall']} last_fall={r['last_fall']}"], ysol)
    for k in range(10):
        A = random_gl(N, vrng)
        b = [vrng.randrange(2) for _ in range(N)]
        eqs = affine(eqs0, N, A, b)
        # y with x = A y + b ; check the transformed system on it and on random points
        xm = b_to_mask(x_sol)
        y = solve_affine_inverse(A, b, xm, N)
        ybits = [(y >> i) & 1 for i in range(N)]
        assert not any(evaluate(eqs, ybits))
        for _ in range(50):     # transformed(y) == original(A y + b) on random y
            yy = vrng.randrange(1 << N)
            yb = [(yy >> i) & 1 for i in range(N)]
            xb = [(bin(A[i] & yy).count("1") & 1) ^ b[i] for i in range(N)]
            assert evaluate(eqs, yb) == evaluate(eqs0, xb)
        r = analyze(N, eqs, dmax=6)
        assert (r["D_solve"], r["last_fall"], r["dims"]) == (base["D_solve"], base["last_fall"], base["dims"]), r
        write_m4gb(f"{OUT}/wd_n16_aff_{k}.in", eqs, N,
                   [f"seed-2025 n=16 instance after random affine change x = A y + b (dense)",
                    f"intrinsic (invariant): #sol={r['nsol']} D_solve={r['D_solve']} first_fall={r['first_fall']} last_fall={r['last_fall']}"], ybits)
    print("perm/aff variants written; intrinsic degrees identical on all 20 (checked)")
    # ---- HKY-literal random-Q instances
    for k in range(10):
        I = make_instance(16, random.Random(16500 + k), "randomQ")
        ns = count_solutions_bilinear(I["nvars"], I["nv"], I["eqs"])
        r = analyze(I["nvars"], I["eqs"], dmax=6, nsol=ns)
        write_m4gb(f"{OUT}/wd_n16_randQ_{k}.in", I["eqs"], I["nvars"],
                   [f"HKY-literal: random nonzero Q, n=16, seed={16500+k}, poly=0x{I['F'].mod:x}, a6=0x{I['a6']:x}",
                    f"intrinsic: #sol={ns} D_solve={r['D_solve']} first_fall={r['first_fall']} last_fall={r['last_fall']}"])
        print(f"  randQ {k}: #sol={ns} D_solve={r['D_solve']} LFD={r['last_fall']}")
    print("done:", sorted(os.listdir(OUT)))
