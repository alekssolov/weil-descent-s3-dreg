"""Odd n: HKY convention (both subspaces of dim ceil(n/2) -> N=n+1) vs balanced ((n+1)/2, (n-1)/2 -> N=n)."""
import random, sys, time
sys.path.insert(0,'.')
from wd_step1 import GF2n, BinCurve, INF
from wd_step2 import random_subspace
from wd_instances import random_factor_base_point, coords_in_subspace
from wd_vspace import analyze

def descent(F, a6, xQ, B0, B1):
    d0, d1 = len(B0), len(B1); N = d0 + d1
    xQ2 = F.sq(xQ)
    # coefficient in F_{2^n} of u_i w_j : v_i^2 v'_j^2 + v_i v'_j xQ ;  u_i : v_i^2 xQ^2 ;  w_j : v'_j^2 xQ^2 ; const a6
    coef = {}
    for i, v in enumerate(B0):
        for j, w in enumerate(B1):
            coef[(i, d0 + j)] = F.add(F.mul(F.sq(v), F.sq(w)), F.mul(F.mul(v, w), xQ))
    lin = {i: F.mul(F.sq(v), xQ2) for i, v in enumerate(B0)}
    lin.update({d0 + j: F.mul(F.sq(w), xQ2) for j, w in enumerate(B1)})
    eqs = []
    for bit in range(F.n):
        q = {k: 1 for k, c in coef.items() if (c >> bit) & 1}
        l = {k: 1 for k, c in lin.items() if (c >> bit) & 1}
        eqs.append((q, l, (a6 >> bit) & 1))
    return eqs, N

def count_sol(N, d0, eqs):
    # brute force (N <= 20 here)
    cnt = 0
    for a in range(1 << N):
        bits = [(a >> i) & 1 for i in range(N)]
        ok = True
        for q, l, c in eqs:
            v = c
            for (i, j) in q: v ^= bits[i] & bits[j]
            for i in l: v ^= bits[i]
            if v: ok = False; break
        cnt += ok
    return cnt

def instance(n, rng, d0, d1):
    F = GF2n(n); a2, a6 = F.rand(rng), (F.rand(rng) or 1); E = BinCurve(F, a2, a6)
    while True:
        B0 = random_subspace(F, d0, rng); B1 = B0 if d1 == d0 else random_subspace(F, d1, rng)
        P0 = random_factor_base_point(F, E, B0, rng); P1 = random_factor_base_point(F, E, B1, rng, exclude_x=(P0[0],))
        Q = E.add(P0, P1)
        if Q is not INF: break
    eqs, N = descent(F, a6, Q[0], B0, B1)
    u = coords_in_subspace(F, B0, P0[0]); w = coords_in_subspace(F, B1, P1[0])
    sol = list(u) + list(w)
    for q, l, c in eqs:
        v = c
        for (i, j) in q: v ^= sol[i] & sol[j]
        for i in l: v ^= sol[i]
        assert v == 0
    return eqs, N, d0

print(f"{'n':>3} {'construction':>12} {'N':>3} {'#sol':>4} {'D_solve':>7} {'LFD':>4}  dimV3  dimV4")
for n in (13, 15, 17, 19):
    for label, d0, d1 in (("HKY ceil", (n+1)//2, (n+1)//2), ("balanced", (n+1)//2, (n-1)//2)):
        for k in range(2):
            eqs, N, d0_ = instance(n, random.Random(500*n + 10*k + (d1 != d0)), d0, d1)
            ns = count_sol(N, d0_, eqs)
            r = analyze(N, eqs, dmax=5, nsol=ns)
            d = r['dims']
            print(f"{n:>3} {label:>12} {N:>3} {ns:>4} {r['D_solve']:>7} {r['last_fall']:>4}  {d.get(3,'-'):>5}  {d.get(4,'-'):>5}", flush=True)
