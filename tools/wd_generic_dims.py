"""Generic dimensions of V_2, V_3, V_4 for Weil-descent systems with subspace dims (d0, d1):
N = d0 + d1 variables, n equations (n = field degree, d0, d1 <= n)."""
import csv, os, random, sys, time
sys.path.insert(0, '.')
from wd_step1 import GF2n, BinCurve, INF
from wd_step2 import random_subspace
from wd_instances import random_factor_base_point, coords_in_subspace
from wd_vspace import analyze
from wd_odd_check import descent

def count_sol_bilinear(eqs, d0, d1):
    N = d0 + d1
    C, D, E, K = [], [], [], []
    for quad, lin, const in eqs:
        c = [0]*d0; d = 0; e = 0
        for (i, j) in quad:
            assert i < d0 <= j
            c[i] |= 1 << (j - d0)
        for i in lin:
            if i < d0: d |= 1 << i
            else: e |= 1 << (i - d0)
        C.append(c); D.append(d); E.append(e); K.append(const)
    total = 0
    for u in range(1 << d0):
        piv = {}; ok = True
        for r in range(len(eqs)):
            row = E[r]; uu = u; i = 0
            while uu:
                if uu & 1: row ^= C[r][i]
                uu >>= 1; i += 1
            rhs = K[r] ^ (bin(D[r] & u).count('1') & 1)
            while row:
                b = row.bit_length() - 1
                if b in piv:
                    row ^= piv[b][0]; rhs ^= piv[b][1]
                else:
                    piv[b] = (row, rhs); break
            if row == 0 and rhs: ok = False; break
        if ok: total += 1 << (d1 - len(piv))
    return total

def instance(n, d0, d1, rng):
    F = GF2n(n); a2, a6 = F.rand(rng), (F.rand(rng) or 1); E = BinCurve(F, a2, a6)
    while True:
        B0 = random_subspace(F, d0, rng); B1 = random_subspace(F, d1, rng)
        P0 = random_factor_base_point(F, E, B0, rng); P1 = random_factor_base_point(F, E, B1, rng, exclude_x=(P0[0],))
        Q = E.add(P0, P1)
        if Q is not INF: break
    eqs, N = descent(F, a6, Q[0], B0, B1)
    return eqs, N

OUT = 'generic_dims.csv'
done = set()
if os.path.exists(OUT):
    for r in csv.DictReader(open(OUT)):
        done.add((int(r['N']), int(r['n']), int(r['seed'])))
else:
    with open(OUT, 'w', newline='') as f:
        csv.writer(f).writerow(['N', 'd0', 'd1', 'n', 'seed', 'nsol', 'dimV2', 'dimV3', 'dimV3_low', 'dimV4', 'dimV4_low', 'std3', 'std4', 'D_solve', 'time'])
grid = []
for N in (16, 18, 20, 22, 24):
    d0 = N // 2; d1 = N - d0
    for n in range(d1, N + 1):
        grid.append((N, d0, d1, n))
for N, d0, d1, n in grid:
    seed = 7000 + 100 * N + n
    if (N, n, seed) in done: continue
    t = time.time()
    eqs, NN = instance(n, d0, d1, random.Random(seed))
    ns = count_sol_bilinear(eqs, d0, d1)
    r = analyze(NN, eqs, dmax=4, nsol=ns)
    d, dl, s = r['dims'], r['dims_low'], r['std']
    row = [N, d0, d1, n, seed, ns, d.get(2), d.get(3), dl.get(3), d.get(4), dl.get(4), s.get(3), s.get(4), r['D_solve'], f"{time.time()-t:.0f}"]
    with open(OUT, 'a', newline='') as f:
        csv.writer(f).writerow(row)
    print(row, flush=True)
