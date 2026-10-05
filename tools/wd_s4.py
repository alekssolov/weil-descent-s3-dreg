"""Weil descent of the fourth summation polynomial S_4 (three-point decomposition) over F_{2^n}.

S_4(X1,X2,X3,X4) = Res_X( S_3(X1,X2,X), S_3(X3,X4,X) ), with S_3(X1,X2,X) = A X^2 + B X + C,
  A = (X1+X2)^2,  B = X1 X2,  C = (X1 X2)^2 + a6      (char 2, a1 = 1)
so that  S_4 = (A C' + A' C)^2 + (A B' + A' B)(B C' + B' C)    (resultant of two quadratics; signs vanish in char 2).

Decomposition Q = P1 + P2 + P3 with x(P_i) in a random F_2-subspace V of dimension d:  S_4(x1,x2,x3,x(Q)) = 0,
x_i = sum_j u_{ij} v_j.  Squaring is F_2-linear, so the descended equations have F_2-degree <= 6.

usage:  python3 wd_s4.py n [seed] [d]      (default d = ceil(n/3), N = 3d variables, n equations)
writes  s4_N{N}_n{n}.txt  in the wdclose input format (monomials as hex bitmasks, one equation per line) and prints statistics.
"""
import sys, random, time, itertools
sys.path.insert(0, '.')
from wd_step1 import GF2n, BinCurve, INF, S3_char2

# ---------------------------------------------------------------- field helpers
def solve_z2z(F, dval):
    """all z with z^2 + z = dval (0 or 2 solutions) by linear algebra over F_2."""
    n = F.n
    cols = [F.add(F.sq(1 << i), 1 << i) for i in range(n)]           # image of basis vectors under z -> z^2+z
    # gaussian elimination on columns: find one preimage, kernel is {0,1}
    rows = []                                                          # (value, combination) pivot list
    piv = {}
    for i, c in enumerate(cols):
        v, comb = c, 1 << i
        for b in sorted(piv, reverse=True):
            if v >> b & 1: pv, pc = piv[b]; v ^= pv; comb ^= pc
        if v: piv[v.bit_length() - 1] = (v, comb)
    v, comb = dval, 0
    for b in sorted(piv, reverse=True):
        if v >> b & 1: pv, pc = piv[b]; v ^= pv; comb ^= pc
    if v: return []
    return [comb, comb ^ 1]

def roots_quadratic(F, A, B, C):
    """roots X of A X^2 + B X + C over F_{2^n}."""
    if A == 0:
        if B == 0: return []          # constant
        return [F.mul(C, F.inv(B))]
    if B == 0:                        # X^2 = C/A, unique square root
        return [F.pow(F.mul(C, F.inv(A)), 1 << (F.n - 1))]
    # X = (B/A) Y :  Y^2 + Y = A C / B^2
    t = F.mul(B, F.inv(A))
    dval = F.mul(F.mul(A, C), F.inv(F.sq(B)))
    return [F.mul(t, y) for y in solve_z2z(F, dval)]

def S4_char2(F, a6, x1, x2, x3, x4):
    A = F.sq(F.add(x1, x2)); B = F.mul(x1, x2); C = F.add(F.sq(B), a6)
    A2 = F.sq(F.add(x3, x4)); B2 = F.mul(x3, x4); C2 = F.add(F.sq(B2), a6)
    t1 = F.add(F.mul(A, C2), F.mul(A2, C))
    t2 = F.add(F.mul(A, B2), F.mul(A2, B))
    t3 = F.add(F.mul(B, C2), F.mul(B2, C))
    return F.add(F.sq(t1), F.mul(t2, t3))

# ---------------------------------------------------------------- Boolean polynomials with F_{2^n} coefficients
class BP:
    """dict: monomial bitmask (over the N F_2-variables) -> coefficient in F_{2^n}; Boolean ring (u^2 = u)."""
    __slots__ = ('F', 'c')
    def __init__(self, F, c=None): self.F = F; self.c = c if c is not None else {}
    def add(self, o):
        r = dict(self.c)
        for m, v in o.c.items():
            w = r.get(m, 0) ^ v
            if w: r[m] = w
            elif m in r: del r[m]
        return BP(self.F, r)
    def mul(self, o):
        F = self.F; r = {}
        for m1, v1 in self.c.items():
            for m2, v2 in o.c.items():
                m = m1 | m2; w = r.get(m, 0) ^ F.mul(v1, v2)
                if w: r[m] = w
                elif m in r: del r[m]
        return BP(F, r)
    def sq(self):     # (sum c_m m)^2 = sum c_m^2 m  in the Boolean ring, char 2
        F = self.F
        return BP(F, {m: F.sq(v) for m, v in self.c.items()})
    def maxdeg(self): return max((bin(m).count('1') for m in self.c), default=0)

def descent_S4(F, a6, q, bases):
    """bases: list of three lists of d field elements (subspace bases for x1,x2,x3). Returns (eqs, N):
    eqs = list of n Boolean polynomials (dict bitmask -> 1) = coordinates of S_4(x1,x2,x3,q)."""
    d = [len(b) for b in bases]; N = sum(d); off = [0, d[0], d[0] + d[1]]
    def lin(k):   # x_k = sum_j u_{kj} v_j
        return BP(F, {1 << (off[k] + j): v for j, v in enumerate(bases[k])})
    x1, x2, x3 = lin(0), lin(1), lin(2)
    const = lambda v: BP(F, {0: v} if v else {})
    A = x1.add(x2).sq(); B = x1.mul(x2); C = B.sq().add(const(a6))
    x4 = const(q)
    A2 = x3.add(x4).sq(); B2 = x3.mul(x4); C2 = B2.sq().add(const(a6))
    t1 = A.mul(C2).add(A2.mul(C))
    t2 = A.mul(B2).add(A2.mul(B))
    t3 = B.mul(C2).add(B2.mul(C))
    S = t1.sq().add(t2.mul(t3))
    n = F.n
    eqs = [dict() for _ in range(n)]
    for m, v in S.c.items():
        for r in range(n):
            if v >> r & 1: eqs[r][m] = 1
    return eqs, N, S

def random_subspace(F, d, rng):
    while True:
        vs = [F.rand(rng) for _ in range(d)]
        # independence over F_2
        piv = {}
        ok = True
        for v in vs:
            w = v
            for b in sorted(piv, reverse=True):
                if w >> b & 1: w ^= piv[b]
            if w == 0: ok = False; break
            piv[w.bit_length() - 1] = w
        if ok: return vs

def coords(vs, x):
    """coordinates of x in span(vs) over F_2 (or None)."""
    piv = {}                              # bit -> (vector, combination)
    for i, v in enumerate(vs):
        w, comb = v, 1 << i
        for b in sorted(piv, reverse=True):
            if w >> b & 1: pv, pc = piv[b]; w ^= pv; comb ^= pc
        piv[w.bit_length() - 1] = (w, comb)
    w, comb = x, 0
    for b in sorted(piv, reverse=True):
        if w >> b & 1: pv, pc = piv[b]; w ^= pv; comb ^= pc
    return comb if w == 0 else None

def span_set(vs):
    S = {0}
    for v in vs: S |= {s ^ v for s in S}
    return S

def count_solutions(F, a6, q, bases):
    """exact number of (x1,x2,x3) in V1 x V2 x V3 with S_4(x1,x2,x3,q) = 0, by exhaustive evaluation.
    (A count through common roots of the two S_3 in F_{2^n} undercounts: the resultant also vanishes when the common
    root lies in F_{2^{2n}} -- solutions of the descended system without a decomposition behind them -- and when
    x1 = x2 and x3 = q. The Boolean ideal sees all of them, so Lemma 2 needs this count.)"""
    V1, V2, V3 = (sorted(span_set(b)) for b in bases)
    sols = []
    for x1 in V1:
        for x2 in V2:
            for x3 in V3:
                if S4_char2(F, a6, x1, x2, x3, q) == 0: sols.append((x1, x2, x3))
    return sols

def evaluate(eq, bits):
    """Boolean polynomial (dict bitmask -> 1) at a 0/1 assignment given as a bitmask."""
    return sum(1 for m in eq if (m & bits) == m) & 1

def export(eqs, N, D, path):
    with open(path, 'w') as f:
        f.write(f"{N} {D}\n{len(eqs)}\n")
        for e in eqs:
            f.write(" ".join(format(m, 'x') for m in sorted(e)) + "\n")

def instance(n, d, rng, verify=True, count=True):
    F = GF2n(n); a2, a6 = F.rand(rng), (F.rand(rng) or 1); E = BinCurve(F, a2, a6)
    if verify:   # S_4 vanishes on x-coordinates of P1+P2+P3 = Q, and (generically) not on random quadruples
        for _ in range(5):
            P = [E.random_point(rng) for _ in range(3)]
            Q = E.add(E.add(P[0], P[1]), P[2])
            if Q is INF: continue
            assert S4_char2(F, a6, P[0][0], P[1][0], P[2][0], Q[0]) == 0, "S_4 does not vanish on a decomposition"
        nz = sum(1 for _ in range(5) if S4_char2(F, a6, F.rand(rng), F.rand(rng), F.rand(rng), F.rand(rng)) != 0)
        assert nz >= 3, "S_4 vanishes on random quadruples"
    while True:
        bases = [random_subspace(F, d, rng) for _ in range(3)]
        # consistent instance: Q = P1 + P2 + P3 with x(P_i) in V_i
        pts = []
        for k in range(3):
            V = list(span_set(bases[k]))
            while True:
                x = rng.choice(V); P = E_point_with_x(E, x, rng)
                if P is not None: pts.append(P); break
        Q = E.add(E.add(pts[0], pts[1]), pts[2])
        if Q is not INF and Q[0] != 0: break
    q = Q[0]
    eqs, N, S = descent_S4(F, a6, q, bases)
    sols = count_solutions(F, a6, q, bases) if count else []
    if verify and count:
        for (x1, x2, x3) in sols:
            bits = coords(bases[0], x1) | (coords(bases[1], x2) << d) | (coords(bases[2], x3) << (2 * d))
            assert all(evaluate(e, bits) == 0 for e in eqs), "descended system does not vanish on a solution"
        # and the known decomposition is among the solutions
        assert (pts[0][0], pts[1][0], pts[2][0]) in set(sols)
    return F, E, q, bases, eqs, N, S, sols

def E_point_with_x(E, x, rng):
    F = E.F
    # y^2 + x y = x^3 + a2 x^2 + a6 ;  x != 0: y = x z, z^2 + z = (x^3 + a2 x^2 + a6)/x^2
    if x == 0: return None
    rhs = F.add(F.add(F.mul(F.sq(x), x), F.mul(E.a2, F.sq(x))), E.a6)
    zs = solve_z2z(F, F.mul(rhs, F.inv(F.sq(x))))
    if not zs: return None
    z = rng.choice(zs)
    return (x, F.mul(x, z))

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--nocount"]; nocount = "--nocount" in sys.argv
    n = int(args[0]); seed = int(args[1]) if len(args) > 1 else 11000 + n
    d = int(args[2]) if len(args) > 2 else -(-n // 3)
    t = time.time()
    F, E, q, bases, eqs, N, S, sols = instance(n, d, random.Random(seed), count=not nocount)
    degs = [max((bin(m).count('1') for m in e), default=0) for e in eqs]
    from collections import Counter
    cnt = Counter(bin(m).count('1') for m in S.c)
    print(f"n={n} d={d} N={N} seed={seed}: {len(eqs)} equations, max degree {max(degs)} (per-equation degrees {sorted(Counter(degs).items())})")
    print(f"  monomials of S_4 over F_2^n by degree: {sorted(cnt.items())}   total {len(S.c)}")
    print(f"  solutions (x1,x2,x3) in V^3: {len(sols) if not nocount else 'not counted (--nocount)'}   time {time.time()-t:.1f}s")
    D = max(degs)
    path = f"s4_N{N}_n{n}.txt"; export(eqs, N, D, path); print(f"  written {path} (run wdclose with D >= {D})")
