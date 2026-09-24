"""Does the choice of V matter?  V = subfield F_{2^{n/2}} (Frobenius-invariant) vs random subspace."""
import random, sys
sys.path.insert(0,'.')
from wd_step1 import GF2n, BinCurve, INF
from wd_instances import random_factor_base_point, coords_in_subspace
from wd_odd_check import descent
from wd_generic_dims import count_sol_bilinear
from wd_vspace import analyze
from math import comb

def subfield_basis(F, k):
    """basis of {x : x^(2^k) = x} in F_{2^n} (n = F.n, k | n) via kernel of x -> x^(2^k)+x."""
    n=F.n
    # matrix of the F_2-linear map L(x)=x^(2^k)+x on the polynomial basis
    cols=[]
    for i in range(n):
        x=1<<i; y=x
        for _ in range(k): y=F.sq(y)
        cols.append(y ^ x)
    # kernel: solve sum c_i cols[i] = 0
    # Gaussian elimination on columns
    piv={}; basis=[]; tags=[]
    for i,c in enumerate(cols):
        v=c; t=1<<i
        while v:
            b=v.bit_length()-1
            if b in piv: pv,pt=piv[b]; v^=pv; t^=pt
            else: piv[b]=(v,t); break
        if v==0: basis.append(t)   # kernel vector (as field element: sum of basis elements with bits in t)
    assert len(basis)==k, (len(basis),k)
    return basis

def instance_V(n, rng, mode):
    F=GF2n(n); a2,a6=F.rand(rng),(F.rand(rng) or 1); E=BinCurve(F,a2,a6)
    d=n//2
    if mode=="subfield": B=subfield_basis(F,d)
    else:
        from wd_step2 import random_subspace
        B=random_subspace(F,d,rng)
    while True:
        P0=random_factor_base_point(F,E,B,rng); P1=random_factor_base_point(F,E,B,rng,exclude_x=(P0[0],))
        Q=E.add(P0,P1)
        if Q is not INF: break
    eqs,N=descent(F,a6,Q[0],B,B)
    return eqs,N,d

print(f"{'n':>3} {'V':>9} {'#sol':>5} {'dimV2':>5} {'dimV3':>6} {'dimV3∩B2':>8} {'dimV4':>6} {'D_solve':>7}   (P3 random-V formula)")
for n in (16,18,20,22,24):
    for mode in ("random","subfield"):
        for k in range(2):
            eqs,N,d=instance_V(n,random.Random(300+10*n+k),mode)
            ns=count_sol_bilinear(eqs,d,d)
            r=analyze(N,eqs,dmax=4,nsol=ns); dm=r['dims']; dl=r['dims_low']
            P3=(3*N-2)*n-comb(N,2)
            print(f"{n:>3} {mode:>9} {ns:>5} {dm.get(2):>5} {dm.get(3):>6} {dl.get(3):>8} {str(dm.get(4)):>6} {str(r['D_solve']):>7}   P3={P3}", flush=True)
