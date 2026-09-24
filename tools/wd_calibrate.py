#!/usr/bin/env python3
"""
wd_calibrate.py -- regenerate the exact random.Random(2025) instances of
wd_m4gb_export.py (same RNG consumption order: n = 12, 14, 16, 18, 20) and
compute their intrinsic degrees with wd_vspace.  Compare with M4GB lcmdegree.
"""
import os
import random
import sys

sys.path.insert(0, ".")
from wd_step2 import build_consistent_instance, DescentSystem
from wd_m4gb_export import poly_to_str
from wd_vspace import analyze

M4GB = {12: 3, 14: 3, 16: 4, 18: 4, 20: 4}   # max lcmdegree with #n>0 (session results)
HKY = {12: 3, 16: 3, 18: 4, 20: 4}           # HKY Table, Magma F4

rng = random.Random(2025)
outdir = "calib_instances"
os.makedirs(outdir, exist_ok=True)
print(f"{'n':>3} {'N':>3} {'#sol':>5} {'FFD':>4} {'LFD':>4} {'D_solve':>8} {'M4GB':>5} {'HKY':>4}  "
      f"{'dim V_2':>7} {'dim V_3':>7} {'dim V_4':>7}  field_poly  a6")
for n in (12, 14, 16, 18, 20):
    F, E, a6, basis, xQ, sol, (P0, P1, R) = build_consistent_instance(n, rng)
    S = DescentSystem(F, a6, xQ, basis)
    with open(os.path.join(outdir, f"wd_n{n}.in"), "w") as f:
        f.write(f"# seed-2025 sequence, n={n}, poly=0x{F.mod:x}, a6=0x{a6:x}\n$fieldsize 2\n$vars {S.nvars} X\n")
        for q, l, c in S.eqs:
            s = poly_to_str(q, l, c)
            if s != "0":
                f.write(s + "\n")
        for i in range(S.nvars):
            f.write(f"X{i}*X{i}+X{i}\n")
    r = analyze(S.nvars, S.eqs, dmax=6)
    d = r["dims"]
    print(f"{n:>3} {S.nvars:>3} {r['nsol']:>5} {r['first_fall']:>4} {r['last_fall']:>4} {r['D_solve']:>8} "
          f"{M4GB.get(n,'-'):>5} {HKY.get(n,'-'):>4}  {d.get(2,'-'):>7} {d.get(3,'-'):>7} {d.get(4,'-'):>7}"
          f"  0x{F.mod:x}  0x{a6:x}   [{r['time']:.1f}s]", flush=True)
