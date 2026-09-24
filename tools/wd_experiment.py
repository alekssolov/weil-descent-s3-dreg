#!/usr/bin/env python3
"""
wd_experiment.py -- intrinsic-degree survey (wd_vspace) over many instances.

Experiments (all rows appended to results.csv):
  A  the 10 seeds of wd_n16_variability.py (seeds 100..109, wd_step2 construction)
  B  n in NS, K instances, mode=consistent (uniform P0,P1 in F_V), fresh curve+V each
  C  n in NS, K instances, mode=randomQ  (HKY: random nonzero Q), fresh curve+V each
  D  n=16,18: fixed curve+V, vary Q (consistent);  fixed curve, vary V
  E  n=22 (K2 instances, both modes), n=24 (K3 instances, degree cap 4)

CSV columns: exp,n,N,seed,mode,nsol,FFD,LFD,D_solve,dimV2,dimV3,dimV4,time,a6,poly
"""
import csv
import random
import sys
import time

sys.path.insert(0, ".")
from wd_step2 import build_consistent_instance, DescentSystem
from wd_instances import make_instance, count_solutions_bilinear
from wd_vspace import analyze

OUT = "results.csv"
NS = (12, 14, 16, 18, 20)
K = 30
K2, K3 = 10, 4


import os
DONE = set()
if os.path.exists(OUT):
    for r in csv.DictReader(open(OUT)):
        if r["exp"] != "exp":
            DONE.add((r["exp"], int(r["n"]), int(r["seed"])))


def run(exp, n, seed, mode, I, dmax=6, extra=""):
    if (exp, n, seed) in DONE:
        return
    t = time.time()
    nsol = count_solutions_bilinear(I["nvars"], I["nv"], I["eqs"])
    r = analyze(I["nvars"], I["eqs"], dmax=dmax, nsol=nsol)
    d = r["dims"]
    row = [exp, n, I["nvars"], seed, mode, nsol, r["first_fall"], r["last_fall"],
           r["D_solve"] if r["D_solve"] is not None else f">{dmax}",
           d.get(2, ""), d.get(3, ""), d.get(4, ""), f"{time.time()-t:.1f}",
           f"0x{I['a6']:x}", f"0x{I['F'].mod:x}", extra]
    with open(OUT, "a", newline="") as f:
        csv.writer(f).writerow(row)
    print(",".join(map(str, row)), flush=True)


def wrap_step2(n, rng):
    F, E, a6, basis, xQ, sol, _ = build_consistent_instance(n, rng)
    S = DescentSystem(F, a6, xQ, basis)
    return {"F": F, "E": E, "basis": basis, "xQ": xQ, "eqs": S.eqs, "nvars": S.nvars,
            "nv": len(basis), "sol": sol, "a2": E.a2, "a6": a6, "mode": "step2"}


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "ABCDE"
    if not os.path.exists(OUT):
        with open(OUT, "a", newline="") as f:
            csv.writer(f).writerow(["exp", "n", "N", "seed", "mode", "nsol", "FFD", "LFD", "D_solve",
                                    "dimV2", "dimV3", "dimV4", "time", "a6", "poly", "extra"])
    if "A" in which:
        for seed in range(100, 110):
            run("A", 16, seed, "step2", wrap_step2(16, random.Random(seed)))
    if "B" in which:
        for n in NS:
            for k in range(K):
                seed = 1000 * n + k
                run("B", n, seed, "consistent", make_instance(n, random.Random(seed), "consistent"))
    if "C" in which:
        for n in NS:
            for k in range(K):
                seed = 1000 * n + 500 + k
                run("C", n, seed, "randomQ", make_instance(n, random.Random(seed), "randomQ"))
    if "D" in which:
        for n in (16, 18):
            base = make_instance(n, random.Random(77 + n), "consistent")
            curve, basis = (base["a2"], base["a6"]), base["basis"]
            for k in range(10):   # same curve, same V, different Q
                seed = 90000 + 100 * n + k
                run("D_sameCV", n, seed, "consistent",
                    make_instance(n, random.Random(seed), "consistent", curve=curve, basis=basis))
            for k in range(10):   # same curve, different V
                seed = 95000 + 100 * n + k
                run("D_sameC", n, seed, "consistent",
                    make_instance(n, random.Random(seed), "consistent", curve=curve))
    if "E" in which:
        for k in range(K2):
            for mode, off in (("consistent", 0), ("randomQ", 500)):
                seed = 22000 + off + k
                run("E", 22, seed, mode, make_instance(22, random.Random(seed), mode), dmax=4)
        for k in range(K3):
            for mode, off in (("consistent", 0), ("randomQ", 500)):
                seed = 24000 + off + k
                run("E", 24, seed, mode, make_instance(24, random.Random(seed), mode), dmax=4)
