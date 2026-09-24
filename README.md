# Supplementary material

**The degree of regularity of Weil-descent systems for the third summation polynomial over F_{2^n}:
exact closure dimensions, thresholds and a growth law** — A. Solovyev, 2026.

Everything needed to regenerate every number in the paper: instance generators, the closure
programs, the Gröbner-basis solver patch, all solver inputs and logs, and the raw result tables.

## Layout

```
tools/        programs (Python 3 + numpy; C for wdclose)
data/         raw results: CSV tables, closure outputs, solver logs and inputs
m4gb/         report to the M4GB authors (PR #14) and the benchmark pack sent to them
paper/        paper.tex, paper.pdf
```

## Tools (tools/)

| file | purpose |
|---|---|
| wd_step1.py | GF(2^n) arithmetic (Ben-Or irreducibility test), ordinary curve y²+xy=x³+a₂x²+a₆, S₃ derivation and its check on random triples |
| wd_step2.py | Weil descent with linear constraints (original HKY-style generator); `python3 wd_step2.py` runs self-tests |
| wd_instances.py | corrected generator: uniform P₀,P₁ ∈ F_V, random-Q (HKY) mode, fixed curve/subspace; bilinear solution counter |
| wd_generic_dims.py | generator for unbalanced systems (subspace dims d₀,d₁; n equations) and the (N,n) grid |
| wd_vspace.py | Python closure V_D (HKY Def. 1) in the Boolean ring; first/last fall degree; D_solve via Lemma 2 (standard-monomial count) |
| wd_verify_gb.py | verification of Lemma 2 against the Buchberger–Möller basis of the solution set (n=16) |
| wd_calibrate.py | regenerates the seed-2025 instances n=12…20 run in M4GB (Table 1) |
| wd_experiment.py | the 378-instance survey at n ≤ 24 (data/small_n/results.csv) |
| wd_odd_check.py | odd n: HKY convention (N=n+1) vs balanced (N=n) — Table 4 (D=3 rows) |
| wd_subfield.py | control: V = subfield F_{2^{n/2}} (Section 8) |
| wd_buchberger.py | Buchberger in the Boolean ring; strategies lcm/sugar/label; product criterion switch |
| wd_m4gb_emul.py | emulation of M4GB's update()/selection() criteria as written and corrected (`--fixed`) |
| wd_m4gb_variants.py | permutation / affine-change variants of the n=16 instance (tool-dependence test) |
| wdclose.c | C closure: `gcc -O3 -march=native -fopenmp -o wdclose wdclose.c`; `./wdclose input.txt D [k] [chunk] [macaulay]` |
| wd_export_c.py | writes wdclose input files: `python3 wd_export_c.py N n D [seed]` (seed default 7000+100N+n) |
| m4gb_gm_patch.py, m4gb_gebauer_moeller.patch | the fix of the Gebauer–Möller F-criterion in M4GB's `src/m4gb.hpp` (2018-01-31 revision) |

## Data → tables of the paper

| paper | data |
|---|---|
| Table 1 (n = 12…40) | data/m4gb_runs/n16/m4gb_results_n12_20.txt (before fix), n16/m4gb_traces_session3.txt (after fix); scan1…scan4/results*.txt and logs/ for n = 26…40; data/small_n/results.csv for the closure values at n ≤ 24 |
| §4, 378 instances | data/small_n/results.csv (columns: exp, n, N, seed, mode, nsol, FFD, LFD, D_solve, dim V₂, V₃, V₄) |
| Table 2, grid N = 16…24 | data/small_n/generic_dims.csv; N = 26…30: data/small_n/p4_check.csv, p4_check2.txt |
| Table 3, level 5 | data/closure_level5/c_N*_n*.out (wdclose outputs; `dim(V_5 cap B_<=k)` lines); second instances: scan7/ |
| Table 4, D = 4 out-of-sample (M4GB) | data/m4gb_runs/scan5, scan6 (inputs with predictions in the header, results.txt, logs/) |
| Table 4, D = 5 (M4GB) | data/m4gb_runs/scan7 (32,19/20), scan8 (34, 36, 40) |
| §7.2, stage-1 ranks R₄, R₅ | data/closure_level5/mac_*.out and results_macaulay.txt (wdclose `macaulay` mode) |
| §3.3, M4GB defect | data/m4gb_runs/n16/ (wd_n16.in, full output before the fix, traces after; 20 permutation/affine variants; 10 random-Q instances), m4gb/M4GB_ISSUE.md |
| §8, controls | tools/wd_subfield.py (subfield), random quadrics: generated inline (see paper) |

Solver inputs (`*.in`) are in M4GB's format (`$fieldsize 2`, `$vars N X`, one polynomial per line,
field equations included). Each header line records N, n, seed, the number of solutions and, where a
prediction was made before the run, the predicted degree and margin.

## Reproducing

* Any instance: `python3 tools/wd_export_c.py N n 5` then `./wdclose c_N{N}_n{n}.txt 5` (closure) or
  `./wdclose … 5 10 4096 macaulay` (stage-1 rank). Times: N=24 4 min, N=28 1 h, N=32 5–13 h (8 threads).
* M4GB: build from github.com/cr-marcstevens/m4gb with `tools/m4gb_gm_patch.py` applied to
  `src/m4gb.hpp`; `make bin/solver_m4gb.exe MAXVARS=<N> FIELDSIZE=2`; run `bin/.libs/solver_m4gb.exe -s -i file.in`
  (the file in `bin/` is a libtool wrapper). Read D_reg as the largest `lcmdegree` of a step with `#n > 0`.
* Seeds: square systems of Table 1 use `random.Random(1000*n + k)` (scan1–4) or the seed-2025 sequence
  (n ≤ 20, `wd_calibrate.py`); grid and level-5 systems use `7000 + 100*N + n`; out-of-sample systems
  `9000 + 100*N + n`.

## Licence

Code: MIT. Data and text: CC BY 4.0.
