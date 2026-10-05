# Update, October 2026 (paper version 2)

New since the first version of the supplementary material:

## tools/
- `wdrref.c` (v2) -- memory-light closure: echelon in full RREF stored over the free columns only; ~3-9x faster than
  wdclose and with memory going to ~0 as the closure fills the truncated ideal. Used for everything below.
- `wdclose.c` (v10) -- progress, checkpoints/resume, `-nsol`, `-tile`, blocked in-batch elimination (reference tool).
- `wd_s4.py` -- Weil descent of the FOURTH summation polynomial (three-point decompositions): S_4 as the resultant of
  two S_3, self-checks on the curve, descent to degree-6 Boolean equations, EXHAUSTIVE solution count over V^3
  (a count through common roots in F_{2^n} misses solutions with a common root in F_{2^{2n}} and the degenerate
  x1=x2, x3=x(Q); the Boolean ideal sees all of them), export for wdclose/wdrref and for M4GB; `--nocount` option.
- `fastcount.py` -- exact F_2-solution count of the bilinear S_3 systems.

## data/s4/   (paper Section 8, Table 6)
Inputs `s4_N{N}_n{n}.txt` for N=18 (n=6..18), N=21 (n=9,12,15,18,21), N=24 (n=12,24); outputs `*_D7.out` (all collapse
at D=7: dim V_7 = |B_<=7| - #sol) and `*_D8.out` for (18,6),(18,7),(18,8),(21,9) (same standard monomials).
`results_m4gb_s4.txt`: M4GB on the 21-variable S_4 systems ends with std::bad_alloc (matrices of 2.1e6 rows).

## data/table1_closures/   (paper Table 1, closure column for n = 28..40)
Inputs `c_N{N}_n{n}.txt` (seed 7000+100N+n), `nsol.txt`, outputs `*_D4.out` for all eleven systems and `*_D5.out`
for n = 35, 37, 38, 39, 40. Degree-4 closures of the non-collapsing systems equal the P_4 formula to the unit
(61425, 72927, 74727, 85761, 87780 at (36,35), (38,37), (38,38), (40,39), (40,40)).

## data/m4gb_runs/scan8/   (paper Table 4, D = 5 boundary at N = 34, 36, 40)
Inputs with the predicted degree in the header, results of the corrected M4GB.

## m4gb/benchmark_pack_2/
Fourteen S_3 systems with 41..52 variables and the predicted degree of regularity recorded in each header
(sent to the M4GB authors for their next solver).

## paper/
paper.tex / paper.pdf, version 2: cryptanalytic framing, Section 8 on S_4, Table 1 complete up to n = 40.
