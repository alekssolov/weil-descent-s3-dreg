M4GB benchmark pack: Weil-descent systems of Semaev's third summation polynomial over GF(2^n)
==========================================================================================

Source: point decomposition x(P0)+x(P1) = x(Q) on an ordinary elliptic curve over GF(2^n), with
x(P0), x(P1) restricted to a random GF(2)-subspace V of dimension ceil(n/2) (the setting of
Huang-Kosters-Yeo, CRYPTO 2015, Sec. 5.2; Petit-Quisquater). After Weil descent this is a system of
n quadratic equations in N = 2*ceil(n/2) variables over GF(2), bilinear in the two coordinate blocks
plus linear terms, together with the field equations X_i*X_i+X_i (included in the files). The files
below use the solver's input format ($fieldsize 2, $vars N X). All systems are consistent (the
number of GF(2) solutions is given; solutions come in pairs by the X0<->X1 symmetry).

Why they may be useful as benchmarks: sparse, highly structured, small degree of regularity (3..6)
at 26..40 variables, so the running time is dominated by many degree-4/5 critical pairs rather
than by a few large matrices -- a different regime from the dense MQ-challenge systems. The exact
maximal step degree is known for every file (independently computed as the smallest D such that the
degree-capped closure V_D of HKY's Definition 1 contains the reduced Groebner basis).

All timings below are wall-clock on one desktop CPU (Cygwin, GCC, -O3 -march=native), with the
build from PR #14 (corrected F-criterion). "D_reg" = largest lcmdegree of a step reporting #n > 0;
the unpatched build reports the same D_reg on all files here except those marked (*), where it
reports one more.

file                vars eqs  #sol  D_reg   time      peak mem   notes
wd_n26_s26000.in     26   26     2    4     26 s       7 MiB    (*) also: the patched build ends with a
                                                                 segmentation fault AFTER the last new
                                                                 polynomials (unpatched build: no crash)
wd_n28_s28000.in     28   28     2    4     71 s       7 MiB
wd_n30_s30000.in     30   30     2    4    214 s       7 MiB
wd_n32_s32000.in     32   32     4    4    470 s       9 MiB
wd_n33_s33000.in     34   33     8    4     18 min   590 MiB    odd n: one extra variable
wd_n34_s34000.in     34   34     4    4     17 min   565 MiB
wd_n35_s35000.in     36   35     2    5     36 min   566 MiB    first degree 5
wd_n36_s36000.in     36   36     2    4     38 min   548 MiB    back to 4 (threshold depends on (N,n))
wd_n37_s37000.in     38   37     4    5     82 min   1.0 GiB
wd_n38_s38000.in     38   38     4    5     78 min   0.8 GiB
wd_n39_s39000.in     40   39     4    5    5.1 h    2.2 GiB
wd_n40_s40000.in     40   40     4    5    6.0 h    1.5 GiB    Kosters-Yeo (Magma, 38 GB) did not finish this size
wd_N30_n24.in        30   24    68    5     12 min   1.5 GiB    fewer equations than variables (subspace dim 15)
wd_N34_n31.in        34   31    15    5     30 min   1.3 GiB
wd_N32_n20.in        32   20  4026    5     52 min   4.4 GiB    large solution set
wd_N32_n19.in        32   19  8154    6     76 min   1.6 GiB    degree 6 needed (and sufficient)

Larger instances (n up to 52, and n = 51 where degree 6 is predicted for a square system) can be
generated on request; the generator is a small Python script and the expected degree follows from
a closed-form criterion we have validated on all of the above.
