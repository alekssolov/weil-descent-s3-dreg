scan8: direct M4GB test of the degree-5/6 boundary predicted by Q5(N,n) = 3n*C(N,3) - C(3n+1,2)*N - n(79N-660)
Degree 5 suffices (D_reg=5) iff Q5 > C(N,5); otherwise D_reg >= 6 (expect productive lcmdegree=6 steps).
(#sol not computed for N>=40: ~2^(N-n) up to a small factor)

wd_N34_n22.in  Q5=  275198  C(N,5)=  278256  margin   -3058 (-1.10%)  -> D_reg = 6
wd_N34_n23.in  Q5=  284188  C(N,5)=  278256  margin   +5932 (+2.13%)  -> D_reg = 5
wd_N36_n24.in  Q5=  367056  C(N,5)=  376992  margin   -9936 (-2.64%)  -> D_reg = 6
wd_N36_n25.in  Q5=  378300  C(N,5)=  376992  margin   +1308 (+0.35%)  -> D_reg = 5
wd_N40_n30.in  Q5=  650400  C(N,5)=  658008  margin   -7608 (-1.16%)  -> D_reg = 6
wd_N40_n31.in  Q5=  666500  C(N,5)=  658008  margin   +8492 (+1.29%)  -> D_reg = 5
wd_N44_n36.in  Q5= 1069992  C(N,5)= 1086008  margin  -16016 (-1.47%)  -> D_reg = 6
wd_N44_n37.in  Q5= 1092388  C(N,5)= 1086008  margin   +6380 (+0.59%)  -> D_reg = 5

RESULTS (patched M4GB, real binaries bin/.libs/m4gb_nXX_fixed.exe):
  wd_N34_n23 -> 5 (1.8 h, 5.6 GB)      wd_N34_n22 -> 6 (3.4 h, 8.1 GB; 6 productive degree-6 steps)
  wd_N36_n25 -> 5 (3.5 h, 6.9 GB)      wd_N36_n24 -> 6 (5.6 h, 7.8 GB; 26 productive degree-6 steps)
  wd_N40_n31 -> 5 (27.8 h, 3.6 GB)     wd_N40_n30 -> 6 (62.8 h, 8.5 GB; 15 productive degree-6 steps)
All six as predicted. N=44 not run (weeks; degree-6 phase likely > 16 GB).
NOTE (reproducibility): bin/solver_m4gb.exe is a 27 KB libtool wrapper; the real binary is bin/.libs/solver_m4gb.exe.
Renaming the wrapper does not preserve a build -- copy bin/.libs/solver_m4gb.exe after every make.
