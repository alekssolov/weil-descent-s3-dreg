#!/bin/bash
# Stage 1 (minutes each): D=4 for all eleven systems. Expected: collapse for n=28,30,32,33,34,36 (measured D_solve=4);
#   NO collapse for n=35,37,38,39,40 (so D_solve >= 5 there).
declare -A NS=( [c_N28_n28]=3 [c_N30_n30]=4 [c_N32_n32]=1 [c_N34_n33]=1 [c_N34_n34]=1 [c_N36_n35]=2 [c_N36_n36]=2 [c_N38_n37]=2 [c_N38_n38]=1 [c_N40_n39]=4 [c_N40_n40]=3 )
for b in c_N28_n28 c_N30_n30 c_N32_n32 c_N34_n33 c_N34_n34 c_N36_n35 c_N36_n36 c_N38_n37 c_N38_n38 c_N40_n39 c_N40_n40; do
  ./wdrref $b.txt 4 -prog 300 -nsol ${NS[$b]} > ${b}_D4.out 2> ${b}_D4.err
  echo "$b D=4: $(head -1 ${b}_D4.out) | $(grep non-pivot ${b}_D4.out)" | tee -a results_table1.txt
done
# Stage 2 (hours, one at a time, watch peak memory): D=5 for the systems that did not collapse at D=4.
#   (36,35): ~6 GB peak expected; (38,37) and (38,38): ~10 GB -- run alone; skip N=40 (~18 GB).
for b in c_N36_n35 c_N38_n37 c_N38_n38; do
  ./wdrref $b.txt 5 -prog 600 -nsol ${NS[$b]} > ${b}_D5.out 2> ${b}_D5.err
  echo "$b D=5: $(head -1 ${b}_D5.out) | $(grep non-pivot ${b}_D5.out)" | tee -a results_table1.txt
done
