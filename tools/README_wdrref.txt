wdrref v2 -- degree-capped closure with the echelon in FULL reduced row echelon form, stored over the free (non-pivot)
columns only.  Same input and output as wdclose; results verified identical on (20,10)D4, (22,11)D4, (24,12)D4/D5,
S4 (12,12)D7 and S4 (18,8)D7 (62162; 3206/11774/30338; 842 non-pivots).

Why: memory ~ nrows x (#free columns) instead of ~ nrows x M/2, and the rows SHRINK as pivots accumulate, so the work
shrinks too.  Reference environment, one thread: (24,12) D=5: 57 s (wdclose 154 s), S4 (18,8) D=7: 119 s (wdclose 1084 s),
rows memory 0.01-0.03 GB instead of 0.19-0.25 GB.

Build:  gcc -O3 -march=native -fopenmp -o wdrref wdrref.c
Usage:  ./wdrref input.txt [D] [K] [chunk] [-tile T] [-prog SEC] [-nsol S] [-compact FRAC]
   defaults: K=5, chunk=4096, tile=128 (each thread reduces its own tiles with its own tables -- no barriers),
   -compact 0.10 (repack rows when 10% of positions are dead), -nsol S stops when |B_<=D| - S pivots are reached.
   Progress line shows pivots, products, live free columns, rows memory and its peak.
Not yet in wdrref: checkpoints (a Ctrl-C ends the run; rerun from scratch).

Test first (compare with your 1m33 for wdclose tiled):   time ./wdrref s4_N18_n18.txt 7 -nsol 1
Then N=24, D=7 (the wdclose echelon would need ~18 GB; here the peak should stay well below), one at a time:
   ./wdrref s4_N24_n24.txt 7 -prog 300 -nsol 1     > s4_N24_n24_D7.out 2> s4_N24_n24_D7.err
   ./wdrref s4_N24_n12.txt 7 -prog 300 -nsol 3881  > s4_N24_n12_D7.out 2> s4_N24_n12_D7.err     # the decisive n/N = 0.5 point
   ./wdrref s4_N24_n15.txt 7 -prog 300 -nsol 549   > s4_N24_n15_D7.out 2> s4_N24_n15_D7.err
   ./wdrref s4_N24_n18.txt 7 -prog 300 -nsol 53    > s4_N24_n18_D7.out 2> s4_N24_n18_D7.err
Watch "rows=...GB (peak ...)" in the .err; the batch buffers add ~0.6 GB at N=24.

v2 (Sep 27): rows stored in blocks of 4096 (no per-row malloc: Cygwin's allocator serialises threads), and the in-batch
elimination is blocked (groups of 8, Gray table) and parallel. Same results (re-verified on all six checks).
