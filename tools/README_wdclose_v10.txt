wdclose v10 -- same results as v6..v9 (verified: (20,10)D4=4995, (22,11)D4=6864, (24,12)D4=9126, (24,12)D5=51270/9222,
S4 (12,12)D7=3300, Macaulay (24,12)D5=33031), faster and with a barrier-free mode for memory-bound machines.

Changes vs v9:
  * blocked in-batch elimination (groups of 8, Gray-code tables) instead of row-by-row: ~15% faster in the reference environment
  * -tile T : the batch is reduced in tiles of T rows; each thread owns whole tiles and its own table, no barriers at all;
              tile rows stay in the core's cache, pivot rows are streamed once per tile => ~25x less RAM traffic, ~2x more XORs.
              Slower on a compute-bound single core, expected faster on an 8-thread machine that sits at 50% CPU (memory-bound).
              Use with a smaller K (8): tables are per thread, 256 x row_length each.
  * K = 12 is NOT recommended (table leaves the cache: 2x slower at N=24).

Build:  gcc -O3 -march=native -fopenmp -o wdclose wdclose.c
Usage:  ./wdclose input.txt D [K] [chunk] [macaulay] [-every SEC] [-prog SEC] [-ckpt FILE] [-fresh] [-nsol S] [-lowfirst] [-tile T]

Timing tests on s4_N18_n18.txt (your v9 result: 10.5 min wall, 45 min CPU):
  time ./wdclose s4_N18_n18.txt 7 10 4096 -every 0 -nsol 1 -fresh                          # v10 default
  time ./wdclose s4_N18_n18.txt 7 8 4096 -every 0 -nsol 1 -fresh -tile 256                 # tiled, K=8
  time OMP_NUM_THREADS=1 ./wdclose s4_N18_n18.txt 7 10 4096 -every 0 -nsol 1 -fresh        # one thread, for the honest speed-up
If tiled wins, also try -tile 128 and -tile 512. Checkpoints from v7..v9 are readable.
