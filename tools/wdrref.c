/*
 * wdrref.c -- degree-capped closure V_D (HKY Def. 1) with a memory-light echelon.
 *
 * Same input, same output as wdclose.c.  The echelon is kept in FULL reduced row echelon form:
 * every stored row has its pivot bit and otherwise only bits at free (non-pivot) columns, and the
 * rows are stored as dense bit vectors over the free columns only (compacted lazily).  Memory is
 * therefore ~ nrows * (#free columns) bits and goes to ~0 when the closure fills the truncated ideal,
 * instead of ~ nrows * M / 2 bits for the block-triangular echelon of wdclose.
 *
 * Products are formed in the full monomial space, reduced against the RREF rows (one pass: pivot rows
 * have no other pivot bits), self-eliminated, and the new pivots are back-substituted into all stored
 * rows.  Gray-code tables (Four Russians) are used in both directions; each thread reduces its own
 * tile of rows with its own tables (no barriers).
 *
 * Build:  gcc -O3 -march=native -fopenmp -o wdrref wdrref.c
 * Usage:  ./wdrref input.txt [D] [K] [chunk] [-tile T] [-prog SEC] [-nsol S] [-compact FRAC]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#else
static int omp_get_max_threads(void) { return 1; }
static int omp_get_thread_num(void) { return 0; }
#endif
typedef uint64_t u64;
static double now(void){ struct timespec ts; clock_gettime(CLOCK_REALTIME,&ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }
static void fmt_time(double s, char *buf, size_t n) { long t = (long)(s + 0.5); snprintf(buf, n, "%ld:%02ld:%02ld", t / 3600, (t / 60) % 60, t % 60); }
static volatile sig_atomic_t stop_req = 0;
static void on_signal(int sig) { (void)sig; if (stop_req) _exit(130); stop_req = 1; }

/* ---------------- monomials ---------------- */
static int N, D, M, WM;              /* M = #monomials of degree <= D, WM = words for a full-space vector */
static int *deg_start; static u64 *mono_mask; static uint8_t *mono_deg; static u64 binom[64][10]; static int32_t *mulvar;
static int col_index(u64 mask) { int d = __builtin_popcountll(mask), j = 1; u64 r = 0; u64 m = mask; while (m) { int v = __builtin_ctzll(m); r += binom[v][j]; j++; m &= m - 1; } return deg_start[d] + (int)r; }

/* ---------------- free-column index space ---------------- */
/* free_pos[c] = position of column c in the free index space, or -1 if c is a pivot (or dead). Positions of pivot
   columns are "dead" until the next compaction (all stored rows have zero there after back-substitution). */
static int *free_pos = NULL; static int *pos_col = NULL;   /* position -> column (or -1 if dead) */
static int NP = 0, WF = 0;            /* number of positions in use (incl. dead), words per stored row */
static long ndead = 0;
static char *is_pivot = NULL;

/* ---------------- echelon ---------------- */
static int *pivcol = NULL; static char *expanded = NULL; static long nrows = 0, cap_rows = 0;
#define RB 4096                       /* rows per storage block; all blocks share the current WF */
static u64 **blocks = NULL; static long nblocks = 0, cap_blocks = 0;
#define ROWP(i) (blocks[(i) >> 12] + (size_t)((i) & (RB - 1)) * WF)
static void ensure_blocks(long need_rows) {
    long nb = (need_rows + RB - 1) / RB;
    if (nb > cap_blocks) { long nc = cap_blocks ? cap_blocks : 64; while (nc < nb) nc *= 2; blocks = realloc(blocks, nc * sizeof(u64 *)); cap_blocks = nc; }
    for (; nblocks < nb; ++nblocks) { blocks[nblocks] = malloc((size_t)RB * (WF ? WF : 1) * sizeof(u64)); if (!blocks[nblocks]) { fprintf(stderr, "out of memory (rows)\n"); exit(1); } }
}
static long *piv_of_col = NULL;      /* column -> row index (or -1) */
static size_t mem_rows = 0, mem_peak = 0;
static void ensure_rows(long need) {
    ensure_blocks(need);
    if (need <= cap_rows) return; long nc = cap_rows ? cap_rows : 65536; while (nc < need) nc *= 2;
    pivcol = realloc(pivcol, nc * sizeof(int)); expanded = realloc(expanded, nc); cap_rows = nc;
}
static inline int getbit(const u64 *r, int c) { return (int)((r[c >> 6] >> (c & 63)) & 1); }
static inline void setbit(u64 *r, int c) { r[c >> 6] ^= 1ull << (c & 63); }
static inline int highest_bit(const u64 *r, int w) { for (int i = w - 1; i >= 0; --i) if (r[i]) return i * 64 + 63 - __builtin_clzll(r[i]); return -1; }

/* compaction: drop dead positions, renumber, repack every row */
static void compact(void) {
    int *newpos = malloc(NP * sizeof(int)); int np = 0;
    for (int p = 0; p < NP; ++p) { if (pos_col[p] >= 0) { newpos[p] = np; pos_col[np] = pos_col[p]; free_pos[pos_col[p]] = np; np++; } else newpos[p] = -1; }
    int newWF = (np + 63) / 64; if (newWF == 0) newWF = 1;
    long nbl = (nrows + RB - 1) / RB;
    #pragma omp parallel for schedule(dynamic, 1)
    for (long b = 0; b < nbl; ++b) {
        u64 *nbk = calloc((size_t)RB * newWF, sizeof(u64));
        long i0 = b * RB, i1 = i0 + RB < nrows ? i0 + RB : nrows;
        for (long i = i0; i < i1; ++i) {
            const u64 *r = blocks[b] + (size_t)(i - i0) * WF; u64 *nr = nbk + (size_t)(i - i0) * newWF;
            for (int w = 0; w < WF; ++w) { u64 x = r[w]; while (x) { int bb = __builtin_ctzll(x); x &= x - 1; int q = newpos[w * 64 + bb]; if (q >= 0) nr[q >> 6] |= 1ull << (q & 63); } }
        }
        free(blocks[b]); blocks[b] = nbk;
    }
    for (long b = nbl; b < nblocks; ++b) { free(blocks[b]); }   /* spare blocks are re-created with the new WF on demand */
    nblocks = nbl;
    mem_rows = (size_t)nrows * newWF * 8;
    NP = np; WF = newWF; ndead = 0; free(newpos);
}

/* ---------------- Gray-code table over a group of stored rows ---------------- */
/* table rows are over the free space (WF words); group = rows g[0..kk) */
static void build_table(u64 *tab, const long *g, int kk, int L) {
    memset(tab, 0, (size_t)L * sizeof(u64));
    for (int j = 0; j < kk; ++j) { const u64 *rj = ROWP(g[j]); int half = 1 << j;
        for (int p = 0; p < half; ++p) { u64 *dst = tab + (size_t)(p | half) * WF; const u64 *src = tab + (size_t)p * WF; memcpy(dst, src, (size_t)L * sizeof(u64)); for (int w = 0; w < L; ++w) dst[w] ^= rj[w]; } }
}

int main(int argc, char **argv) {
    const char *pos[8] = {0}; int npos = 0; int TILE = 128, K = 5, CH = 4096; double prog_every = 30, compact_frac = 0.10; long nsol = -1;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-tile") && i + 1 < argc) TILE = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-prog") && i + 1 < argc) prog_every = atof(argv[++i]);
        else if (!strcmp(argv[i], "-nsol") && i + 1 < argc) nsol = atol(argv[++i]);
        else if (!strcmp(argv[i], "-compact") && i + 1 < argc) compact_frac = atof(argv[++i]);
        else if (npos < 8) pos[npos++] = argv[i];
    }
    if (npos < 1) { fprintf(stderr, "usage: %s input.txt [D] [K] [chunk] [-tile T] [-prog SEC] [-nsol S] [-compact FRAC]\n", argv[0]); return 1; }
    FILE *f = fopen(pos[0], "r"); if (!f) { perror("input"); return 1; }
    int Dfile, n; if (fscanf(f, "%d %d", &N, &Dfile) != 2 || fscanf(f, "%d", &n) != 1) { fprintf(stderr, "bad header\n"); return 1; }
    D = npos > 1 ? atoi(pos[1]) : Dfile; if (npos > 2) K = atoi(pos[2]); if (npos > 3) CH = atoi(pos[3]);
    if (N > 62 || D > 9 || K > 12) { fprintf(stderr, "N<=62, D<=9, K<=12\n"); return 1; }
    signal(SIGINT, on_signal); signal(SIGTERM, on_signal);
    for (int i = 0; i < 64; ++i) { binom[i][0] = 1; for (int j = 1; j < 10; ++j) binom[i][j] = i ? binom[i-1][j-1] + binom[i-1][j] : 0; }
    deg_start = malloc((D + 2) * sizeof(int)); deg_start[0] = 0; for (int d = 0; d <= D; ++d) deg_start[d + 1] = deg_start[d] + (int)binom[N][d];
    M = deg_start[D + 1]; WM = (M + 63) / 64;
    mono_mask = malloc((size_t)M * sizeof(u64)); mono_deg = malloc(M);
    for (int d = 0; d <= D; ++d) { u64 m = d ? ((1ull << d) - 1) : 0; if (d == 0) { mono_mask[0] = 0; mono_deg[0] = 0; continue; }
        while (m < (1ull << N)) { int c = col_index(m); mono_mask[c] = m; mono_deg[c] = (uint8_t)d; u64 t = m | (m - 1); m = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzll(m) + 1)); } }
    mulvar = malloc((size_t)N * M * sizeof(int32_t));
    for (int a = 0; a < N; ++a) for (int c = 0; c < M; ++c) { u64 m = mono_mask[c]; mulvar[(size_t)a * M + c] = (mono_deg[c] < D || (m >> a & 1)) ? col_index(m | (1ull << a)) : -1; }
    int nthr = omp_get_max_threads();
    fprintf(stderr, "wdrref: N=%d D=%d monomials=%d threads=%d K=%d tile=%d chunk=%d\n", N, D, M, nthr, K, TILE, CH);

    /* free space: initially all columns, position = column */
    free_pos = malloc(M * sizeof(int)); pos_col = malloc(M * sizeof(int)); is_pivot = calloc(M, 1); piv_of_col = malloc(M * sizeof(long));
    for (int c = 0; c < M; ++c) { free_pos[c] = c; pos_col[c] = c; piv_of_col[c] = -1; }
    NP = M; WF = WM;

    /* batch buffers: full-space product vectors, then (after conversion) free-space residuals */
    u64 *full = malloc((size_t)CH * WM * sizeof(u64));          /* full-space vectors of the batch */
    u64 *res = malloc((size_t)CH * WM * sizeof(u64));           /* free-space residuals (WF <= WM words) */
    int B = 0;
    /* generators into the batch (full space) */
    memset(full, 0, (size_t)CH * WM * sizeof(u64));
    for (int i = 0; i < n; ++i) {
        u64 *r = full + (size_t)B * WM; int ch, any = 0; u64 m;
        while ((ch = fgetc(f)) != EOF) { if (ch == '\n') { if (any) break; else continue; } if (ch == ' ' || ch == '\t' || ch == '\r') continue; ungetc(ch, f);
            if (fscanf(f, "%llx", (unsigned long long *)&m) != 1) break; any = 1; if (__builtin_popcountll(m) <= D) setbit(r, col_index(m)); }
        if (highest_bit(r, WM) >= 0) B++;
    }
    fclose(f);
    /* per-thread tables over the free space */
    u64 *tabs = malloc((size_t)nthr * ((size_t)1 << K) * WM * sizeof(u64));
    int WMR = (M + 63) / 64;                                               /* words for a bitmap over row indices (nrows <= M) */
    u64 *hitbuf = malloc((size_t)nthr * TILE * WMR * sizeof(u64));
    if (!full || !res || !tabs || !hitbuf) { fprintf(stderr, "out of memory for batch buffers\n"); return 1; }
    long *sel = malloc(sizeof(long) * (CH / N + 2)); long *newrows = malloc(sizeof(long) * CH); int *newcol = malloc(sizeof(int) * CH);
    long *stack_rows = malloc(sizeof(long) * (CH + 1));
    long next_expand = 0, total_products = 0; double t0 = now(), last_prog = now(); int done = 0;
    for (;;) {
        /* ===== 1. reduce the batch: full-space vectors -> free-space residuals =====
           A product's bits at pivot columns decide which stored rows are XORed in (no cascading: RREF rows carry no other
           pivot bits), so per product we build a hit bitmap over row indices and a free-space residual, then sweep the stored
           rows in groups of K with a Gray-code table, reading the K-bit pattern straight from the hit bitmap. */
        if (B > 0) {
            long nr_now = nrows; int WR = (int)((nr_now + 63) / 64);
            #pragma omp parallel for schedule(dynamic, 1)
            for (int t = 0; t < (B + TILE - 1) / TILE; ++t) {
                int tid = omp_get_thread_num(); u64 *tab = tabs + (size_t)tid * ((size_t)1 << K) * WM;
                u64 *hits = hitbuf + (size_t)tid * TILE * WMR;
                int i0 = t * TILE, i1 = i0 + TILE < B ? i0 + TILE : B, T = i1 - i0;
                for (int i = 0; i < T; ++i) {
                    const u64 *fv = full + (size_t)(i0 + i) * WM; u64 *rv = res + (size_t)(i0 + i) * WM; u64 *hv = hits + (size_t)i * WMR;
                    memset(rv, 0, (size_t)WF * sizeof(u64)); if (WR) memset(hv, 0, (size_t)WR * sizeof(u64));
                    for (int w = 0; w < WM; ++w) { u64 x = fv[w]; while (x) { int b = __builtin_ctzll(x); x &= x - 1; int c = w * 64 + b;
                        if (is_pivot[c]) { long ri = piv_of_col[c]; hv[ri >> 6] |= 1ull << (ri & 63); } else { int q = free_pos[c]; rv[q >> 6] |= 1ull << (q & 63); } } }
                }
                long g[16];
                for (long b0 = 0; b0 < nr_now; b0 += K) {
                    int kk = (int)(nr_now - b0 < K ? nr_now - b0 : K);
                    /* does any product of the tile hit this group?  (cheap test before building the table) */
                    int any = 0;
                    for (int i = 0; i < T && !any; ++i) { const u64 *hv = hits + (size_t)i * WMR; for (int j = 0; j < kk; ++j) if (hv[(b0 + j) >> 6] >> ((b0 + j) & 63) & 1) { any = 1; break; } }
                    if (!any) continue;
                    for (int j = 0; j < kk; ++j) g[j] = b0 + j;
                    build_table(tab, g, kk, WF);
                    for (int i = 0; i < T; ++i) {
                        const u64 *hv = hits + (size_t)i * WMR; int p = 0;
                        for (int j = 0; j < kk; ++j) p |= (int)((hv[(b0 + j) >> 6] >> ((b0 + j) & 63)) & 1) << j;
                        if (p) { u64 *rv = res + (size_t)(i0 + i) * WM; const u64 *tr = tab + (size_t)p * WF; for (int w = 0; w < WF; ++w) rv[w] ^= tr[w]; }
                    }
                }
            }
            /* ===== 2. self-elimination of residuals (free space), blocked: groups of up to KS new pivots, Gray table,
                        applied to all other residual rows in parallel (per-thread copies not needed: table is read-only) ===== */
            int np = 0;
            {
                const int KS = 8; static u64 *stab = NULL; static size_t stab_cap = 0;
                if (stab_cap < ((size_t)1 << KS) * WM) { stab_cap = ((size_t)1 << KS) * WM; stab = realloc(stab, stab_cap * sizeof(u64)); }
                int i = 0, gidx[16], gq[16];
                for (;;) {
                    int g = 0;
                    for (; i < B && g < KS; ++i) {
                        u64 *r = res + (size_t)i * WM;
                        for (int j = 0; j < g; ++j) if (getbit(r, gq[j])) { const u64 *sr = res + (size_t)gidx[j] * WM; for (int w = 0; w < WF; ++w) r[w] ^= sr[w]; }
                        int h = highest_bit(r, WF); if (h < 0) continue;
                        gidx[g] = i; gq[g] = h; g++;
                    }
                    if (g == 0) break;
                    for (int j = 1; j < g; ++j) for (int e = 0; e < j; ++e) { u64 *re = res + (size_t)gidx[e] * WM; if (getbit(re, gq[j])) { const u64 *sr = res + (size_t)gidx[j] * WM; for (int w = 0; w < WF; ++w) re[w] ^= sr[w]; } }
                    memset(stab, 0, (size_t)WF * sizeof(u64));
                    for (int j = 0; j < g; ++j) { const u64 *rj = res + (size_t)gidx[j] * WM; int half = 1 << j;
                        for (int p = 0; p < half; ++p) { u64 *dst = stab + (size_t)(p | half) * WF; const u64 *src = stab + (size_t)p * WF; memcpy(dst, src, (size_t)WF * sizeof(u64)); for (int w = 0; w < WF; ++w) dst[w] ^= rj[w]; } }
                    #pragma omp parallel for schedule(static)
                    for (int sidx = 0; sidx < B; ++sidx) {
                        int ing = 0; for (int j = 0; j < g; ++j) if (sidx == gidx[j]) { ing = 1; break; }
                        if (ing) continue;
                        u64 *r = res + (size_t)sidx * WM; int p = 0; for (int j = 0; j < g; ++j) p |= getbit(r, gq[j]) << j;
                        if (p) { const u64 *tr = stab + (size_t)p * WF; for (int w = 0; w < WF; ++w) r[w] ^= tr[w]; }
                    }
                    for (int j = 0; j < g; ++j) { newrows[np] = gidx[j]; newcol[np] = pos_col[gq[j]]; np++; }
                    if (i >= B) break;
                }
            }
            /* ===== 3. back-substitute the new pivots into all stored rows (tiled, per-thread tables over the NEW rows) ===== */
            if (np > 0 && nrows > 0) {
                #pragma omp parallel for schedule(dynamic, 1)
                for (long t = 0; t < (nrows + 1023) / 1024; ++t) {
                    int tid = omp_get_thread_num(); u64 *tab = tabs + (size_t)tid * ((size_t)1 << K) * WM;
                    long i0 = t * 1024, i1 = i0 + 1024 < nrows ? i0 + 1024 : nrows;
                    for (int b0 = 0; b0 < np; b0 += K) {
                        int kk = np - b0 < K ? np - b0 : K; int qs[16];
                        for (int j = 0; j < kk; ++j) qs[j] = free_pos[newcol[b0 + j]];
                        memset(tab, 0, (size_t)WF * sizeof(u64));
                        for (int j = 0; j < kk; ++j) { const u64 *rj = res + (size_t)newrows[b0 + j] * WM; int half = 1 << j;
                            for (int p = 0; p < half; ++p) { u64 *dst = tab + (size_t)(p | half) * WF; const u64 *src = tab + (size_t)p * WF; memcpy(dst, src, (size_t)WF * sizeof(u64)); for (int w = 0; w < WF; ++w) dst[w] ^= rj[w]; } }
                        for (long i = i0; i < i1; ++i) { u64 *r = ROWP(i); int p = 0; for (int j = 0; j < kk; ++j) p |= getbit(r, qs[j]) << j;
                            if (p) { const u64 *tr = tab + (size_t)p * WF; for (int w = 0; w < WF; ++w) r[w] ^= tr[w]; } }
                    }
                }
            }
            /* ===== 4. store the new rows; their pivot positions become dead ===== */
            ensure_rows(nrows + np);
            for (int j = 0; j < np; ++j) {
                int c = newcol[j], q = free_pos[c]; u64 *src = res + (size_t)newrows[j] * WM;
                src[q >> 6] &= ~(1ull << (q & 63));                     /* the pivot bit itself is implicit */
                u64 *r = ROWP(nrows); memcpy(r, src, (size_t)WF * sizeof(u64));
                pivcol[nrows] = c; expanded[nrows] = 0; piv_of_col[c] = nrows; is_pivot[c] = 1; nrows++;
                pos_col[q] = -1; free_pos[c] = -1; ndead++; mem_rows += (size_t)WF * 8;
            }
            if (mem_rows > mem_peak) mem_peak = mem_rows;
            if (ndead > compact_frac * NP && ndead > 64) compact();
            B = 0;
        }
        if (nsol >= 0 && nrows >= M - nsol) { fprintf(stderr, "[rref] collapsed: dim V_D = |B_<=D| - %ld; stopping early\n", nsol); done = 1; }
        if (stop_req) { fprintf(stderr, "[rref] interrupted (no checkpointing in wdrref yet)\n"); return 3; }
        double tn = now();
        if (tn - last_prog >= prog_every) { char te[32]; fmt_time(tn - t0, te, sizeof te);
            fprintf(stderr, "[rref] pivots=%ld  expanded %ld/%ld  products=%ld  free cols=%d (dead %ld)  rows=%.2fGB (peak %.2f)  elapsed %s\n", nrows, next_expand, nrows, total_products, NP - (int)ndead, ndead, mem_rows / 1e9, mem_peak / 1e9, te); last_prog = tn; }
        if (done) break;
        /* ===== 5. refill: products x_a * row for unexpanded rows of degree <= D-1 ===== */
        int nsel = 0;
        while ((nsel + 1) * N <= CH && next_expand < nrows) {
            long i = next_expand++;
            if (mono_deg[pivcol[i]] > D - 1) { expanded[i] = 1; continue; }
            expanded[i] = 1; sel[nsel++] = i;
        }
        if (nsel == 0) break;
        memset(full, 0, (size_t)nsel * N * WM * sizeof(u64));
        #pragma omp parallel for schedule(dynamic, 4)
        for (long k = 0; k < (long)nsel * N; ++k) {
            long i = sel[k / N]; int a = (int)(k % N); u64 *fv = full + (size_t)k * WM;
            const int32_t *mv = mulvar + (size_t)a * M;
            int c0 = mv[pivcol[i]]; if (c0 >= 0) setbit(fv, c0);            /* x_a * pivot monomial */
            const u64 *r = ROWP(i);
            for (int w = 0; w < WF; ++w) { u64 x = r[w]; while (x) { int b = __builtin_ctzll(x); x &= x - 1; int c = pos_col[w * 64 + b]; if (c >= 0) { int cc = mv[c]; if (cc >= 0) setbit(fv, cc); } } }
        }
        /* compact non-zero products to the front */
        for (long k = 0; k < (long)nsel * N; ++k) { u64 *fv = full + (size_t)k * WM; if (highest_bit(fv, WM) >= 0) { if (k != B) memcpy(full + (size_t)B * WM, fv, (size_t)WM * sizeof(u64)); B++; total_products++; } }
    }
    /* ===== report ===== */
    long *cnt = calloc(D + 1, sizeof(long)); for (long i = 0; i < nrows; ++i) cnt[mono_deg[pivcol[i]]]++;
    printf("N=%d D=%d dim(V_D)=%ld  products=%ld  time=%.0fs  rows_mem=%.2fGB peak=%.2fGB  (rref)\n", N, D, nrows, total_products, now() - t0, mem_rows / 1e9, mem_peak / 1e9);
    long cum = 0; for (int d = 0; d <= D; ++d) { cum += cnt[d]; printf("  dim(V_%d cap B_<=%d) = %ld\n", D, d, cum); }
    long *npv = calloc(D + 1, sizeof(long)); long nptot = 0; for (int c = 0; c < M; ++c) if (!is_pivot[c]) { npv[mono_deg[c]]++; nptot++; }
    printf("  non-pivot monomials: total %ld  by degree:", nptot); for (int d = 0; d <= D; ++d) printf(" %d:%ld", d, npv[d]); printf("\n");
    return 0;
}
