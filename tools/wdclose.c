/*
 * wdclose.c -- degree-capped closure V_D (HKY Def. 1) of a set of Boolean
 * polynomials over GF(2), i.e. the smallest space containing the generators
 * (degree <= D) and closed under g -> x_a*g whenever deg g <= D-1.
 *
 * Output: dim(V_D cap B_{<=k}) for k = 1..D.
 *
 * Build:  gcc -O3 -march=native -o wdclose wdclose.c
 * Usage:  ./wdclose input.txt [D] [k_table] [chunk]
 *
 * Input format (written by wd_export_c.py):
 *   N D
 *   n
 *   <n lines: space-separated monomials as hex bitmasks over the N variables; "0" = constant 1>
 *
 * Linear algebra: rows are bit vectors over the squarefree monomials of degree
 * <= D (degree-major order, colex inside a degree).  The echelon is kept
 * block-triangular: every pivot row has zeros at the pivot columns of all
 * earlier chunks and the rows of one chunk are in RREF among themselves, so a
 * new row is reduced by sweeping the blocks in insertion order with Gray-code
 * tables of k rows (Method of Four Russians).  Pivot = highest set bit, hence
 * the pivot column also gives the row's degree.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
static double now(void){ struct timespec ts; clock_gettime(CLOCK_REALTIME,&ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }
#include <signal.h>
#include <unistd.h>
static volatile sig_atomic_t stop_req = 0;
static void on_signal(int sig) { (void)sig; if (stop_req) _exit(130); stop_req = 1; }
static void fmt_time(double s, char *buf, size_t n) { long t = (long)(s + 0.5); snprintf(buf, n, "%ld:%02ld:%02ld", t / 3600, (t / 60) % 60, t % 60); }

typedef uint64_t u64;
static int N, D, W, M;
static int *deg_start;           /* first column index of each degree */
static u64 *mono_mask;           /* column -> monomial bitmask */
static uint8_t *mono_deg;
static u64 binom[64][10];
static int32_t *mulvar;          /* [a*M + idx] -> column of x_a * monomial(idx) (deg < D or contains a) */

static int col_index(u64 mask) {   /* degree-major, colex within degree */
    int d = __builtin_popcountll(mask), j = 1; u64 r = 0; u64 m = mask;
    while (m) { int v = __builtin_ctzll(m); r += binom[v][j]; j++; m &= m - 1; }
    return deg_start[d] + (int)r;
}

/* ---------------- echelon storage ---------------- */
#define RB 4096                   /* rows per storage block */
static u64 **rowptr = NULL;       /* row i -> pointer to its (variable-length) words */
static int *rowlen = NULL;        /* number of words stored for row i = (pivcol>>6)+1 */
static long cap_rows = 0, nrows = 0; static size_t mem_rows = 0;
static int *pivcol = NULL;        /* pivot column of each stored row */
static char *expanded = NULL;     /* row already multiplied by all variables */
static int *chunk_end = NULL;     /* chunk boundaries: rows [chunk_start, chunk_end) share a chunk */
static int nchunks = 0, cap_chunks = 0;
static char *is_pivot_col = NULL;
static char ckpt_path[4096] = ""; static double ckpt_every = 600, prog_every = 30;
static u64 CKPT_MAGIC = 0x57444B5054000007ull;   /* "WDKPT" v7: per-degree expansion queues */
/* expansion order: lowest degree first. deg_rows[d] = indices of stored rows of degree d (insertion order); qpos[d] = next to expand */
static long **deg_rows = NULL, *deg_len = NULL, *deg_cap = NULL, *qpos = NULL, bucketed = 0;
static int lowfirst = 0;            /* -lowfirst: expand rows of lowest degree first; default: insertion order */
static long next_expand = 0;        /* insertion-order cursor */
static void bucket_rows(void) {
    for (; bucketed < nrows; ++bucketed) {
        int d = mono_deg[pivcol[bucketed]];
        if (deg_len[d] >= deg_cap[d]) { deg_cap[d] = deg_cap[d] ? 2 * deg_cap[d] : 1024; deg_rows[d] = realloc(deg_rows[d], deg_cap[d] * sizeof(long)); }
        deg_rows[d][deg_len[d]++] = bucketed;
    }
}
static long expanded_count(void) { if (lowfirst) { long s = 0; for (int d = 0; d <= D; ++d) s += qpos[d]; return s; } long s = 0; for (long i = 0; i < nrows; ++i) s += expanded[i] ? 1 : 0; return s; }

static void ensure_rows(long need) {
    if (need <= cap_rows) return;
    long nc = cap_rows + RB;
    while (nc < need) nc += RB;
    rowptr = realloc(rowptr, nc * sizeof(u64 *));
    rowlen = realloc(rowlen, nc * sizeof(int));
    pivcol = realloc(pivcol, nc * sizeof(int));
    expanded = realloc(expanded, nc);
    if (!rowptr || !rowlen || !pivcol || !expanded) { fprintf(stderr, "out of memory (%ld row headers)\n", nc); exit(1); }
    cap_rows = nc;
}
static inline u64 *ROW(long i) { return rowptr[i]; }
static void store_row(long i, const u64 *r, int h) {
    int len = (h >> 6) + 1;
    u64 *p = malloc((size_t)len * sizeof(u64));
    if (!p) { fprintf(stderr, "out of memory storing row %ld (%.2f GB in rows)\n", i, mem_rows / 1e9); exit(1); }
    memcpy(p, r, (size_t)len * sizeof(u64));
    rowptr[i] = p; rowlen[i] = len; mem_rows += (size_t)len * sizeof(u64);
}
static inline void xorrow_len(u64 *restrict a, const u64 *restrict b, int len) { for (int w = 0; w < len; ++w) a[w] ^= b[w]; }
static inline int highest_bit(const u64 *r) {
    for (int w = W - 1; w >= 0; --w) if (r[w]) return w * 64 + 63 - __builtin_clzll(r[w]);
    return -1;
}
static inline int getbit(const u64 *r, int c) { return (r[c >> 6] >> (c & 63)) & 1; }
static inline void xorrow(u64 *restrict a, const u64 *restrict b) { for (int w = 0; w < W; ++w) a[w] ^= b[w]; }

/* ---------------- Four-Russians reduction of a batch against the echelon ---------------- */
static int K = 10;
static u64 *table = NULL;   /* (1<<K) * W */
static int TILE = 0;        /* -tile T: reduce the batch in tiles of T rows, each thread owns a tile and its own table (no barriers) */
static u64 *ttables = NULL; static int nthr_alloc = 0;

static void reduce_tile(u64 *rows, int T, u64 *tab) {
    int cs = 0;
    for (int c = 0; c < nchunks; ++c) {
        int ce = chunk_end[c];
        for (int b0 = cs; b0 < ce; b0 += K) {
            int kk = ce - b0 < K ? ce - b0 : K;
            int cols[16]; int L = 0;
            for (int j = 0; j < kk; ++j) { cols[j] = pivcol[b0 + j]; if (rowlen[b0 + j] > L) L = rowlen[b0 + j]; }
            memset(tab, 0, (size_t)L * sizeof(u64));
            for (int j = 0; j < kk; ++j) {
                const u64 *rj = ROW(b0 + j); const int lj = rowlen[b0 + j]; const int half = 1 << j;
                for (int p = 0; p < half; ++p) {
                    u64 *dst = tab + (size_t)(p | half) * W; const u64 *src = tab + (size_t)p * W;
                    memcpy(dst, src, (size_t)L * sizeof(u64));
                    for (int w = 0; w < lj; ++w) dst[w] ^= rj[w];
                }
            }
            for (int i = 0; i < T; ++i) {
                u64 *r = rows + (size_t)i * W;
                int p = 0;
                for (int j = 0; j < kk; ++j) p |= getbit(r, cols[j]) << j;
                if (p) { const u64 *t = tab + (size_t)p * W; for (int w = 0; w < L; ++w) r[w] ^= t[w]; }
            }
        }
        cs = ce;
    }
}

static void reduce_batch(u64 *batch, int B) {
    if (TILE > 0 && TILE < B) {
        /* tiled mode: each thread reduces whole tiles independently with its own table -- no synchronisation at all;
           the tile rows stay in the core's cache while the pivot rows are streamed once per tile */
        int ntiles = (B + TILE - 1) / TILE;
        #pragma omp parallel for schedule(dynamic, 1)
        for (int t = 0; t < ntiles; ++t) {
            int tid = omp_get_thread_num();
            int T = (t + 1) * TILE <= B ? TILE : B - t * TILE;
            reduce_tile(batch + (size_t)t * TILE * W, T, ttables + (size_t)tid * ((size_t)1 << K) * W);
        }
        return;
    }
    /* one parallel region for the whole batch; worksharing loops with barriers inside (fork/join per chunk is far too
       expensive on some OpenMP runtimes, e.g. libgomp under Cygwin). All threads compute the per-group bookkeeping
       (cols, L) redundantly from read-only data, so nothing needs to be shared explicitly. */
    #pragma omp parallel
    {
        int cs = 0;
        for (int c = 0; c < nchunks; ++c) {
            int ce = chunk_end[c];
            for (int b0 = cs; b0 < ce; b0 += K) {
                int kk = ce - b0 < K ? ce - b0 : K;
                int cols[16]; int L = 0;
                for (int j = 0; j < kk; ++j) { cols[j] = pivcol[b0 + j]; if (rowlen[b0 + j] > L) L = rowlen[b0 + j]; }
                #pragma omp for schedule(static)
                for (int w = 0; w < L; ++w) table[w] = 0;
                for (int j = 0; j < kk; ++j) {                       /* Gray-code table by doubling; barrier after each step */
                    const u64 *rj = ROW(b0 + j); const int lj = rowlen[b0 + j]; const int half = 1 << j;
                    #pragma omp for schedule(static)
                    for (int p = 0; p < half; ++p) {
                        u64 *dst = table + (size_t)(p | half) * W; const u64 *src = table + (size_t)p * W;
                        memcpy(dst, src, (size_t)L * sizeof(u64));
                        for (int w = 0; w < lj; ++w) dst[w] ^= rj[w];
                    }
                }
                #pragma omp for schedule(static)
                for (int i = 0; i < B; ++i) {
                    u64 *r = batch + (size_t)i * W;
                    int p = 0;
                    for (int j = 0; j < kk; ++j) p |= getbit(r, cols[j]) << j;
                    if (p) { const u64 *t = table + (size_t)p * W; for (int w = 0; w < L; ++w) r[w] ^= t[w]; }
                }
            }
            cs = ce;
        }
    }
}

/* self-elimination inside a batch: returns number of new pivot rows, moved to the front, RREF among themselves */
static int KS = 8;                 /* group size for the blocked in-batch elimination */
static u64 *stable = NULL;         /* (1<<KS) * W */
static int self_eliminate(u64 *batch, int B, int *piv) {
    /* Blocked (Four-Russians) elimination inside the batch. Result: pivot rows listed in piv[0..np), in reduced echelon
       form among themselves (zeros at each other's pivot columns); all other rows have zeros at the pivot columns. */
    int np = 0, i = 0;
    int gidx[16], gcol[16];
    for (;;) {
        /* 1) collect the next group of up to KS pivots, reducing each candidate by the earlier members of the group */
        int g = 0;
        for (; i < B && g < KS; ++i) {
            u64 *r = batch + (size_t)i * W;
            for (int j = 0; j < g; ++j) if (getbit(r, gcol[j])) xorrow(r, batch + (size_t)gidx[j] * W);
            int h = highest_bit(r);
            if (h < 0) continue;
            gidx[g] = i; gcol[g] = h; g++; piv[np++] = i;
        }
        if (g == 0) break;
        /* 2) reduced echelon form inside the group: clear later pivot columns from earlier rows */
        for (int j = 1; j < g; ++j) for (int e = 0; e < j; ++e) { u64 *re = batch + (size_t)gidx[e] * W; if (getbit(re, gcol[j])) xorrow(re, batch + (size_t)gidx[j] * W); }
        /* 3) Gray-code table over the group and application to every other row of the batch */
        int L = 0; for (int j = 0; j < g; ++j) { int l = highest_bit(batch + (size_t)gidx[j] * W) / 64 + 1; if (l > L) L = l; }
        #pragma omp parallel
        {
            #pragma omp for schedule(static)
            for (int w = 0; w < L; ++w) stable[w] = 0;
            for (int j = 0; j < g; ++j) {
                const u64 *rj = batch + (size_t)gidx[j] * W; const int half = 1 << j;
                #pragma omp for schedule(static)
                for (int p = 0; p < half; ++p) {
                    u64 *dst = stable + (size_t)(p | half) * W; const u64 *src = stable + (size_t)p * W;
                    memcpy(dst, src, (size_t)L * sizeof(u64));
                    for (int w = 0; w < L; ++w) dst[w] ^= rj[w];
                }
            }
            #pragma omp for schedule(static)
            for (int s = 0; s < B; ++s) {
                int ing = 0; for (int j = 0; j < g; ++j) if (s == gidx[j]) { ing = 1; break; }
                if (ing) continue;
                u64 *r = batch + (size_t)s * W;
                int p = 0; for (int j = 0; j < g; ++j) p |= getbit(r, gcol[j]) << j;
                if (p) { const u64 *t = stable + (size_t)p * W; for (int w = 0; w < L; ++w) r[w] ^= t[w]; }
            }
        }
        if (i >= B) break;
    }
    return np;
}

/* ---------- checkpoints: the whole echelon + loop position, written atomically (tmp + rename) ---------- */
static int save_ckpt(int phase, long next_expand, long total_products, double elapsed, long nW, u64 **Wrows, int *Wlen, long wpos, int cpos) {
    if (!ckpt_path[0]) return 0;
    char tmp[4200]; snprintf(tmp, sizeof tmp, "%s.tmp", ckpt_path);
    FILE *g = fopen(tmp, "wb"); if (!g) { perror("checkpoint open"); return -1; }
    int hdr[6] = { N, D, M, phase, nchunks, cpos }; long lng[6] = { nrows, next_expand, total_products, nW, wpos, lowfirst };
    int ok = fwrite(&CKPT_MAGIC, 8, 1, g) == 1 && fwrite(hdr, sizeof hdr, 1, g) == 1 && fwrite(lng, sizeof lng, 1, g) == 1 && fwrite(&elapsed, 8, 1, g) == 1;
    ok = ok && fwrite(qpos, sizeof(long), D + 1, g) == (size_t)(D + 1);
    ok = ok && (nrows == 0 || fwrite(expanded, 1, nrows, g) == (size_t)nrows);
    ok = ok && (nrows == 0 || (fwrite(pivcol, sizeof(int), nrows, g) == (size_t)nrows && fwrite(rowlen, sizeof(int), nrows, g) == (size_t)nrows));
    ok = ok && (nchunks == 0 || fwrite(chunk_end, sizeof(int), nchunks, g) == (size_t)nchunks);
    for (long i = 0; ok && i < nrows; ++i) ok = fwrite(rowptr[i], 8, rowlen[i], g) == (size_t)rowlen[i];
    if (phase == 2) { ok = ok && fwrite(Wlen, sizeof(int), nW, g) == (size_t)nW; for (long w = 0; ok && w < nW; ++w) ok = fwrite(Wrows[w], 8, Wlen[w], g) == (size_t)Wlen[w]; }
    if (fclose(g) != 0) ok = 0;
    if (!ok) { fprintf(stderr, "checkpoint write failed (disk full?)\n"); remove(tmp); return -1; }
    if (rename(tmp, ckpt_path) != 0) { perror("checkpoint rename"); return -1; }
    return 0;
}
/* returns phase (1 or 2) on success, 0 if no checkpoint, -1 on mismatch */
static int load_ckpt(long *next_expand, long *total_products, double *elapsed, long *nW, u64 ***Wrows, int **Wlen, long *wpos, int *cpos) {
    if (!ckpt_path[0]) return 0;
    FILE *g = fopen(ckpt_path, "rb"); if (!g) return 0;
    u64 magic; int hdr[6]; long lng[6];
    if (fread(&magic, 8, 1, g) != 1 || (magic != CKPT_MAGIC && magic != 0x57444B5054000006ull)) { fclose(g); fprintf(stderr, "checkpoint %s is from another version of wdclose (use -fresh to start over)\n", ckpt_path); return -1; }
    int oldfmt = (magic == 0x57444B5054000006ull);      /* v6/v7 checkpoint: insertion-order cursor only */
    if (fread(hdr, sizeof hdr, 1, g) != 1 || fread(lng, sizeof lng, 1, g) != 1 || fread(elapsed, 8, 1, g) != 1) { fclose(g); return -1; }
    long nr0 = lng[0]; ensure_rows(nr0);
    if (!oldfmt) {
        if (fread(qpos, sizeof(long), D + 1, g) != (size_t)(D + 1)) { fclose(g); return -1; }
        if (nr0 > 0 && fread(expanded, 1, nr0, g) != (size_t)nr0) { fclose(g); return -1; }
        lowfirst = (int)lng[5];
    } else { memset(qpos, 0, sizeof(long) * (D + 1)); lowfirst = 0; }
    *next_expand = lng[1];
    if (hdr[0] != N || hdr[1] != D || hdr[2] != M) { fclose(g); fprintf(stderr, "checkpoint %s is for N=%d D=%d, not for this run\n", ckpt_path, hdr[0], hdr[1]); return -1; }
    int phase = hdr[3]; nchunks = hdr[4]; *cpos = hdr[5];
    long nr = lng[0]; *next_expand = lng[1]; *total_products = lng[2]; *nW = lng[3]; *wpos = lng[4];
    ensure_rows(nr);
    int *pc = malloc(sizeof(int) * (nr + 1)), *rl = malloc(sizeof(int) * (nr + 1));
    if (nr > 0 && (fread(pc, sizeof(int), nr, g) != (size_t)nr || fread(rl, sizeof(int), nr, g) != (size_t)nr)) { fclose(g); return -1; }
    cap_chunks = nchunks + 1024; chunk_end = realloc(chunk_end, cap_chunks * sizeof(int));
    if (nchunks > 0 && fread(chunk_end, sizeof(int), nchunks, g) != (size_t)nchunks) { fclose(g); return -1; }
    for (long i = 0; i < nr; ++i) {
        u64 *p = malloc((size_t)rl[i] * 8); if (!p || fread(p, 8, rl[i], g) != (size_t)rl[i]) { fclose(g); return -1; }
        rowptr[i] = p; rowlen[i] = rl[i]; pivcol[i] = pc[i]; mem_rows += (size_t)rl[i] * 8; is_pivot_col[pc[i]] = 1;
        if (oldfmt) expanded[i] = (i < *next_expand && mono_deg[pc[i]] <= ((phase == 2) ? -1 : D - 1)) ? 1 : 0;
    }
    nrows = nr; free(pc); free(rl);
    bucket_rows();
    if (phase == 2) {
        *Wlen = malloc(sizeof(int) * (*nW + 1)); *Wrows = malloc(sizeof(u64 *) * (*nW + 1));
        if (fread(*Wlen, sizeof(int), *nW, g) != (size_t)*nW) { fclose(g); return -1; }
        for (long w = 0; w < *nW; ++w) { (*Wrows)[w] = malloc((size_t)(*Wlen)[w] * 8); if (fread((*Wrows)[w], 8, (*Wlen)[w], g) != (size_t)(*Wlen)[w]) { fclose(g); return -1; } }
    }
    fclose(g); return phase;
}

int main(int argc, char **argv) {
    const char *pos[8] = {0}; int npos = 0, fresh = 0; char ckpt_arg[4096] = ""; long nsol = -1;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-every") && i + 1 < argc) ckpt_every = atof(argv[++i]);
        else if (!strcmp(argv[i], "-prog") && i + 1 < argc) prog_every = atof(argv[++i]);
        else if (!strcmp(argv[i], "-ckpt") && i + 1 < argc) snprintf(ckpt_arg, sizeof ckpt_arg, "%s", argv[++i]);
        else if (!strcmp(argv[i], "-fresh")) fresh = 1;
        else if (!strcmp(argv[i], "-lowfirst")) lowfirst = 1;
        else if (!strcmp(argv[i], "-tile") && i + 1 < argc) TILE = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-nsol") && i + 1 < argc) nsol = atol(argv[++i]);   /* early exit when V_D = J_<=D (all but #sol monomials are pivots) */
        else if (npos < 8) pos[npos++] = argv[i];
    }
    if (npos < 1) {
        fprintf(stderr, "usage: %s input.txt [D] [k] [chunk] [macaulay] [-every SEC] [-prog SEC] [-ckpt FILE] [-fresh] [-nsol S] [-lowfirst] [-tile T]\n"
                        "  progress line every -prog seconds (default 30); checkpoint every -every seconds (default 600, 0 = only on Ctrl-C)\n"
                        "  Ctrl-C (once) writes a checkpoint and exits with code 3; rerunning the same command resumes from it;\n"
                        "  the checkpoint file is <input>.D<D>[.mac].ckpt unless -ckpt is given, and is removed when the run completes; -fresh ignores it\n", argv[0]);
        return 1;
    }
    FILE *f = fopen(pos[0], "r"); if (!f) { perror("input"); return 1; }
    int Dfile, n;
    if (fscanf(f, "%d %d", &N, &Dfile) != 2 || fscanf(f, "%d", &n) != 1) { fprintf(stderr, "bad header\n"); return 1; }
    D = npos > 1 ? atoi(pos[1]) : Dfile;
    if (npos > 2) K = atoi(pos[2]);
    int CH = npos > 3 ? atoi(pos[3]) : 4096;
    int MACAULAY = (npos > 4 && strcmp(pos[4], "macaulay") == 0);   /* rank of {t*q : q in W_2, deg t <= D-2} */
    if (ckpt_arg[0]) snprintf(ckpt_path, sizeof ckpt_path, "%s", ckpt_arg); else snprintf(ckpt_path, sizeof ckpt_path, "%s.D%d%s.ckpt", pos[0], D, MACAULAY ? ".mac" : "");
    if (fresh) remove(ckpt_path);
    signal(SIGINT, on_signal); signal(SIGTERM, on_signal);
    if (N > 62 || D > 9 || K > 14) { fprintf(stderr, "N<=62, D<=9, k<=14\n"); return 1; }
    for (int i = 0; i < 64; ++i) { binom[i][0] = 1; for (int j = 1; j < 10; ++j) binom[i][j] = i ? binom[i-1][j-1] + binom[i-1][j] : 0; }
    deg_start = calloc(D + 2, sizeof(int));
    M = 0;
    for (int d = 0; d <= D; ++d) { deg_start[d] = M; M += (int)binom[N][d]; }
    deg_start[D + 1] = M;
    W = (M + 63) / 64;
    mono_mask = malloc((size_t)M * sizeof(u64)); mono_deg = malloc(M);
    for (int d = 0; d <= D; ++d) {                       /* enumerate d-subsets of N variables (Gosper's hack) */
        if (d == 0) { mono_mask[0] = 0; mono_deg[0] = 0; continue; }
        u64 m = (1ull << d) - 1, last = m << (N - d);
        for (;;) {
            int c = col_index(m); mono_mask[c] = m; mono_deg[c] = (uint8_t)d;
            if (m == last) break;
            u64 lo = m & -m, ripple = m + lo, ones = ((m ^ ripple) >> 2) / lo;   /* next subset in colex order */
            m = ripple | ones;
        }
    }
    mulvar = malloc((size_t)N * M * sizeof(int32_t));
    for (int a = 0; a < N; ++a) for (int c = 0; c < M; ++c) {
        u64 m = mono_mask[c];
        mulvar[(size_t)a * M + c] = (mono_deg[c] < D || (m >> a & 1)) ? col_index(m | (1ull << a)) : -1;
    }
    is_pivot_col = calloc(M, 1);
    table = malloc(((size_t)1 << K) * W * sizeof(u64));
    if (TILE > 0) { nthr_alloc = omp_get_max_threads(); ttables = malloc((size_t)nthr_alloc * ((size_t)1 << K) * W * sizeof(u64)); if (!ttables) { fprintf(stderr, "cannot allocate %d per-thread tables\n", nthr_alloc); return 1; } }
    stable = malloc(((size_t)1 << KS) * W * sizeof(u64));
    u64 *batch = malloc((size_t)CH * W * sizeof(u64));
    int *piv = malloc(CH * sizeof(int));
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    fprintf(stderr, "N=%d D=%d monomials=%d words/row=%d  threads=%d  (rows stored to their pivot; worst case %.2f GB)\n", N, D, M, W, nthreads, (double)M * W * 4 / 1e9);

    deg_rows = calloc(D + 1, sizeof(long *)); deg_len = calloc(D + 1, sizeof(long)); deg_cap = calloc(D + 1, sizeof(long)); qpos = calloc(D + 1, sizeof(long));
    /* resume from a checkpoint if one exists for this input and D */
    long total_products = 0, nW = 0, wpos = 0; int cpos = 0; double elapsed0 = 0;
    u64 **Wrows = NULL; int *Wlen = NULL;
    int phase = load_ckpt(&next_expand, &total_products, &elapsed0, &nW, &Wrows, &Wlen, &wpos, &cpos);
    if (phase < 0) { fprintf(stderr, "cannot use checkpoint %s (remove it or use -fresh)\n", ckpt_path); return 1; }
    if (phase > 0) { char tb[32]; fmt_time(elapsed0, tb, sizeof tb); fprintf(stderr, "resumed from %s: phase %d, pivots=%ld, expanded %ld rows, products=%ld, %s elapsed before\n", ckpt_path, phase, nrows, next_expand, total_products, tb); }
    /* read generators into the first batch */
    int B = 0;
    memset(batch, 0, (size_t)CH * W * sizeof(u64));
    for (int i = 0; i < n && phase == 0; ++i) {
        u64 *r = batch + (size_t)B * W;
        int ch; u64 m;
        /* read tokens until end of line */
        int any = 0;
        while ((ch = fgetc(f)) != EOF) {
            if (ch == '\n') { if (any) break; else continue; }
            if (ch == ' ' || ch == '\t' || ch == '\r') continue;
            ungetc(ch, f);
            if (fscanf(f, "%llx", (unsigned long long *)&m) != 1) break;
            any = 1;
            if (__builtin_popcountll(m) <= D) { int c = col_index(m); r[c >> 6] ^= 1ull << (c & 63); }
        }
        if (highest_bit(r) >= 0) B++;
    }
    fclose(f);

    double t0 = now() - elapsed0, last_ckpt = now(), last_prog = now();
    long iter = 0;                        /* next_expand = first stored row not yet expanded */
    if (phase != 2) for (;;) {
        /* --- process current batch --- */
        if (B > 0) {
            reduce_batch(batch, B);
            int np = self_eliminate(batch, B, piv);
            if (np > 0) {
                ensure_rows(nrows + np);
                if (nchunks >= cap_chunks) { cap_chunks = cap_chunks ? cap_chunks * 2 : 1024; chunk_end = realloc(chunk_end, cap_chunks * sizeof(int)); }
                for (int j = 0; j < np; ++j) {
                    u64 *r = batch + (size_t)piv[j] * W;
                    int h = highest_bit(r);
                    store_row(nrows, r, h);
                    pivcol[nrows] = h; is_pivot_col[h] = 1; expanded[nrows] = 0;
                    nrows++;
                }
                chunk_end[nchunks++] = (int)nrows;
            }
            B = 0;
        }
        {   /* --- progress line and checkpoint (batch is empty here, all rows < next_expand fully processed) --- */
            double tn = now();
            if (tn - last_prog >= prog_every) {
                char te[32], tr[32]; long ex = expanded_count(); double frac = nrows ? (double)ex / nrows : 0; fmt_time(tn - t0, te, sizeof te);
                if (frac > 0.02) fmt_time((tn - t0) * (1 - frac) / frac, tr, sizeof tr); else snprintf(tr, sizeof tr, "?");
                fprintf(stderr, "[closure] pivots=%ld  expanded %ld/%ld rows (%.1f%%)  products=%ld  rows=%.2fGB  elapsed %s  remaining ~%s (rough; rows still grow)\n",
                        nrows, ex, nrows, 100 * frac, total_products, mem_rows / 1e9, te, tr);
                last_prog = tn;
            }
            if (stop_req || (ckpt_every > 0 && tn - last_ckpt >= ckpt_every)) {
                if (save_ckpt(1, next_expand, total_products, tn - t0, 0, NULL, NULL, 0, 0) == 0) fprintf(stderr, "[checkpoint] %s written (%.2f GB)%s\n", ckpt_path, mem_rows / 1e9, stop_req ? " -- interrupted; rerun the same command to resume" : "");
                last_ckpt = now();
                if (stop_req) return 3;
            }
        }
        /* --- refill batch with products x_a * (unexpanded rows of degree <= D-1) --- */
        {
            /* 1) pick the rows to expand (sequential bookkeeping), 2) generate all their products in ONE parallel region, 3) compact */
            static long *sel = NULL; static int selcap = 0;
            int nsel = 0, cap = MACAULAY ? 2 : D - 1;
            if (selcap < CH / N + 2) { selcap = CH / N + 2; sel = realloc(sel, selcap * sizeof(long)); }
            while ((nsel + 1) * N <= CH - B) {
                bucket_rows();
                long i = -1;
                if (lowfirst) {
                    int dsel = -1;
                    for (int d = 0; d <= cap; ++d) if (qpos[d] < deg_len[d]) { dsel = d; break; }
                    if (dsel < 0) break;
                    i = deg_rows[dsel][qpos[dsel]++];
                } else {
                    while (next_expand < nrows && (expanded[next_expand] || mono_deg[pivcol[next_expand]] > cap)) { if (mono_deg[pivcol[next_expand]] > cap) expanded[next_expand] = 1; next_expand++; }
                    if (next_expand >= nrows) break;
                    i = next_expand++;
                }
                expanded[i] = 1; sel[nsel++] = i;
            }
            if (nsel > 0) {
                u64 *base = batch + (size_t)B * W;          /* nsel*N consecutive slots are free */
                #pragma omp parallel for schedule(dynamic, 4)
                for (long k = 0; k < (long)nsel * N; ++k) {
                    long i = sel[k / N]; int a = (int)(k % N);
                    const u64 *src = ROW(i); const int slen = rowlen[i];
                    u64 *r = base + (size_t)k * W;
                    memset(r, 0, (size_t)W * sizeof(u64));
                    const int32_t *mv = mulvar + (size_t)a * M;
                    for (int w = 0; w < slen; ++w) {
                        u64 x = src[w];
                        while (x) { int b = __builtin_ctzll(x); x &= x - 1; int c = mv[w * 64 + b]; r[c >> 6] ^= 1ull << (c & 63); }
                    }
                }
                for (long k = 0; k < (long)nsel * N; ++k) {   /* compact: drop zero products */
                    u64 *r = base + (size_t)k * W;
                    if (highest_bit(r) >= 0) {
                        if (r != batch + (size_t)B * W) memcpy(batch + (size_t)B * W, r, (size_t)W * sizeof(u64));
                        B++; total_products++;
                    }
                }
            }
        }
        if (B == 0) {
            bucket_rows(); int more = 0, cap = MACAULAY ? 2 : D - 1;
            if (lowfirst) { for (int d = 0; d <= cap; ++d) if (qpos[d] < deg_len[d]) more = 1; }
            else { for (long r = next_expand; r < nrows; ++r) if (!expanded[r] && mono_deg[pivcol[r]] <= cap) { more = 1; break; } }
            if (!more) break;
        }
        if (nsol >= 0 && nrows >= M - nsol) { fprintf(stderr, "[closure] collapsed: dim V_D = |B_<=D| - %ld, nothing more can be added; stopping early\n", nsol); B = 0; break; }
        if (++iter % 16 == 0)
            fprintf(stderr, "  chunk %ld: pivots=%ld products=%ld  rows=%.2fGB  %.0fs\n", iter, nrows, total_products, mem_rows / 1e9, (now() - t0));
    }
    if (MACAULAY) {
        /* Stage 1 only: take W_2 = all pivot rows of degree <= 2 of the (full) closure just computed,
           reset the echelon, and insert every product t*w with deg t <= D-2 -- no expansion of falls. */
        long nprod = total_products; B = 0;
        if (phase != 2) {
            nW = 0; Wrows = malloc(sizeof(u64*) * (nrows + 1)); Wlen = malloc(sizeof(int) * (nrows + 1));
            for (long i = 0; i < nrows; ++i) if (mono_deg[pivcol[i]] <= 2) { Wrows[nW] = rowptr[i]; Wlen[nW] = rowlen[i]; nW++; }
            fprintf(stderr, "macaulay mode: |W_2| = %ld rows of degree <= 2; multipliers of degree <= %d\n", nW, D - 2);
            /* reset echelon (keep W rows alive: they are separate mallocs referenced by Wrows) */
            nrows = 0; nchunks = 0; memset(is_pivot_col, 0, M); mem_rows = 0;
            nprod = 0; wpos = 0; cpos = 0;
        }
        long ncols_t = deg_start[D - 1];
        for (long w = wpos; w < nW; ++w) {
            for (int c = (w == wpos ? cpos : 0); c < ncols_t; ++c) {          /* all monomials t of degree <= D-2 */
                u64 t = mono_mask[c];
                u64 *r = batch + (size_t)B * W;
                memset(r, 0, (size_t)W * sizeof(u64));
                for (int k = 0; k < Wlen[w]; ++k) r[k] = Wrows[w][k];
                u64 tt = t;
                while (tt) {                                       /* multiply by the variables of t, one at a time */
                    int a = __builtin_ctzll(tt); tt &= tt - 1;
                    const int32_t *mv = mulvar + (size_t)a * M;
                    u64 *r2 = batch + (size_t)(B + 1) * W;         /* scratch slot (CH >= N+2 guaranteed below) */
                    memset(r2, 0, (size_t)W * sizeof(u64));
                    for (int k = 0; k < W; ++k) { u64 x = r[k]; while (x) { int b = __builtin_ctzll(x); x &= x - 1; int cc = mv[k * 64 + b]; r2[cc >> 6] ^= 1ull << (cc & 63); } }
                    memcpy(r, r2, (size_t)W * sizeof(u64));
                }
                if (highest_bit(r) >= 0) { B++; nprod++; }
                if (B >= CH - 2) {
                    reduce_batch(batch, B);
                    int np = self_eliminate(batch, B, piv);
                    if (np > 0) {
                        ensure_rows(nrows + np);
                        if (nchunks >= cap_chunks) { cap_chunks = cap_chunks ? cap_chunks * 2 : 1024; chunk_end = realloc(chunk_end, cap_chunks * sizeof(int)); }
                        for (int j = 0; j < np; ++j) { u64 *q = batch + (size_t)piv[j] * W; int h = highest_bit(q); store_row(nrows, q, h); pivcol[nrows] = h; is_pivot_col[h] = 1; nrows++; }
                        chunk_end[nchunks++] = (int)nrows;
                    }
                    B = 0;
                    double tn = now();
                    if (tn - last_prog >= prog_every) {
                        char te[32], tr[32]; double frac = ((double)w * ncols_t + c + 1) / ((double)nW * ncols_t); fmt_time(tn - t0, te, sizeof te);
                        if (frac > 0.01) fmt_time((tn - t0) * (1 - frac) / frac, tr, sizeof tr); else snprintf(tr, sizeof tr, "?");
                        fprintf(stderr, "[macaulay] pivots=%ld  products %ld  done %.1f%% of W_2 x monomials  rows=%.2fGB  elapsed %s  remaining ~%s\n", nrows, nprod, 100 * frac, mem_rows / 1e9, te, tr);
                        last_prog = tn;
                    }
                    if (stop_req || (ckpt_every > 0 && tn - last_ckpt >= ckpt_every)) {
                        if (save_ckpt(2, next_expand, nprod, tn - t0, nW, Wrows, Wlen, w, c + 1) == 0) fprintf(stderr, "[checkpoint] %s written (%.2f GB)%s\n", ckpt_path, mem_rows / 1e9, stop_req ? " -- interrupted; rerun the same command to resume" : "");
                        last_ckpt = now();
                        if (stop_req) return 3;
                    }
                }
            }
        }
        if (B > 0) {
            reduce_batch(batch, B);
            int np = self_eliminate(batch, B, piv);
            if (np > 0) {
                ensure_rows(nrows + np);
                if (nchunks >= cap_chunks) { cap_chunks = cap_chunks ? cap_chunks * 2 : 1024; chunk_end = realloc(chunk_end, cap_chunks * sizeof(int)); }
                for (int j = 0; j < np; ++j) { u64 *q = batch + (size_t)piv[j] * W; int h = highest_bit(q); store_row(nrows, q, h); pivcol[nrows] = h; is_pivot_col[h] = 1; nrows++; }
                chunk_end[nchunks++] = (int)nrows;
            }
        }
        total_products = nprod;
        printf("MACAULAY (stage 1) ");
    }
    if (ckpt_path[0]) remove(ckpt_path);
    if (getenv("WD_DEBUG")) { bucket_rows(); for (int d = 0; d <= D; ++d) fprintf(stderr, "  deg %d: rows %ld expanded %ld\n", d, deg_len[d], qpos[d]); }
    /* profile */
    long *cnt = calloc(D + 1, sizeof(long));
    for (long i = 0; i < nrows; ++i) cnt[mono_deg[pivcol[i]]]++;
    long cum = 0;
    printf("N=%d D=%d dim(V_D)=%ld  products=%ld  time=%.0fs  rows_mem=%.2fGB\n", N, D, nrows, total_products, (now() - t0), mem_rows / 1e9);
    for (int d = 0; d <= D; ++d) { cum += cnt[d]; printf("  dim(V_%d cap B_<=%d) = %ld\n", D, d, cum); }
    /* non-pivot columns = standard monomials (deg <= D) w.r.t. LM(V_D); if all have degree <= D-1 and their number equals #solutions, the reduced GB lies in V_D */
    long *np = calloc(D + 1, sizeof(long)); long nptot = 0;
    for (int c = 0; c < M; ++c) if (!is_pivot_col[c]) { np[mono_deg[c]]++; nptot++; }
    printf("  non-pivot monomials: total %ld  by degree:", nptot);
    for (int d = 0; d <= D; ++d) printf(" %d:%ld", d, np[d]);
    printf("\n");
    return 0;
}
