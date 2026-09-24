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

typedef uint64_t u64;
static int N, D, W, M;
static int *deg_start;           /* first column index of each degree */
static u64 *mono_mask;           /* column -> monomial bitmask */
static uint8_t *mono_deg;
static u64 binom[64][8];
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

static void reduce_batch(u64 *batch, int B) {
    int cs = 0;
    for (int c = 0; c < nchunks; ++c) {
        int ce = chunk_end[c];
        for (int b0 = cs; b0 < ce; b0 += K) {
            int kk = ce - b0 < K ? ce - b0 : K;
            int cols[16];
            for (int j = 0; j < kk; ++j) cols[j] = pivcol[b0 + j];
            /* Gray-code table: table[p] = XOR of rows with bits of p */
            memset(table, 0, (size_t)W * sizeof(u64));
            for (int p = 1; p < (1 << kk); ++p) {
                int j = __builtin_ctz(p);
                memcpy(table + (size_t)p * W, table + (size_t)(p ^ (1 << j)) * W, (size_t)W * sizeof(u64));
                xorrow_len(table + (size_t)p * W, ROW(b0 + j), rowlen[b0 + j]);
            }
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < B; ++i) {
                u64 *r = batch + (size_t)i * W;
                int p = 0;
                for (int j = 0; j < kk; ++j) p |= getbit(r, cols[j]) << j;
                if (p) xorrow(r, table + (size_t)p * W);
            }
        }
        cs = ce;
    }
}

/* self-elimination inside a batch: returns number of new pivot rows, moved to the front, RREF among themselves */
static int self_eliminate(u64 *batch, int B, int *piv) {
    int np = 0;
    for (int i = 0; i < B; ++i) {
        u64 *r = batch + (size_t)i * W;
        int h = highest_bit(r);
        if (h < 0) continue;
        /* r becomes a pivot with column h; clear h from all other nonzero rows (later rows and earlier pivots) */
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < B; ++j) if (j != i) {
            u64 *s = batch + (size_t)j * W;
            if (getbit(s, h)) xorrow(s, r);
        }
        piv[np] = i; np++;
        (void)h;
    }
    return np;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s input.txt [D] [k] [chunk]\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1], "r"); if (!f) { perror("input"); return 1; }
    int Dfile, n;
    if (fscanf(f, "%d %d", &N, &Dfile) != 2 || fscanf(f, "%d", &n) != 1) { fprintf(stderr, "bad header\n"); return 1; }
    D = argc > 2 ? atoi(argv[2]) : Dfile;
    if (argc > 3) K = atoi(argv[3]);
    int CH = argc > 4 ? atoi(argv[4]) : 4096;
    int MACAULAY = (argc > 5 && strcmp(argv[5], "macaulay") == 0);   /* rank of {t*q : q in W_2, deg t <= D-2} */
    if (N > 62 || D > 7 || K > 14) { fprintf(stderr, "N<=62, D<=7, k<=14\n"); return 1; }
    for (int i = 0; i < 64; ++i) { binom[i][0] = 1; for (int j = 1; j < 8; ++j) binom[i][j] = i ? binom[i-1][j-1] + binom[i-1][j] : 0; }
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
    u64 *batch = malloc((size_t)CH * W * sizeof(u64));
    int *piv = malloc(CH * sizeof(int));
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    fprintf(stderr, "N=%d D=%d monomials=%d words/row=%d  threads=%d  (rows stored to their pivot; worst case %.2f GB)\n", N, D, M, W, nthreads, (double)M * W * 4 / 1e9);

    /* read generators into the first batch */
    int B = 0;
    memset(batch, 0, (size_t)CH * W * sizeof(u64));
    for (int i = 0; i < n; ++i) {
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

    double t0 = now();
    long total_products = 0;
    long next_expand = 0, iter = 0;      /* first stored row not yet expanded */
    for (;;) {
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
        /* --- refill batch with products x_a * (unexpanded rows of degree <= D-1) --- */
        while (next_expand < nrows && B + N <= CH) {
            long i = next_expand++;
            if (mono_deg[pivcol[i]] > (MACAULAY ? 2 : D - 1)) continue;   /* macaulay mode: only the cap-3 closure (to get W_2) */
            expanded[i] = 1;
            const u64 *src = ROW(i); const int slen = rowlen[i];
            u64 *base = batch + (size_t)B * W;          /* N consecutive slots are free (B + N <= CH) */
            #pragma omp parallel for schedule(static)
            for (int a = 0; a < N; ++a) {
                u64 *r = base + (size_t)a * W;
                memset(r, 0, (size_t)W * sizeof(u64));
                const int32_t *mv = mulvar + (size_t)a * M;
                for (int w = 0; w < slen; ++w) {
                    u64 x = src[w];
                    while (x) { int b = __builtin_ctzll(x); x &= x - 1; int c = mv[w * 64 + b]; r[c >> 6] ^= 1ull << (c & 63); }
                }
            }
            for (int a = 0; a < N; ++a) {               /* compact: drop zero products */
                u64 *r = base + (size_t)a * W;
                if (highest_bit(r) >= 0) {
                    if (r != batch + (size_t)B * W) memcpy(batch + (size_t)B * W, r, (size_t)W * sizeof(u64));
                    B++; total_products++;
                }
            }
        }
        if (B == 0 && next_expand >= nrows) break;
        if (++iter % 16 == 0)
            fprintf(stderr, "  chunk %ld: pivots=%ld products=%ld  rows=%.2fGB  %.0fs\n", iter, nrows, total_products, mem_rows / 1e9, (now() - t0));
    }
    if (MACAULAY) {
        /* Stage 1 only: take W_2 = all pivot rows of degree <= 2 of the (full) closure just computed,
           reset the echelon, and insert every product t*w with deg t <= D-2 -- no expansion of falls. */
        long nW = 0; u64 **Wrows = malloc(sizeof(u64*) * (nrows + 1)); int *Wlen = malloc(sizeof(int) * (nrows + 1));
        for (long i = 0; i < nrows; ++i) if (mono_deg[pivcol[i]] <= 2) { Wrows[nW] = rowptr[i]; Wlen[nW] = rowlen[i]; nW++; }
        fprintf(stderr, "macaulay mode: |W_2| = %ld rows of degree <= 2; multipliers of degree <= %d\n", nW, D - 2);
        /* reset echelon (keep W rows alive: they are separate mallocs referenced by Wrows) */
        nrows = 0; nchunks = 0; memset(is_pivot_col, 0, M);
        long nprod = 0; B = 0;
        for (long w = 0; w < nW; ++w) {
            for (int c = 0; c < deg_start[D - 1]; ++c) {          /* all monomials t of degree <= D-2 */
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
                    if ((nprod / (CH - 2)) % 16 == 0) fprintf(stderr, "  macaulay: pivots=%ld products=%ld  %.0fs\n", nrows, nprod, now() - t0);
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
    /* profile */
    long *cnt = calloc(D + 1, sizeof(long));
    for (long i = 0; i < nrows; ++i) cnt[mono_deg[pivcol[i]]]++;
    long cum = 0;
    printf("N=%d D=%d dim(V_D)=%ld  products=%ld  time=%.0fs  rows_mem=%.2fGB\n", N, D, nrows, total_products, (now() - t0), mem_rows / 1e9);
    for (int d = 0; d <= D; ++d) { cum += cnt[d]; printf("  dim(V_%d cap B_<=%d) = %ld\n", D, d, cum); }
    return 0;
}
