import sys, random
sys.path.insert(0, '.')
from wd_generic_dims import instance

def count(eqs, d0, d1):
    """number of F2-solutions of the bilinear system: for each u in F2^d0 the system is affine in w."""
    n = len(eqs); masks = []; linw = []; linu = []; cst = []
    for quad, lin, const in eqs:
        m = [0] * d0; lw = 0; lu = 0
        for (i, j) in quad:
            if i < d0 <= j: m[i] |= 1 << (j - d0)
            elif j < d0 <= i: m[j] |= 1 << (i - d0)
            else: raise ValueError("non-bilinear term")
        for v in lin:
            if v < d0: lu |= 1 << v
            else: lw |= 1 << (v - d0)
        masks.append(m); linw.append(lw); linu.append(lu); cst.append(1 if const else 0)
    total = 0
    for u in range(1 << d0):
        piv = {}; ok = True
        for r in range(n):
            row = linw[r]; x = u; i = 0
            while x:
                if x & 1: row ^= masks[r][i]
                x >>= 1; i += 1
            rhs = (bin(u & linu[r]).count('1') & 1) ^ cst[r]
            while row:
                b = row.bit_length() - 1
                if b in piv: pr, prhs = piv[b]; row ^= pr; rhs ^= prhs
                else: piv[b] = (row, rhs); break
            if row == 0 and rhs: ok = False; break
        if ok: total += 1 << (d1 - len(piv))
    return total

if __name__ == "__main__":
    N, n = int(sys.argv[1]), int(sys.argv[2]); seed = 7000 + 100 * N + n
    eqs, NN = instance(n, N // 2, N - N // 2, random.Random(seed))
    print(f"c_N{N}_n{n}.txt nsol={count(eqs, N // 2, N - N // 2)}", flush=True)
