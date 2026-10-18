"""Exact verification of inequality (*) for cones over complete graphs.

cone_s = K_M plus a vertex u adjacent to s of the M clique vertices (N = M+1 vertices).
(*)  for h = N-2-eps (eps in {0,1}) and 1 <= d <= h-3:
       sum_{i=1}^{d+1} C(d+1,i) c_{j-i}(cone_{h-d})  >=  c_{j-1}(cone_h)   for every j,
where c_j(G) = number of connected spanning subgraphs of G with j edges.  Coefficientwise (*)
implies the polynomial inequality (1-q^{d+1}) Rel(cone_{h-d}) >= p Rel(cone_h) on [0,1].

The connected-subgraph polynomials are computed by the subset recursion
    c_S(x) = (1+x)^{e(S)} - sum_{T : v0 in T, T proper subset of S} c_T(x) (1+x)^{e(S\\T)}
grouped by orbit type under the symmetry group of cone_s: a vertex subset S is described by
(iota, a, b) = ([u in S], |S cap N(u)|, |S minus N(u) minus u|).  Polynomials are stored as Python
integers (Kronecker substitution, base 2^B) so that products use big-integer multiplication.

Usage: python3 cone_star.py Nmin Nmax [B]
Prints, for each N and eps, the largest d for which coefficientwise (*) holds, and whether it holds
for all 1 <= d <= N-5-eps (the range needed by Theorem 4.6 of the paper).
"""
import sys
from math import comb
from functools import lru_cache

Nmin, Nmax = int(sys.argv[1]), int(sys.argv[2])
B = int(sys.argv[3]) if len(sys.argv) > 3 else None

def run(N, eps):
    M = N - 1; h = N - 2 - eps
    results = {}
    for s in sorted(set([h] + [h - d for d in range(1, h - 2)])):
        results[s] = cone_poly(M, s)
    out = {}
    for d in range(1, h - 2):
        ch, cl = results[h], results[h - d]
        mh = len(ch) - 1
        ok = True
        for j in range(0, mh + 2):
            rhs = ch[j - 1] if 0 <= j - 1 <= mh else 0
            lhs = sum(comb(d + 1, i) * (cl[j - i] if 0 <= j - i < len(cl) else 0) for i in range(1, d + 2))
            if lhs < rhs:
                ok = False; break
        out[d] = ok
    return out

def cone_poly(M, s):
    """Connected-subgraph polynomial (list of ints, index = #edges) of cone_s over K_M."""
    maxe = comb(M, 2) + s
    bits = B if B else (maxe + 8)          # coefficients < 2^(maxe) * (small factor)
    mask = (1 << bits) - 1
    def to_int(coeffs):
        v = 0
        for c in reversed(coeffs): v = (v << bits) | c
        return v
    def to_list(v, deg):
        return [(v >> (bits * i)) & mask for i in range(deg + 1)]
    # (1+x)^e as integer
    @lru_cache(maxsize=None)
    def binpoly(e):
        return to_int([comb(e, j) for j in range(e + 1)])
    def edges(i, a, b):      # edges inside a set of type (iota,a,b)
        return comb(a + b, 2) + (a if i else 0)
    # c[(i,a,b)] as big integers; process by increasing size
    c = {}
    for size in range(1, M + 2):
        for i in (0, 1):
            for a in range(0, s + 1):
                b = size - i - a
                if b < 0 or b > M - s: continue
                e = edges(i, a, b)
                total = binpoly(e)
                acc = 0
                # choose v0: u if i==1, else a neighbour-type vertex if a>=1, else a non-neighbour
                if i == 1:
                    # T contains u: type (1, a', b') with a'<=a, b'<=b, (a',b') != (a,b)
                    for a2 in range(0, a + 1):
                        for b2 in range(0, b + 1):
                            if (a2, b2) == (a, b): continue
                            mult = comb(a, a2) * comb(b, b2)
                            rest = edges(0, a - a2, b - b2)
                            acc += mult * c[(1, a2, b2)] * binpoly(rest)
                elif a >= 1:
                    # T contains a fixed neighbour vertex: type (0, a2, b2) with 1<=a2<=a
                    for a2 in range(1, a + 1):
                        for b2 in range(0, b + 1):
                            if (a2, b2) == (a, b): continue
                            mult = comb(a - 1, a2 - 1) * comb(b, b2)
                            rest = edges(0, a - a2, b - b2)
                            acc += mult * c[(0, a2, b2)] * binpoly(rest)
                else:
                    for b2 in range(1, b + 1):
                        if b2 == b: continue
                        mult = comb(b - 1, b2 - 1)
                        rest = edges(0, 0, b - b2)
                        acc += mult * c[(0, 0, b2)] * binpoly(rest)
                c[(i, a, b)] = total - acc
    full = c[(1, s, M - s)]
    return to_list(full, maxe)

if __name__ == "__main__":
    # self-check against the generic routine for small N
    if Nmin <= 8:
        from multigraph_lib import conn_poly_multigraph
        def cone_generic(N, s):
            mult = [[1 if u != v else 0 for v in range(N)] for u in range(N)]
            for v in range(1, N):
                if v > s: mult[0][v] = mult[v][0] = 0
            return conn_poly_multigraph(N, mult)
        for N in (5, 6, 7):
            for s in range(1, N):
                assert cone_poly(N - 1, s) == cone_generic(N, s), (N, s)
        print("self-check against generic subset recursion passed (N=5,6,7)")
    for N in range(Nmin, Nmax + 1):
        for eps in (0, 1):
            out = run(N, eps)
            need = N - 5 - eps
            ok_needed = all(out.get(d, False) for d in range(1, need + 1))
            maxd = max([d for d, v in out.items() if v] + [0])
            print("N=%2d (n=%2d) eps=%d: coefficientwise (*) holds for d <= %d (needed d <= %d): %s" % (
                N, N + 1, eps, maxd, need, "OK" if ok_needed else "FAIL"), flush=True)
