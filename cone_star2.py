"""Second, faster implementation of the coefficientwise check of inequality (*) (Proposition 4.8(a)).

For the cone graph C_N(s) (vertex u joined to s vertices of K_M, M = N-1),
    sum_j c_j(C_N(s)) x^j = sum_{a=1}^{s} C(s,a) x^a gamma_M(a; x),
where gamma_M(a; x) = sum_l gamma_l x^l and gamma_l is the number of l-edge spanning subgraphs of K_M
in which every component meets a fixed a-set A (choose the a cone edges that are present, then every
component of the clique part must contain one of their endpoints).  Decomposing an arbitrary graph on
M labelled vertices by the union U of its components avoiding A gives the recursion
    gamma_M(a) = (1+x)^{C(M,2)} - sum_{u=1}^{M-a} C(M-a,u) (1+x)^{C(u,2)} gamma_{M-u}(a),
with gamma_a(a) = (1+x)^{C(a,2)}.  O(M^3) big-integer polynomial products in total.
Polynomials with nonnegative integer coefficients are stored as single Python integers (base 2^bits).

Usage: python3 cone_star2.py Nmin Nmax
"""
import sys
from math import comb

def check(N, eps, cache):
    M = N - 1; h = N - 2 - eps
    maxe = comb(M, 2) + M
    bits = maxe + 10
    mask = (1 << bits) - 1
    def to_int(coeffs):
        v = 0
        for c in reversed(coeffs): v = (v << bits) | c
        return v
    def to_list(v, deg):
        return [(v >> (bits * i)) & mask for i in range(deg + 1)]
    key = (M, bits)
    if key not in cache:
        binpoly = {}
        def bp(e):
            if e not in binpoly: binpoly[e] = to_int([comb(e, j) for j in range(e + 1)])
            return binpoly[e]
        gamma = {}
        for a in range(1, M + 1):
            gamma[(a, a)] = bp(comb(a, 2))
            for m in range(a + 1, M + 1):
                acc = 0
                for u in range(1, m - a + 1):
                    acc += comb(m - a, u) * bp(comb(u, 2)) * gamma[(m - u, a)]
                gamma[(m, a)] = bp(comb(m, 2)) - acc
                if gamma[(m, a)] < 0:
                    raise ArithmeticError("Negative root-set polynomial")
        cache[key] = (gamma, to_list, bits)
    gamma, to_list, bits = cache[key]
    def cone(s):
        v = 0
        for a in range(1, s + 1):
            v += comb(s, a) * (gamma[(M, a)] << (bits * a))     # multiply by x^a
        return to_list(v, comb(M, 2) + s)
    def cone_int(s):
        v = 0
        for a in range(1, s + 1):
            v += comb(s, a) * (gamma[(M, a)] << (bits * a))
        return v
    # coefficientwise comparison of ((1+x)^{d+1} - 1) * cone_{h-d}  with  x * cone_h, both as big integers
    # (coefficients stay below 2^bits: the left side is bounded by 2^{d+1} * 2^{C(M,2)+h-d} <= 2^{maxe+1})
    bp_cache = {}
    def bp(e):
        if e not in bp_cache: bp_cache[e] = to_int([comb(e, j) for j in range(e + 1)])
        return bp_cache[e]
    deg = comb(M, 2) + h + 1
    rhs = to_list(cone_int(h) << bits, deg)
    out = {}
    for d in range(1, h - 2):
        lhs = to_list((bp(d + 1) - 1) * cone_int(h - d), deg)
        out[d] = all(l >= r for l, r in zip(lhs, rhs))
    return out, cone

if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: cone_star2.py Nmin Nmax")
    Nmin, Nmax = int(sys.argv[1]), int(sys.argv[2])
    if not 7 <= Nmin <= Nmax <= 50:
        raise SystemExit("Certified replay range: 7 <= Nmin <= Nmax <= 50")
    failed = False
    cache = {}
    if Nmin <= 10:
        # validate the cone polynomials against the generic subset recursion for N <= 10
        from multigraph_lib import conn_poly_multigraph
        def cone_generic(N, s):
            mult = [[1 if u != v else 0 for v in range(N)] for u in range(N)]
            for v in range(1, N):
                if v > s: mult[0][v] = mult[v][0] = 0
            return conn_poly_multigraph(N, mult)
        for N in range(5, 11):
            _, cone = check(N, 0, cache)
            for s in range(1, N):
                if cone(s) != cone_generic(N, s):
                    raise ArithmeticError((N, s, "Generic cone mismatch"))
        print("validated against the generic subset recursion for 5 <= N <= 10, all s", flush=True)
    import time
    for N in range(Nmin, Nmax + 1):
        t0 = time.time()
        for eps in (0, 1):
            out, _ = check(N, eps, cache)
            need = N - 5 - eps
            ok_needed = all(out.get(d, False) for d in range(1, need + 1))
            failed = failed or not ok_needed or set(out) != set(range(1, need + 1))
            maxd = max([d for d, v in out.items() if v] + [0])
            print("N=%2d (n=%2d) eps=%d: coefficientwise (*) holds for d <= %d (needed d <= %d): %s   [%.1f s]" % (
                N, N + 1, eps, maxd, need, "OK" if ok_needed else "FAIL", time.time() - t0), flush=True)
        cache.clear()
    raise SystemExit(1 if failed else 0)
