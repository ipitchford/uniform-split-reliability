"""Coefficientwise test of (*): sum_{i=1}^{d+1} C(d+1,i) c_{j-i}(cone_{h-d}) >= c_{j-1}(cone_h)
for all j, where cone_s = K_{N-1} plus a vertex adjacent to s clique vertices, h = N-2-eps."""
import sys
from math import comb
from multigraph_lib import conn_poly_multigraph
def cone(N, s):
    mult = [[1 if u != v else 0 for v in range(N)] for u in range(N)]
    for v in range(1, N):      # vertex 0 = u; keep edges to 1..s only
        if v > s: mult[0][v] = mult[v][0] = 0
    return conn_poly_multigraph(N, mult)
N = int(sys.argv[1])
cache = {}
def c(s):
    if s not in cache: cache[s] = cone(N, s)
    return cache[s]
for eps in (0, 1):
    h = N - 2 - eps
    for d in range(1, h - 2):
        ch, cl = c(h), c(h - d)
        m = len(ch) - 1  # = C(N-1,2)+h ; cl has C(N-1,2)+h-d+1 entries
        ok = True; worst = None
        for j in range(0, m + 2):
            rhs = ch[j - 1] if 0 <= j - 1 <= m else 0
            lhs = sum(comb(d + 1, i) * (cl[j - i] if 0 <= j - i < len(cl) else 0) for i in range(1, d + 2))
            if lhs < rhs:
                ok = False; worst = (j, lhs, rhs); break
        print("N=%d eps=%d d=%d coefficientwise (*) holds: %s %s" % (N, eps, d, ok, worst or ""))
