"""Test the crude vertex-conditioning bound for the boundary ULR claim.
Need: (1 - q^{d+1}) Psi(K_{N-1} - K_{1,1+eps+d}) >= p Psi(K_{N-1} - K_{1,1+eps}) on [0,1],
i.e. Rel(K_N - K_{1,s'}) with the cone of degree ... ; equivalently in terms of H = z v Gamma:
  (1-q^{d+1}) Rel(K_N - K_{1,1+eps+d}) >= p Rel(K_N - K_{1,1+eps}),  N = n-1, d = deg_Gamma(v) >= 1.
Checks coefficientwise (in the basis p^i q^{M-i}) and numerically on a grid."""
import sys
from math import comb
from multigraph_lib import conn_poly_multigraph
import numpy as np

def star_minus(N, s):
    mult = [[1 if u != v else 0 for v in range(N)] for u in range(N)]
    for v in range(1, s + 1):
        mult[0][v] = mult[v][0] = 0
    return mult

def rel_poly(N, s, ps):
    c = conn_poly_multigraph(N, star_minus(N, s))
    M = len(c) - 1
    return np.array([sum(c[i] * p**i * (1-p)**(M-i) for i in range(M+1)) for p in ps]), c

n = int(sys.argv[1]); N = n - 1
ps = np.linspace(0.001, 0.999, 999)
for eps in (0, 1):
    R0, _ = rel_poly(N, 1 + eps, ps)
    for d in range(1, N - 2 - eps):
        R1, _ = rel_poly(N, 1 + eps + d, ps)
        lhs = (1 - (1-ps)**(d+1)) * R1; rhs = ps * R0
        ratio = lhs / rhs
        print("n=%d eps=%d d=%d  min ratio LHS/RHS = %.4f at p=%.3f  holds=%s" % (n, eps, d, ratio.min(), ps[ratio.argmin()], (ratio >= 1).all()))
