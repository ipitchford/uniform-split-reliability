"""Numerical test of the crude bound route to (B1') for k = n-1 and k = n.
For q = n-1 and every graph F' on q vertices with f edges (f = delta for k=n-1, f = delta+1 for k=n),
check  S_delta(p) * Rel(K_q - F') >= Rel(K_q - K_{1,c})   (c=1 for k=n-1, c=2 for k=n)
with S_delta(p) = 1 + (1-p) + ... + (1-p)^{delta-1}, on a grid of p in (0,1).
delta ranges over 2..q-2.  Uses all graphs on q vertices via geng (as complements)."""
import sys, subprocess, itertools
import numpy as np
import networkx as nx
from multigraph_lib import conn_poly_multigraph, simple_mult
q = int(sys.argv[1])
def rel_poly(adj_edges, q):
    adj = [0]*q
    for u,v in adj_edges: adj[u] |= 1<<v; adj[v] |= 1<<u
    c = conn_poly_multigraph(q, simple_mult(q, adj))
    m = len(c)-1
    ps = np.linspace(0.001, 0.999, 400)
    return sum(c[i]*ps**i*(1-ps)**(m-i) for i in range(m+1)), ps
allE = [(u,v) for u in range(q) for v in range(u+1,q)]
def comp_rel(F):
    E = [e for e in allE if e not in set(F)]
    return rel_poly(E, q)
R_e, ps = comp_rel([(0,1)])
R_star2, _ = comp_rel([(0,1),(0,2)])
worst = {}
out = subprocess.run(['./bin/geng','-q',str(q)], capture_output=True, text=True).stdout.split()
for g6 in out:
    F = nx.from_graph6_bytes(g6.encode())
    f = F.number_of_edges()
    if f < 2 or f > q-2+1: continue
    Fe = [tuple(sorted(e)) for e in F.edges()]
    R, _ = comp_rel(Fe)
    for case, c, target in (("k=n-1", f, R_e), ("k=n", f-1, R_star2)):
        delta = c
        if delta < 2 or delta > q-2: continue
        S = sum((1-ps)**j for j in range(delta))
        ratio = (S*R/target).min()
        key = (case, delta)
        if key not in worst or ratio < worst[key][0]: worst[key] = (ratio, g6)
for k in sorted(worst): print(k, worst[k])
