"""Check Lemma B' at k = n-1 and k = n: c(H/st) (loop included) <= c(K_{n-1}) resp. c(K_{n-1}-e)."""
import sys, networkx as nx
from collections import Counter
from splitlib import dominates
from multigraph_lib import conn_poly_multigraph, identify, simple_mult
n = int(sys.argv[1]); g6file = sys.argv[2]
p = n-1
K = [[1 if u!=v else 0 for v in range(p)] for u in range(p)]
cK = conn_poly_multigraph(p, K)
Ke = [row[:] for row in K]; Ke[0][1]=Ke[1][0]=0
cKe = conn_poly_multigraph(p, Ke)
stats = Counter()
for line in open(g6file):
    G = nx.from_graph6_bytes(line.strip().encode()); m = G.number_of_edges(); k = n*(n-1)//2 - m
    if k not in (n-1, n): continue
    adj = [0]*n
    for u,v in G.edges(): adj[u] |= 1<<v; adj[v] |= 1<<u
    target = cK if k == n-1 else cKe
    for s in range(n):
        for t in range(s+1, n):
            X = conn_poly_multigraph(p, identify(n, adj, s, t))
            if (adj[s]>>t)&1:
                X = X + [0]; X = [X[i] + (X[i-1] if i else 0) for i in range(len(X))]
            L = max(len(X), len(target))
            ok = dominates(target + [0]*(L-len(target)), X + [0]*(L-len(X)))
            stats[(k, ok)] += 1
print(stats)
