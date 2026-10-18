"""Step claim for Lemma B': for H with k<=n-2 missing edges, st not an edge,
xy a missing edge inside W = V-{s,t}, and v in W with sv in E(H):
H' = H + xy - sv.  Is c(X_{H'}) >= c(X_H) coefficientwise (X = H/st multigraph)?
Also records whether the claim holds when additionally restricting to a_v = 2 (tv also an edge)."""
import sys
from collections import Counter
import networkx as nx
from splitlib import dominates
from multigraph_lib import conn_poly_multigraph, identify
n = int(sys.argv[1]); kmax = int(sys.argv[2]); g6file = sys.argv[3]
stats = Counter(); fails = {}
cache = {}
def cX(adj, s, t):
    key = (tuple(adj), s, t)
    if key not in cache: cache[key] = conn_poly_multigraph(n-1, identify(n, adj, s, t))
    return cache[key]
for line in open(g6file):
    G = nx.from_graph6_bytes(line.strip().encode())
    m = G.number_of_edges(); k = n*(n-1)//2 - m
    if k > kmax or k > n-2: continue
    adj = [0]*n
    for u, v in G.edges(): adj[u] |= 1<<v; adj[v] |= 1<<u
    for s in range(n):
        for t in range(n):
            if s == t or (adj[s]>>t)&1: continue
            base = cX(adj, s, t)
            W = [v for v in range(n) if v not in (s, t)]
            for x in W:
                for y in W:
                    if x >= y or (adj[x]>>y)&1: continue
                    for v in W:
                        if not (adj[s]>>v)&1: continue
                        adj2 = list(adj)
                        adj2[x] |= 1<<y; adj2[y] |= 1<<x; adj2[s] &= ~(1<<v); adj2[v] &= ~(1<<s)
                        if adj2[s] == 0: continue
                        new = cX(adj2, s, t)
                        av = 1 + ((adj[t]>>v)&1)
                        cat = ("k=%d" % k, "a_v=%d" % av, "v in xy" if v in (x,y) else "v notin xy")
                        ok = dominates(new, base)
                        stats[(cat, ok)] += 1
                        if not ok and cat not in fails: fails[cat] = (line.strip(), s, t, (x,y), v, base, new)
for key in sorted(stats, key=str): print(key, stats[key])
for c, f in fails.items(): print(c, f)
