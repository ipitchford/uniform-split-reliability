"""Empirical test of candidate monotone moves for split reliability in the dense range.

Move M1 ("concentrate a missing edge at s"): H' = H - sx + vw with x in N(s), deg(s) >= 2,
vw a non-edge, s not in {v,w}.  Checks whether N(H') >= N(H) coefficientwise, broken
down by the position of t.
"""
import sys
from collections import Counter
import networkx as nx
from splitlib import split_vector, adj_from_edges, complement_edges, dominates

n = int(sys.argv[1]); kmax = int(sys.argv[2])
g6file = sys.argv[3]
stats = Counter(); fails = []
cache = {}

def vec(adj, s, t):
    key = (tuple(adj), s, t)
    if key not in cache:
        cache[key] = split_vector(adj, s, t)
    return cache[key]

for line in open(g6file):
    G = nx.from_graph6_bytes(line.strip().encode())
    m = G.number_of_edges(); k = n * (n - 1) // 2 - m
    if k > kmax:
        continue
    adj = [0] * n
    for u, v in G.edges():
        adj[u] |= 1 << v; adj[v] |= 1 << u
    nonedges = complement_edges(adj)
    for s in range(n):
        if bin(adj[s]).count("1") < 2:
            continue
        for t in range(n):
            if t == s:
                continue
            base = vec(adj, s, t)
            for x in range(n):
                if not (adj[s] >> x) & 1:
                    continue
                for (v, w) in nonedges:
                    if s in (v, w):
                        continue
                    adj2 = list(adj)
                    adj2[s] &= ~(1 << x); adj2[x] &= ~(1 << s)
                    adj2[v] |= 1 << w; adj2[w] |= 1 << v
                    # connectivity check
                    H2 = nx.Graph(); H2.add_nodes_from(range(n))
                    H2.add_edges_from([(a, b) for a in range(n) for b in range(a + 1, n) if (adj2[a] >> b) & 1])
                    if not nx.is_connected(H2):
                        stats[("disconnected",)] += 1
                        continue
                    new = vec(adj2, s, t)
                    cat = ("x==t" if x == t else "x!=t",
                           "t in vw" if t in (v, w) else "t notin vw",
                           "st edge" if (adj[s] >> t) & 1 else "st nonedge")
                    ok = dominates(new, base)
                    stats[(cat, ok)] += 1
                    if not ok and len(fails) < 20:
                        fails.append((line.strip(), s, t, x, (v, w), base, new))
for key in sorted(stats, key=str):
    print(key, stats[key])
print("example failures:")
for f in fails[:10]:
    print(f)
