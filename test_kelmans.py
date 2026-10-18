"""Test Kelmans shifts a->b (move every edge aw with w in N(a) minus N[b] to bw) for split
reliability monotonicity on dense two-terminal graphs.  Tabulates by the roles of a,b.
"""
import sys
from collections import Counter
import networkx as nx
from splitlib import split_vector, dominates

n = int(sys.argv[1]); kmax = int(sys.argv[2]); g6file = sys.argv[3]
stats = Counter(); fails = {}
cache = {}

def vec(adj, s, t):
    key = (tuple(adj), s, t)
    if key not in cache:
        cache[key] = split_vector(adj, s, t)
    return cache[key]

def role(x, s, t):
    return "s" if x == s else ("t" if x == t else "o")

for line in open(g6file):
    G = nx.from_graph6_bytes(line.strip().encode())
    m = G.number_of_edges(); k = n * (n - 1) // 2 - m
    if k > kmax:
        continue
    adj = [0] * n
    for u, v in G.edges():
        adj[u] |= 1 << v; adj[v] |= 1 << u
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            moved = adj[a] & ~adj[b] & ~(1 << b)
            if moved == 0:
                continue
            adj2 = list(adj)
            for w in range(n):
                if (moved >> w) & 1:
                    adj2[a] &= ~(1 << w); adj2[w] &= ~(1 << a)
                    adj2[b] |= 1 << w; adj2[w] |= 1 << b
            H2 = nx.Graph(); H2.add_nodes_from(range(n))
            H2.add_edges_from([(p, q) for p in range(n) for q in range(p + 1, n) if (adj2[p] >> q) & 1])
            connected = nx.is_connected(H2)
            for s in range(n):
                for t in range(n):
                    if s == t:
                        continue
                    cat = (role(a, s, t), role(b, s, t), "ab-edge" if (adj[a] >> b) & 1 else "ab-nonedge")
                    if not connected:
                        stats[(cat, "disc")] += 1
                        continue
                    ok = dominates(vec(adj2, s, t), vec(adj, s, t))
                    stats[(cat, ok)] += 1
                    if not ok and cat not in fails:
                        fails[cat] = (line.strip(), a, b, s, t, vec(adj, s, t), vec(adj2, s, t))
for key in sorted(stats, key=str):
    print(key, stats[key])
print("example failures per category:")
for c, f in fails.items():
    print(c, f)
