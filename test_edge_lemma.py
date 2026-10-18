"""Test Lemma (II): for every two-terminal graph H+ with k-1 missing edges (k<=kmax) and
every non-bridge edge e of H+ not incident to s, the number of i-edge split states of H+
containing e is >= the number of i-edge split states of G^(k-1) = K_n - S_{k-1}(s)
containing a given s-edge sz (z not t), for every i.

N^e(H+) = N(H+) - N(H+ - e);  N^{sz}(G^(k-1)) = N(G^(k-1)) - N(G^(k)).
"""
import sys
from collections import Counter
import networkx as nx
from splitlib import split_vector, dominates

n = int(sys.argv[1]); kmax = int(sys.argv[2]); g6file = sys.argv[3]

def star_complement(n, k, s=0, t=1):
    """K_n minus a star at s with k leaves including t (k>=1); K_n if k==0."""
    adj = [0] * n
    leaves = [t] + [v for v in range(n) if v not in (s, t)][: k - 1] if k >= 1 else []
    for u in range(n):
        for v in range(u + 1, n):
            if (u == s and v in leaves) or (v == s and u in leaves):
                continue
            adj[u] |= 1 << v; adj[v] |= 1 << u
    return adj

Gvec = {k: split_vector(star_complement(n, k), 0, 1) for k in range(0, n - 1)}
Nsz = {k: [a - b for a, b in zip(Gvec[k - 1], Gvec[k])] for k in range(1, n - 1)}
for k in Nsz:
    assert min(Nsz[k]) >= 0

stats = Counter(); fails = {}
cache = {}
def vec(adj, s, t):
    key = (tuple(adj), s, t)
    if key not in cache: cache[key] = split_vector(adj, s, t)
    return cache[key]

for line in open(g6file):
    G = nx.from_graph6_bytes(line.strip().encode())
    m = G.number_of_edges(); km1 = n * (n - 1) // 2 - m   # H+ has k-1 missing edges
    k = km1 + 1
    if k > kmax or k > n - 2:
        continue
    adj = [0] * n
    for u, v in G.edges(): adj[u] |= 1 << v; adj[v] |= 1 << u
    bridges = set(tuple(sorted(e)) for e in nx.bridges(G))
    for s in range(n):
        for t in range(n):
            if s == t: continue
            base = vec(adj, s, t)
            for (v, w) in G.edges():
                if s in (v, w): continue
                if tuple(sorted((v, w))) in bridges: continue
                adj2 = list(adj); adj2[v] &= ~(1 << w); adj2[w] &= ~(1 << v)
                Ne = [a - b for a, b in zip(base, vec(adj2, s, t))]
                cat = ("k=%d" % k, "e at t" if t in (v, w) else "e bulk")
                ok = dominates(Ne, Nsz[k])
                stats[(cat, ok)] += 1
                if not ok and cat not in fails:
                    fails[cat] = (line.strip(), s, t, (v, w), Ne[n - 2:], Nsz[k][n - 2:])
for key in sorted(stats, key=str): print(key, stats[key])
print("example failures:")
for c, f in fails.items(): print(c, f)
