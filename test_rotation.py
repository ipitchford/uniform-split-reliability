"""Test the rotation move: H' = H + vw - sv  (vw a non-edge, sv an edge, deg(s)>=2).
Rotating the missing edge vw about v so that it becomes the missing edge sv.
Tabulate monotonicity (N(H') >= N(H) coefficientwise) by role of w and t."""
import sys
from collections import Counter
import networkx as nx
from splitlib import split_vector, complement_edges, dominates
n = int(sys.argv[1]); kmax = int(sys.argv[2]); g6file = sys.argv[3]
stats = Counter(); fails = {}
cache = {}
def vec(adj, s, t):
    key = (tuple(adj), s, t)
    if key not in cache: cache[key] = split_vector(adj, s, t)
    return cache[key]
for line in open(g6file):
    G = nx.from_graph6_bytes(line.strip().encode())
    m = G.number_of_edges(); k = n*(n-1)//2 - m
    if k > kmax: continue
    adj = [0]*n
    for u, v in G.edges(): adj[u] |= 1<<v; adj[v] |= 1<<u
    ne = complement_edges(adj)
    for s in range(n):
        if bin(adj[s]).count('1') < 2: continue
        for t in range(n):
            if t == s: continue
            base = None
            for (v, w) in ne:
                for (vv, ww) in ((v, w), (w, v)):
                    if vv == s or ww == s: continue
                    if not (adj[s] >> vv) & 1: continue
                    adj2 = list(adj)
                    adj2[s] &= ~(1<<vv); adj2[vv] &= ~(1<<s); adj2[vv] |= 1<<ww; adj2[ww] |= 1<<vv
                    H2 = nx.Graph(); H2.add_nodes_from(range(n))
                    H2.add_edges_from([(a,b) for a in range(n) for b in range(a+1,n) if (adj2[a]>>b)&1])
                    if not nx.is_connected(H2): stats[('disc',)] += 1; continue
                    if base is None: base = vec(adj, s, t)
                    new = vec(adj2, s, t)
                    cat = ("v==t" if vv == t else "v!=t", "w==t" if ww == t else "w!=t",
                           "st-edge" if (adj[s]>>t)&1 else "st-nonedge", "k=%d" % k)
                    ok = dominates(new, base)
                    stats[(cat, ok)] += 1
                    if not ok and cat not in fails:
                        fails[cat] = (line.strip(), s, t, vv, ww, base[n-2:], new[n-2:], ne)
for key in sorted(stats, key=str): print(key, stats[key])
print("example failures:")
for c, f in fails.items(): print(c, f)
