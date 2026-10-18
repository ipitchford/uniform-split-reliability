"""For 'stuck' missing bulk edges xy (x,y both non-adjacent to s and to t), test the step
H' = H + xy - sv for every v in N(s): is c(X_{H'}) >= c(X_H)?  Also test the alternative
'direct' route: is N(H) <= N(G) already implied by deg_{V-t}(s) <= n-1-k?"""
import sys
from collections import Counter
import networkx as nx
from splitlib import dominates, split_vector
from multigraph_lib import conn_poly_multigraph, identify
n = int(sys.argv[1]); kmax = int(sys.argv[2]); g6file = sys.argv[3]
stats = Counter(); fails = {}
cache = {}
def cX(adj, s, t):
    key = (tuple(adj), s, t)
    if key not in cache: cache[key] = conn_poly_multigraph(n-1, identify(n, adj, s, t))
    return cache[key]
nstuck = 0
for line in open(g6file):
    G = nx.from_graph6_bytes(line.strip().encode())
    m = G.number_of_edges(); k = n*(n-1)//2 - m
    if k > kmax or k > n-2: continue
    adj = [0]*n
    for u, v in G.edges(): adj[u] |= 1<<v; adj[v] |= 1<<u
    for s in range(n):
        for t in range(n):
            if s == t or (adj[s]>>t)&1: continue
            W = [v for v in range(n) if v not in (s, t)]
            stuck = [(x,y) for x in W for y in W if x < y and not (adj[x]>>y)&1
                     and not (adj[s]>>x)&1 and not (adj[s]>>y)&1 and not (adj[t]>>x)&1 and not (adj[t]>>y)&1]
            if not stuck: continue
            nstuck += 1
            base = cX(adj, s, t)
            degs_nt = bin(adj[s] & ~(1<<t)).count('1')
            direct = degs_nt <= n-1-k
            for (x,y) in stuck:
                for v in W:
                    if not (adj[s]>>v)&1: continue
                    adj2 = list(adj); adj2[x] |= 1<<y; adj2[y] |= 1<<x; adj2[s] &= ~(1<<v); adj2[v] &= ~(1<<s)
                    if adj2[s] == 0: continue
                    ok = dominates(cX(adj2, s, t), base)
                    av = 1 + ((adj[t]>>v)&1)
                    cat = ("k=%d"%k, "a_v=%d"%av, "direct" if direct else "nodirect")
                    stats[(cat, ok)] += 1
                    if not ok and cat not in fails: fails[cat] = (line.strip(), s, t, (x,y), v)
print("two-terminal graphs with a stuck missing edge:", nstuck)
for key in sorted(stats, key=str): print(key, stats[key])
for c, f in fails.items(): print(c, f)
