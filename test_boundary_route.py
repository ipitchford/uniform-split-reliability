"""For k = n-1 (m = C(n-1,2)) and each two-terminal H with deg(s),deg(t) >= 2, look for a non-edge f
such that N^f(H+f) >= c^e(K_{n-1}) coefficientwise, where N^f(H+f)_i = N_i(H+f) - N_i(H) and
c^e_i(K_{n-1}) = c_i(K_{n-1}) - c_i(K_{n-1}-e).  Report whether some f always works and which."""
import sys, networkx as nx
from collections import Counter
from splitlib import split_vector, dominates
from multigraph_lib import conn_poly_multigraph, simple_mult
n = int(sys.argv[1]); g6file = sys.argv[2]
p = n-1
K = [[1 if u!=v else 0 for v in range(p)] for u in range(p)]; cK = conn_poly_multigraph(p, K)
Ke = [r[:] for r in K]; Ke[0][1]=Ke[1][0]=0; cKe = conn_poly_multigraph(p, Ke)
ce = [a-b for a,b in zip(cK, cKe + [0])]   # c^e_i(K_{n-1}), i = 0..C(p,2)
stats = Counter(); fails = []
for line in open(g6file):
    G = nx.from_graph6_bytes(line.strip().encode()); m = G.number_of_edges(); k = n*(n-1)//2 - m
    if k != n-1: continue
    adj=[0]*n
    for u,v in G.edges(): adj[u]|=1<<v; adj[v]|=1<<u
    deg=[bin(a).count('1') for a in adj]
    for s in range(n):
        for t in range(s+1,n):
            if deg[s] < 2 or deg[t] < 2: continue
            base = split_vector(adj, s, t)
            found = []
            for a in range(n):
                for b in range(a+1,n):
                    if (adj[a]>>b)&1: continue
                    adj2=list(adj); adj2[a]|=1<<b; adj2[b]|=1<<a
                    new = split_vector(adj2, s, t)
                    Nf = [x-y for x,y in zip(new, base+[0])]
                    L = max(len(Nf), len(ce))
                    if dominates(Nf+[0]*(L-len(Nf)), ce+[0]*(L-len(ce))):
                        cat = "f at s or t" if (a in (s,t) or b in (s,t)) else "bulk f"
                        found.append((a,b,cat))
            if found:
                stats[("found", tuple(sorted(set(c for _,_,c in found))))] += 1
            else:
                stats["none"] += 1; fails.append((line.strip(), s, t))
print(stats); print("failures:", fails[:5])
