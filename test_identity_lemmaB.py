"""(1) Verify N_i(H) = c_i(H/st) - c_i(H) on all connected n-vertex graphs.
(2) Test Lemma B': c(H/st) <= c(G/st) coefficientwise for all H with k missing edges, k<=n-2,
    where G = K_n - S_k(s), t a leaf; and c(H) >= c(G) (PSS least-reliable)."""
import sys
from collections import Counter
import networkx as nx
from splitlib import split_vector, dominates
from multigraph_lib import conn_poly_multigraph, identify, simple_mult
n = int(sys.argv[1]); kmax = int(sys.argv[2]); g6file = sys.argv[3]

def star_complement(n, k, s=0, t=1):
    adj = [0]*n
    leaves = ([t] + [v for v in range(n) if v not in (s, t)][:k-1]) if k >= 1 else []
    for u in range(n):
        for v in range(u+1, n):
            if (u == s and v in leaves) or (v == s and u in leaves): continue
            adj[u] |= 1<<v; adj[v] |= 1<<u
    return adj
Gc = {}; Gcst = {}
for k in range(0, n-1):
    adj = star_complement(n, k)
    Gc[k] = conn_poly_multigraph(n, simple_mult(n, adj))
    Gcst[k] = conn_poly_multigraph(n-1, identify(n, adj, 0, 1))
stats = Counter(); fails = {}
checked = 0
for line in open(g6file):
    G = nx.from_graph6_bytes(line.strip().encode())
    m = G.number_of_edges(); k = n*(n-1)//2 - m
    if k > kmax: continue
    adj = [0]*n
    for u, v in G.edges(): adj[u] |= 1<<v; adj[v] |= 1<<u
    cH = conn_poly_multigraph(n, simple_mult(n, adj))
    for s in range(n):
        for t in range(s+1, n):
            cHst = conn_poly_multigraph(n-1, identify(n, adj, s, t))
            # identity check
            sv = split_vector(adj, s, t)
            diff = [a - b for a, b in zip(cHst, cH)]
            if (adj[s]>>t)&1:   # st edge becomes a loop: c_i(X)+c_{i-1}(X)
                X = cHst + [0]
                loop = [X[i] + (X[i-1] if i >= 1 else 0) for i in range(len(X))]
                diff = [a - b for a, b in zip(loop, cH)]
            assert diff == sv, (line, s, t, diff, sv)
            checked += 1
            if 1 <= k <= n-2:
                Xl = cHst if not (adj[s]>>t)&1 else loop
                L = max(len(Xl), len(Gcst[k]))
                okB = dominates(Gcst[k] + [0]*(L-len(Gcst[k])), Xl + [0]*(L-len(Xl)))
                okP = dominates(cH, Gc[k])
                stats[("k=%d" % k, "B'", okB)] += 1
                stats[("k=%d" % k, "PSS", okP)] += 1
                if not okB and ("B'", k) not in fails: fails[("B'", k)] = (line.strip(), s, t, cHst, Gcst[k])
                if not okP and ("PSS", k) not in fails: fails[("PSS", k)] = (line.strip(), s, t, cH, Gc[k])
print("identity checked on", checked, "two-terminal graphs")
for key in sorted(stats, key=str): print(key, stats[key])
for c, f in fails.items(): print(c, f)
