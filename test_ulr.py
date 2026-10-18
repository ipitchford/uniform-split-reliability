"""Check the uniformly-least-reliable claim c_i(H) >= c_i(B) coefficientwise for all connected
H with n vertices and m edges, m = C(n-1,2)+r for r in a range (r may be 0 or -1), where B is the
balloon (K_{n-2}+cone of degree n-3-? ... built recursively) -- here we take B = the graph among
all with minimum vector and report whether a unique coefficientwise-minimum exists and what it is."""
import sys
from collections import defaultdict
import networkx as nx
from multigraph_lib import conn_poly_multigraph, simple_mult
n = int(sys.argv[1]); g6file = sys.argv[2]; rs = [int(x) for x in sys.argv[3:]]
C = lambda a,b: a*(a-1)//2 if b == 2 else None
targets = {C(n-1,2) + r: r for r in rs}
bym = defaultdict(list)
for line in open(g6file):
    G = nx.from_graph6_bytes(line.strip().encode())
    m = G.number_of_edges()
    if m not in targets: continue
    adj = [0]*n
    for u, v in G.edges(): adj[u] |= 1<<v; adj[v] |= 1<<u
    bym[m].append((tuple(conn_poly_multigraph(n, simple_mult(n, adj))), line.strip()))
for m in sorted(bym):
    rows = bym[m]; L = len(rows[0][0])
    mn = tuple(min(r[0][j] for r in rows) for j in range(L))
    winners = [g6 for vec, g6 in rows if vec == mn]
    print("n=%d m=%d (k=%d): %d graphs; coefficientwise-min attained by %s" % (n, m, n*(n-1)//2-m, len(rows), winners))
    for g6 in winners:
        G = nx.from_graph6_bytes(g6.encode())
        print("   degrees", sorted(d for _, d in G.degree()), "complement edges", sorted(nx.complement(G).edges()))
