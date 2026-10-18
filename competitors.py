import sys
from collections import defaultdict
import networkx as nx
fn = sys.argv[1]; mwant = int(sys.argv[2]); topn = int(sys.argv[3]) if len(sys.argv) > 3 else 6
rows = []
for line in open(fn):
    p = line.split()
    if int(p[3]) != mwant: continue
    rows.append((tuple(int(x) for x in p[4:]), p[0], int(p[1]), int(p[2])))
L = len(rows[0][0]); mx = tuple(max(r[0][j] for r in rows) for j in range(L))
print("max vector", mx)
# rank competitors by max over j of N_j/mx_j (closeness)
seen = {}
for vec, g6, s, t in rows:
    if vec == mx: continue
    if vec in seen: continue
    seen[vec] = (g6, s, t)
score = sorted(seen.items(), key=lambda kv: -max(kv[0][j]/mx[j] for j in range(L) if mx[j] > 0))
for vec, (g6, s, t) in score[:topn]:
    G = nx.from_graph6_bytes(g6.encode()); deg = dict(G.degree())
    comp = sorted(nx.complement(G).edges())
    ratios = ["%.3f" % (vec[j]/mx[j]) if mx[j] else "-" for j in range(L)]
    print(g6, "s=%d t=%d deg(s)=%d deg(t)=%d" % (s, t, deg[s], deg[t]), "complement", comp)
    print("   ratios N_j/max_j:", ratios)
