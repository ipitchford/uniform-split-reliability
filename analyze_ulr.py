"""Uniformly-least-reliable check from `splitrel --conn` output (lines: graph6 m c_0 ... c_m).
For each m: coefficientwise minimum vector over all graphs; report whether it is attained, by which
graphs, and whether every attaining graph is a pendant vertex on a balloon core
(K_{n-1} minus a star K_{1,k-n+2}, k = C(n,2)-m)."""
import sys, json
from collections import defaultdict
import networkx as nx
fn, outjson = sys.argv[1], sys.argv[2]
bym = defaultdict(list); n = None
for line in open(fn):
    p = line.split(); g6 = p[0]; m = int(p[1]); vec = tuple(int(x) for x in p[2:])
    if n is None: n = nx.from_graph6_bytes(g6.encode()).number_of_nodes()
    bym[m].append((vec, g6))
res = {}
for m in sorted(bym):
    k = n*(n-1)//2 - m
    rows = bym[m]; L = len(rows[0][0])
    mn = tuple(min(r[0][j] for r in rows) for j in range(L))
    winners = [g6 for vec, g6 in rows if vec == mn]
    info = []
    ok_all = len(winners) > 0
    for g6 in winners:
        G = nx.from_graph6_bytes(g6.encode())
        pend = [v for v in G if G.degree(v) == 1]
        isballoon = False
        if len(pend) >= 1:
            v = pend[0]; core = G.copy(); core.remove_node(v)
            Fc = nx.complement(core)
            degs = sorted((d for _, d in Fc.degree()), reverse=True)
            # core complement must be a star with k-(n-2) edges (possibly 0 edges)
            e = Fc.number_of_edges()
            isballoon = (e == k - (n - 2)) and (e == 0 or (degs[0] == e and (len(degs) < 2 or degs[1] <= 1)))
        info.append({"graph6": g6, "pendant": len(pend) >= 1, "balloon_core": isballoon})
        ok_all = ok_all and isballoon
    res[m] = {"m": m, "k": k, "graphs": len(rows), "min_attained": len(winners) > 0,
              "all_minimisers_are_pendant_on_balloon_core": ok_all, "minimisers": info}
    print("n=%d m=%d k=%d: %d graphs, min attained=%s, minimisers=%s, all pendant+balloon core=%s" % (
        n, m, k, len(rows), len(winners) > 0, winners, ok_all))
json.dump({"n": n, "results": res}, open(outjson, "w"), indent=1)
