"""Rigorous existence analysis of uniformly most split reliable (UMSR) two-terminal graphs.

Input: splitrel output (one line per (graph6, s, t): m and N_{n-2..m}).
For each m:
  EXISTS  if some vector equals the componentwise maximum (then that two-terminal graph is UMSR,
          since SR(p) = sum N_i p^i (1-p)^{m-i} with nonnegative basis on [0,1]).
  NONE    certified if the set of vectors attaining max N_{n-2} (needed for optimality near p=0)
          is disjoint from the set attaining the lexicographic maximum of (N_{m-1}, N_{m-2}, ...)
          (needed for optimality near p=1).  Both are necessary conditions for optimality on all of
          [0,1], so disjointness proves that no UMSR graph exists.
Writes a JSON summary and prints a table.
"""
import sys, json
from collections import defaultdict
import networkx as nx

fn, outjson = sys.argv[1], sys.argv[2]
bym = defaultdict(dict)   # m -> vector -> list of (g6, s, t)
n = None
with open(fn) as f:
    for line in f:
        p = line.split()
        g6, s, t, m = p[0], int(p[1]), int(p[2]), int(p[3])
        vec = tuple(int(x) for x in p[4:])
        if n is None:
            n = nx.from_graph6_bytes(g6.encode()).number_of_nodes()
        bym[m].setdefault(vec, []).append((g6, s, t))

def describe(g6, s, t):
    G = nx.from_graph6_bytes(g6.encode())
    deg = dict(G.degree())
    return {"graph6": g6, "s": s, "t": t, "deg_s": deg[s], "deg_t": deg[t],
            "complement_edges": sorted(map(list, nx.complement(G).edges())),
            "bridges": len(list(nx.bridges(G)))}

summary = {"n": n, "results": {}}
print("n = %d" % n)
for m in sorted(bym):
    vecs = list(bym[m])
    L = len(vecs[0])
    mx = tuple(max(v[j] for v in vecs) for j in range(L))
    winners = [v for v in vecs if v == mx]
    k = n * (n - 1) // 2 - m
    if winners:
        reps = [describe(*bym[m][winners[0]][0])]
        for (g6, s, t) in bym[m][winners[0]][1:]:
            if g6 not in [r["graph6"] for r in reps]:
                reps.append(describe(g6, s, t))
        summary["results"][m] = {"m": m, "k": k, "exists": True, "max_vector": list(mx),
                                 "umsr_two_terminal_graphs": reps}
        print("m=%2d k=%2d EXISTS  winners: %s" % (m, k, ", ".join("%s(%d,%d)" % (r["graph6"], r["s"], r["t"]) for r in reps)))
    else:
        bottom = max(v[0] for v in vecs)            # max N_{n-2}
        S0 = {v for v in vecs if v[0] == bottom}
        toplex = max(vecs, key=lambda v: tuple(reversed(v)))
        S1 = {v for v in vecs if tuple(reversed(v)) == tuple(reversed(toplex))}
        certified = S0.isdisjoint(S1)
        cand = describe(*bym[m][toplex][0])
        beater = describe(*bym[m][next(iter(S0))][0])
        summary["results"][m] = {"m": m, "k": k, "exists": False, "certified_nonexistence": certified,
                                 "max_vector": list(mx), "top_lex_candidate": cand,
                                 "top_lex_vector": list(toplex), "max_N_{n-2}": bottom,
                                 "max_N_{n-2}_attained_by": beater}
        print("m=%2d k=%2d NONE    certified=%s  top-lex cand %s(%d,%d) N_{n-2}=%d < max %d (e.g. %s(%d,%d))" % (
            m, k, certified, cand["graph6"], cand["s"], cand["t"], toplex[0], bottom, beater["graph6"], beater["s"], beater["t"]))
ex = [m for m in sorted(bym) if summary["results"][m]["exists"]]
print("SUMMARY n=%d: UMSR exists for m in %s" % (n, ex))
print("   predicted (m = n-1 or m >= C(n,2)-n): %s" % ([n - 1] + list(range(n * (n - 1) // 2 - n, n * (n - 1) // 2 + 1))))
json.dump(summary, open(outjson, "w"), indent=1)
