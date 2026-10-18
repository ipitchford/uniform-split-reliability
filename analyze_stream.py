"""Streaming version of analyze2.py for large splitrel outputs (two passes, O(#m) memory).

Same logic as analyze2.py:
  EXISTS  if some vector equals the coefficientwise maximum over the class;
  NONE    certified if the vector attaining the lexicographic maximum of (N_m, N_{m-1}, ...)
          has N_{n-2} strictly below the maximum of N_{n-2} over the class.
Usage: python3 analyze_stream.py splitrel_output.txt out.json
"""
import sys, json
from collections import defaultdict
import networkx as nx

fn, outjson = sys.argv[1], sys.argv[2]

def parse(line):
    p = line.split()
    return p[0], int(p[1]), int(p[2]), int(p[3]), tuple(int(x) for x in p[4:])

# pass 1: per m, componentwise max, max N_{n-2}, top-lex vector (+ one attaining line)
mx = {}; bottom = {}; toplex = {}; toplex_who = {}; counts = defaultdict(int); n = None
with open(fn) as f:
    for line in f:
        g6, s, t, m, vec = parse(line)
        if n is None:
            n = nx.from_graph6_bytes(g6.encode()).number_of_nodes()
        counts[m] += 1
        if m not in mx:
            mx[m] = list(vec); bottom[m] = vec[0]; toplex[m] = vec; toplex_who[m] = (g6, s, t)
            continue
        cur = mx[m]
        for j, x in enumerate(vec):
            if x > cur[j]: cur[j] = x
        if vec[0] > bottom[m]: bottom[m] = vec[0]
        if tuple(reversed(vec)) > tuple(reversed(toplex[m])):
            toplex[m] = vec; toplex_who[m] = (g6, s, t)

# pass 2: collect winners (vector == max) and one line attaining max N_{n-2}
winners = defaultdict(list); bottom_who = {}
mxt = {m: tuple(v) for m, v in mx.items()}
with open(fn) as f:
    for line in f:
        g6, s, t, m, vec = parse(line)
        if vec == mxt[m]:
            if all(w[0] != g6 for w in winners[m]):
                winners[m].append((g6, s, t))
        if vec[0] == bottom[m] and m not in bottom_who:
            bottom_who[m] = (g6, s, t)

def describe(g6, s, t):
    G = nx.from_graph6_bytes(g6.encode())
    deg = dict(G.degree())
    return {"graph6": g6, "s": s, "t": t, "deg_s": deg[s], "deg_t": deg[t],
            "complement_edges": sorted(map(list, nx.complement(G).edges())),
            "bridges": len(list(nx.bridges(G)))}

summary = {"n": n, "results": {}}
print("n = %d" % n)
for m in sorted(mx):
    k = n * (n - 1) // 2 - m
    if winners[m]:
        reps = [describe(*w) for w in winners[m]]
        summary["results"][m] = {"m": m, "k": k, "pairs": counts[m], "exists": True, "max_vector": mxt[m],
                                 "umsr_two_terminal_graphs": reps}
        print("m=%2d k=%2d EXISTS  (%d two-terminal graphs) winners: %s" % (m, k, counts[m],
              ", ".join("%s(%d,%d)" % (r["graph6"], r["s"], r["t"]) for r in reps)))
    else:
        certified = toplex[m][0] < bottom[m]
        cand = describe(*toplex_who[m]); beater = describe(*bottom_who[m])
        summary["results"][m] = {"m": m, "k": k, "pairs": counts[m], "exists": False,
                                 "certified_nonexistence": certified, "max_vector": mxt[m],
                                 "top_lex_candidate": cand, "top_lex_vector": list(toplex[m]),
                                 "max_N_{n-2}": bottom[m], "max_N_{n-2}_attained_by": beater}
        print("m=%2d k=%2d NONE    certified=%s  top-lex cand %s(%d,%d) N_{n-2}=%d < max %d (e.g. %s(%d,%d))" % (
            m, k, certified, cand["graph6"], cand["s"], cand["t"], toplex[m][0], bottom[m],
            beater["graph6"], beater["s"], beater["t"]))
ex = [m for m in sorted(mx) if summary["results"][m]["exists"]]
print("SUMMARY n=%d: UMSR exists for m in %s" % (n, ex))
print("   predicted (m = n-1 or m >= C(n,2)-n): %s" % ([n - 1] + list(range(n * (n - 1) // 2 - n, n * (n - 1) // 2 + 1))))
json.dump(summary, open(outjson, "w"), indent=1)
