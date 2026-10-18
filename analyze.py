"""Analyse splitrel output: for each m, does a uniformly most split reliable
two-terminal graph exist (a coefficient vector equal to the componentwise maximum)?

Usage: python3 analyze.py out.txt [--verbose]
Prints one line per m: EXISTS/NONE, the winning (graph6,s,t) list (up to split-
equivalence the vectors coincide), and for NONE the coefficient indices at which
the top candidates lose and to whom.
"""
import sys
from collections import defaultdict
import networkx as nx

fn = sys.argv[1]
verbose = "--verbose" in sys.argv
bym = defaultdict(list)
n = None
with open(fn) as f:
    for line in f:
        parts = line.split()
        g6, s, t, m = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
        vec = tuple(int(x) for x in parts[4:])
        if n is None:
            n = nx.from_graph6_bytes(g6.encode()).number_of_nodes()
        bym[m].append((vec, g6, s, t))

def describe(g6, s, t):
    G = nx.from_graph6_bytes(g6.encode())
    deg = dict(G.degree())
    comp = nx.complement(G)
    return (f"deg(s)={deg[s]} deg(t)={deg[t]} dist={nx.shortest_path_length(G,s,t)} "
            f"bridges={len(list(nx.bridges(G)))} complement_edges={sorted(comp.edges())}")

print(f"n = {n}")
summary = {}
for m in sorted(bym):
    rows = bym[m]
    L = len(rows[0][0])
    mx = tuple(max(r[0][j] for r in rows) for j in range(L))
    winners = [r for r in rows if r[0] == mx]
    if winners:
        summary[m] = True
        reps = {}
        for vec, g6, s, t in winners:
            reps.setdefault(g6, []).append((s, t))
        print(f"m={m}: EXISTS  (vector = {mx})")
        for g6, pairs in reps.items():
            print(f"    {g6} pairs={pairs[:4]}{'...' if len(pairs)>4 else ''}  {describe(g6, *pairs[0])}")
    else:
        summary[m] = False
        # candidates: vectors that are maximal (not dominated)
        vecs = {}
        for vec, g6, s, t in rows:
            vecs.setdefault(vec, []).append((g6, s, t))
        maximal = []
        allv = list(vecs)
        for v in allv:
            dominated = any(w != v and all(w[j] >= v[j] for j in range(L)) for w in allv)
            if not dominated:
                maximal.append(v)
        print(f"m={m}: NONE    max-vector = {mx}; {len(maximal)} maximal vectors")
        # report the lexicographically-from-the-top best (locally most reliable near p=1)
        top = max(maximal, key=lambda v: tuple(reversed(v)))
        g6, s, t = vecs[top][0]
        loses = [n - 2 + j for j in range(L) if top[j] < mx[j]]
        print(f"    top-lex (near p=1) candidate {g6} s={s} t={t}: {describe(g6,s,t)}")
        print(f"      its vector {top}; loses at coefficients N_i for i in {loses}")
        for j in range(L):
            if top[j] < mx[j]:
                beaters = [(g, a, b) for v in maximal if v[j] == mx[j] for (g, a, b) in vecs[v][:1]]
                print(f"      N_{n-2+j}: max {mx[j]} vs {top[j]}, e.g. {beaters[:2]} -> {describe(*beaters[0])}")
        if verbose:
            for v in maximal:
                g6, s, t = vecs[v][0]
                print("      maximal:", v, g6, s, t, describe(g6, s, t))
print("SUMMARY n=%d: exists for m in %s; none for m in %s" % (
    n, [m for m in sorted(summary) if summary[m]], [m for m in sorted(summary) if not summary[m]]))
