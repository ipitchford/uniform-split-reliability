"""Brute-force reference for split-reliability coefficients (small graphs only).

For each graph6 line read from stdin and each unordered pair {s,t}, enumerate all
2^m edge subsets and count, for each size i, those whose spanning subgraph has
exactly two components with s and t in different components.  Output format
matches splitrel.c so the two can be diffed.
"""
import sys, itertools
import networkx as nx

def components(n, edges):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
    return [find(x) for x in range(n)]

for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    G = nx.from_graph6_bytes(line.encode())
    n, m = G.number_of_nodes(), G.number_of_edges()
    E = list(G.edges())
    # precompute component labels for every subset
    counts = {}
    for mask in range(1 << m):
        sub = [E[i] for i in range(m) if (mask >> i) & 1]
        lab = components(n, sub)
        ncomp = len(set(lab))
        if ncomp != 2:
            continue
        i = bin(mask).count("1")
        for s in range(n):
            for t in range(s + 1, n):
                if lab[s] != lab[t]:
                    counts[(s, t, i)] = counts.get((s, t, i), 0) + 1
    for s in range(n):
        for t in range(s + 1, n):
            vec = [counts.get((s, t, i), 0) for i in range(n - 2, m + 1)]
            print(line, s, t, m, *vec)
