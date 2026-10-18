"""Exploration: (1) verify the full Kelmans shift is coefficientwise reliability-non-increasing
(all-terminal) on all connected graphs with n vertices; (2) list the connected quasi-threshold
graphs with k = n-1, n missing edges and minimum degree >= 2, their rooted-tree shapes, and the
coefficientwise order of their connectivity vectors relative to the balloon.
Usage: python3 explore_qt.py n g6file
"""
import sys
import networkx as nx
from multigraph_lib import conn_poly_multigraph, simple_mult

n = int(sys.argv[1]); g6file = sys.argv[2]

def adj_of(G):
    adj = [0] * n
    for u, v in G.edges():
        adj[u] |= 1 << v; adj[v] |= 1 << u
    return adj

def cvec(adj):
    return conn_poly_multigraph(n, simple_mult(n, adj))

def is_qt(G):
    for u, v in G.edges():
        Nu = set(G[u]) | {u}; Nv = set(G[v]) | {v}
        if not (Nu <= Nv or Nv <= Nu):
            return False
    return True

def tree_shape(G):
    """Rooted tree of a connected quasi-threshold graph: parent of x = the vertex y != x with
    N[x] subset N[y] and |N[y]| minimal.  Returns depth list and a parenthesised shape."""
    N = {x: set(G[x]) | {x} for x in G}
    parent = {}
    ordk = lambda y: (-len(N[y]), y)
    for x in G:
        cands = [y for y in G if y != x and N[x] <= N[y] and ordk(y) < ordk(x)]
        parent[x] = max(cands, key=ordk) if cands else None
    def sub(x):
        ch = sorted([c for c in G if parent[c] == x], key=lambda c: (-len(N[c]), c))
        return "(" + "".join(sub(c) for c in ch) + ")" if ch else "."
    root = [x for x in G if parent[x] is None][0]
    # depth of every node and a 'spine' description: follow the largest child; list hanging subtrees
    depth = {root: 0}; order = [root]
    for x in order:
        for c in G:
            if parent[c] == x: depth[c] = depth[x] + 1; order.append(c)
    size = {}
    for x in reversed(order):
        size[x] = 1 + sum(size[c] for c in G if parent[c] == x)
    desc = []; x = root
    while True:
        ch = sorted([c for c in G if parent[c] == x], key=lambda c: -size[c])
        if not ch: break
        for c in ch[1:]:
            desc.append("d%d:%s" % (depth[c], sub(c)))
        x = ch[0]
    return sub(root) + "  spine_len=%d hanging=[%s]" % (depth[x] + 1, ", ".join(desc))

mode = sys.argv[3] if len(sys.argv) > 3 else "qt"
if mode == "kelmans":
    bad = 0; total = 0
    for line in open(g6file):
        G = nx.from_graph6_bytes(line.strip().encode())
        adj = adj_of(G); c0 = cvec(adj)
        for a in range(n):
            for b in range(n):
                if a == b: continue
                moved = adj[a] & ~adj[b] & ~(1 << b)
                if not moved: continue
                adj2 = list(adj)
                for w in range(n):
                    if (moved >> w) & 1:
                        adj2[a] &= ~(1 << w); adj2[w] &= ~(1 << a); adj2[b] |= 1 << w; adj2[w] |= 1 << b
                c2 = cvec(adj2); total += 1
                if any(x > y for x, y in zip(c2, c0)):
                    bad += 1
                    if bad <= 3: print("VIOLATION", line.strip(), a, b, c0, c2)
    print("Kelmans full shifts tested:", total, "violations:", bad)
    sys.exit()

rows = []
for line in open(g6file):
    G = nx.from_graph6_bytes(line.strip().encode())
    m = G.number_of_edges(); k = n * (n - 1) // 2 - m
    if k not in (n - 1, n) or not nx.is_connected(G): continue
    if min(d for _, d in G.degree()) < 2 and not (k == n - 1 or k == n): continue
    if not is_qt(G): continue
    rows.append((k, line.strip(), tree_shape(G), min(d for _, d in G.degree()), cvec(adj_of(G))))
for k in (n - 1, n):
    R = [r for r in rows if r[0] == k]
    print("k = %d: %d connected quasi-threshold graphs" % (k, len(R)))
    # balloon = the one with a pendant vertex
    B = [r for r in R if r[3] == 1]
    print("  balloon representatives:", [(r[1], r[2]) for r in B])
    bvec = B[0][4]
    for r in sorted(R, key=lambda r: r[4]):
        dom = all(x >= y for x, y in zip(r[4], bvec))
        print("  %s delta=%d shape=%s  >=balloon:%s  c=%s" % (r[1], r[3], r[2], dom, r[4][n-1:n+6]))
    # poset of coefficientwise domination among the delta>=2 graphs: list covering relations
    R2 = [r for r in R if r[3] >= 2]
    def le(a, b): return all(x <= y for x, y in zip(a[4], b[4]))
    print("  coefficientwise comparabilities among delta>=2 graphs:")
    for a in R2:
        for b in R2:
            if a is b: continue
            if le(a, b) and not any(le(a, c) and le(c, b) and c is not a and c is not b for c in R2):
                print("    %s <= %s" % (a[2], b[2]))
