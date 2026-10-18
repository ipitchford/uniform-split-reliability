"""Connected-spanning-subgraph polynomial of a multigraph and the identity
SR_{s,t}(H) = Rel(H/st) - Rel(H), i.e. N_i(H) = c_i(H/st) - c_i(H)."""
from math import comb


def conn_poly_multigraph(n, mult):
    """mult[u][v] = multiplicity of edge uv (symmetric, zero diagonal).
    Returns list c[0..m] with c[i] = number of connected spanning subgraphs with i edges."""
    full = 1 << n
    ecount = [0] * full
    for S in range(1, full):
        v = (S & -S).bit_length() - 1
        rest = S & ~(1 << v)
        ecount[S] = ecount[rest] + sum(mult[v][u] for u in range(n) if (rest >> u) & 1)
    cpoly = [None] * full
    for S in range(1, full):
        e = ecount[S]
        c = [comb(e, j) for j in range(e + 1)]
        v = (S & -S).bit_length() - 1
        rest = S & ~(1 << v)
        if rest:
            sub = (rest - 1) & rest
            while True:
                T = sub | (1 << v)
                U = S & ~T
                eU = ecount[U]
                for i, ci in enumerate(cpoly[T]):
                    if ci:
                        for k in range(eU + 1):
                            c[i + k] -= ci * comb(eU, k)
                if sub == 0:
                    break
                sub = (sub - 1) & rest
        cpoly[S] = c
    return cpoly[full - 1]


def identify(n, adj, s, t):
    """Multigraph H/st on n-1 vertices: vertices other than s,t keep their order;
    the merged vertex is placed last.  Any st edge is dropped (it becomes a loop)."""
    others = [v for v in range(n) if v not in (s, t)]
    idx = {v: i for i, v in enumerate(others)}
    z = n - 2
    mult = [[0] * (n - 1) for _ in range(n - 1)]
    for u in others:
        for v in others:
            if u < v and (adj[u] >> v) & 1:
                mult[idx[u]][idx[v]] += 1; mult[idx[v]][idx[u]] += 1
        m_uz = ((adj[s] >> u) & 1) + ((adj[t] >> u) & 1)
        mult[idx[u]][z] += m_uz; mult[z][idx[u]] += m_uz
    return mult


def simple_mult(n, adj):
    return [[1 if (u != v and (adj[u] >> v) & 1) else 0 for v in range(n)] for u in range(n)]
