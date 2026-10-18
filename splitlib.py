"""Pure-Python split-reliability coefficients (exact integers), same method as splitrel.c.

split_vector(adj, s, t) -> list N[0..m] (N[i] = number of split states with i edges),
where adj is a list of neighbour-bitmasks for vertices 0..n-1.
"""
from functools import lru_cache
from math import comb


def split_vector(adj, s, t):
    n = len(adj)
    full = 1 << n
    ecount = [0] * full
    for S in range(1, full):
        v = (S & -S).bit_length() - 1
        rest = S & ~(1 << v)
        ecount[S] = ecount[rest] + bin(adj[v] & rest).count("1")
    m = ecount[full - 1]
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
                cT = cpoly[T]
                for i, ci in enumerate(cT):
                    if ci:
                        for k in range(eU + 1):
                            c[i + k] -= ci * comb(eU, k)
                if sub == 0:
                    break
                sub = (sub - 1) & rest
        cpoly[S] = c
    sr = [0] * (m + 1)
    for A in range(full):
        if not (A >> s) & 1 or (A >> t) & 1:
            continue
        B = (full - 1) & ~A
        ca, cb = cpoly[A], cpoly[B]
        for i, ci in enumerate(ca):
            if ci:
                for k, ck in enumerate(cb):
                    sr[i + k] += ci * ck
    return sr


def adj_from_edges(n, edges):
    adj = [0] * n
    for u, v in edges:
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return adj


def edges_from_adj(adj):
    n = len(adj)
    return [(u, v) for u in range(n) for v in range(u + 1, n) if (adj[u] >> v) & 1]


def complement_edges(adj):
    n = len(adj)
    return [(u, v) for u in range(n) for v in range(u + 1, n) if not (adj[u] >> v) & 1]


def dominates(a, b):
    """a >= b coefficientwise (vectors of equal length)."""
    return all(x >= y for x, y in zip(a, b))
