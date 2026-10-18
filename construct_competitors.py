"""Explicit constructions used in the nonexistence proof (Theorem A), checked against exact
split-reliability vectors.

For n >= 7 and k = C(n,2) - m with n+1 <= k <= 3n-10 we build
  * the two-terminal balloon candidate G_{n,m} (Romero's locally most split reliable graph), and
  * the competitor H' of the proof,
compute their split vectors with splitlib, and verify
  (1) N_{n-2}(G) equals the closed form of the proof,
  (2) N_{n-2}(H') equals the closed form (A1) or exceeds N_{n-2}(G) (A2),
  (3) for n <= 9, G's vector equals the lexicographic-from-the-top maximum found by exhaustive search
      (i.e. G is the unique candidate), and N_{n-2}(H') > N_{n-2}(G).
Usage: python3 construct_competitors.py n [splitrel_output_file]
"""
import sys, json
from math import comb
from splitlib import split_vector, adj_from_edges

n = int(sys.argv[1]); datafile = sys.argv[2] if len(sys.argv) > 2 else None
C2 = lambda a: a * (a - 1) // 2

def complete_minus(n, missing):
    edges = [(u, v) for u in range(n) for v in range(u + 1, n) if (u, v) not in missing and (v, u) not in missing]
    return adj_from_edges(n, edges)

def tau_star(p, a):            # tau(K_p - K_{1,a}), a >= 1
    return (p - a - 1) * (p - 1) ** (a - 1) * p ** (p - a - 2)
def tau_star_plus_edge(p, a):  # tau(K_p - (K_{1,a} u K_2)), a >= 1, disjoint
    return (p - a - 1) * (p - 1) ** (a - 1) * (p - 2) * p ** (p - a - 3)

# load exhaustive data if given: for each m the top-lex max vector and the max N_{n-2}
data = {}
if datafile:
    from collections import defaultdict
    bym = defaultdict(list)
    for line in open(datafile):
        p_ = line.split(); m_ = int(p_[3]); bym[m_].append(tuple(int(x) for x in p_[4:]))
    for m_, vecs in bym.items():
        data[m_] = (max(vecs, key=lambda v: tuple(reversed(v))), max(v[0] for v in vecs))

rows = []
for k in range(n + 1, 3 * n - 9):
    m = C2(n) - k
    s, t = 0, 1
    if k <= 2 * n - 6:                           # regime (A1): one pendant terminal
        r = k - n + 2                            # 3 <= r <= n-4 missing edges in the core at u
        u = 2; leaves = [t] + list(range(3, 3 + r - 1))        # t = x is a leaf of u's missing star
        missing_G = {(s, v) for v in range(n) if v not in (s, u)} | {(u, l) for l in leaves}
        G = complete_minus(n, missing_G)
        # competitor: missing star at u with r-1 leaves (t among them) plus a disjoint missing edge bc
        leaves2 = [t] + list(range(3, 3 + r - 2)); b, c = n - 2, n - 1
        assert b not in leaves2 and c not in leaves2 and b != u and c != u
        missing_H = {(s, v) for v in range(n) if v not in (s, u)} | {(u, l) for l in leaves2} | {(b, c)}
        H = complete_minus(n, missing_H)
        formula_G = tau_star(n - 1, r - 1); formula_H = tau_star_plus_edge(n - 1, r - 2)
        regime = "A1"
    else:                                        # regime (A2): pendant path of length 2
        rpp = k - 2 * n + 5                      # 0 <= r'' <= n-5 missing edges in the core at u
        w1 = 2; u = 3; x = t
        # vertices: s=0, t=x=1, w1=2, u=3, core = {1,3,...,n-1} (K_{n-2} on t,u,4..n-1)
        core_missing = {(u, t)} | {(u, v) for v in range(4, 4 + rpp - 1)} if rpp >= 1 else set()
        missing_G = {(s, v) for v in range(n) if v not in (s, w1)} | {(w1, v) for v in range(n) if v not in (s, w1, u)} | core_missing
        G = complete_minus(n, missing_G)
        # competitor from the exchange lemma: in Gamma = K_{n-1} - Fhat_G (vertices V minus s; w = w1 plays 'w')
        # Fhat_G = star at w (all except u, x) + star at u (x and r''-1 shared leaves); take a private leaf l'
        # of w (adjacent in Gamma to everything but w) and another leaf a of w; F' = Fhat_G - {aw} + {a l'}.
        shared = list(range(4, 4 + max(rpp - 1, 0)))
        private = [v for v in range(4, n) if v not in shared]   # leaves of w not adjacent to u in F
        assert len(private) >= 2
        lp, a = private[0], private[1]
        Fprime = ({(w1, v) for v in range(n) if v not in (s, w1, u, a)} | {(u, t)} | {(u, v) for v in shared} | {(a, lp)}) if rpp >= 1 else \
                 ({(w1, v) for v in range(n) if v not in (s, w1, u, a)} | {(a, lp)})
        # realise as single pendant: s pendant on w1, terminal t = x (w1 x not in F'), core = K_{n-1} - F' - {w1 x}
        missing_H = {(s, v) for v in range(n) if v not in (s, w1)} | Fprime | {(w1, t)}
        H = complete_minus(n, missing_H)
        formula_G = None; formula_H = None
        regime = "A2"
    NG = split_vector(G, s, t); NH = split_vector(H, s, t)
    mG, mH = sum(bin(a_).count("1") for a_ in G) // 2, sum(bin(a_).count("1") for a_ in H) // 2
    assert mG == m == mH, (mG, m, mH)
    row = {"k": k, "m": m, "regime": regime, "N_{n-2}(G)": NG[n - 2], "N_{n-2}(H')": NH[n - 2],
           "formula_G": formula_G, "formula_H'": formula_H, "H'_beats_G": NH[n - 2] > NG[n - 2]}
    if formula_G is not None:
        assert formula_G == NG[n - 2] and formula_H == NH[n - 2], (formula_G, NG[n - 2], formula_H, NH[n - 2])
    if m in data:
        toplex, maxbottom = data[m]
        row["G_is_toplex_max"] = (tuple(NG[n - 2:]) == toplex)
        row["max_N_{n-2}_exhaustive"] = maxbottom
    rows.append(row)
    print(row)
json.dump({"n": n, "rows": rows}, open("results/competitors_n%d.json" % n, "w"), indent=1)
print("all competitors beat the candidate:", all(r["H'_beats_G"] for r in rows))
if data:
    print("candidate equals exhaustive top-lex max in all cases:", all(r.get("G_is_toplex_max", True) for r in rows))
