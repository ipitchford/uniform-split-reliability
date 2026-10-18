# AI index: uniform split reliability

Read [STATUS.md](STATUS.md) first. The authoritative statements and proofs are in
[paper/split_reliability.md](paper/split_reliability.md). This unrefereed candidate
does not solve the unrestricted-order boundary conjecture.

## Claim map

| ID | Statement and manuscript locator | Evidence and dependencies |
|---|---|---|
| UMSR-EXISTENCE | All-order classification outside two edge counts; main theorem, nonexistence/dense-existence sections | Written competitor and exchange proofs; Romero's candidate theorem and dense least-reliability theorem |
| UMSR-FINITE51 | Complete classification for 7<=n<=51; boundary reduction and computation sections | `cone_star2.py`, `verification/replay_cones.py`, base cases plus written reduction |
| UMSR-IDENTITY | Split equals identified minus original reliability; preliminaries | Written proof and loop-aware `verification/independent_sanity.py` |
| UMSR-INTERVALS | All N>=7 for p>=1/2; explicit small-p interval in target family; analytic interval theorem | Analytic cut/tree bounds and displayed R,L; finite tests supplementary |
| UMSR-EXCEPTIONS | Five first-excluded triples satisfy coefficientwise comparison; small certificate | `verification/replay_cones.py`, small-case independent enumeration |
| UMSR-OPEN | Two edge counts remain conditional at unrestricted order; discussion | Unresolved middle probability interval; no full-interval proof |

## Exact meanings and traps

`n` is original graph order; cone order is `N=n-1` in the reduction. `k=binom(n,2)-m`.
`p` is retention, `q=1-p`. `C_N(s)` is a clique on N-1 vertices plus one vertex
joined to s clique vertices. Ordinary balloon B and two-terminal balloon G differ.
The sufficient cone inequality is
`(1-q^(d+1)) Rel(C_N(h-d);p) >= p Rel(C_N(h);p)`, with
`h=N-2-epsilon`, epsilon in {0,1}, and `1<=d<=N-5-epsilon`.
It is a sufficient route, not an asserted equivalence to the whole conjecture.

## Replay and adversarial checks

See [REPLAY.md](REPLAY.md) and [ENVIRONMENT.txt](ENVIRONMENT.txt). The cone algorithm
packs integer polynomial coefficients. `verification/crosscheck_roots.py` compares
root-set counts with edge-subset enumeration. `verification/replay_cones.py`
checks exact parameter coverage and fails on corruption/missing cases. Rebuild
`splitrel.c`; do not trust a supplied binary. Respect its supported integer range.

Challenge terminal-loop conventions, degree elevation, sufficient versus
necessary reductions, strict versus weak exchange, zero-star construction,
required versus excluded parameters, and any all-order claim crossing the
unresolved interval. Historical `logs/` and `results/` are not automatically fresh
replay. Quick CI samples do not replace the full finite theorem. Corpus completeness
and arithmetic replay are separate obligations. Universal proof, source inspection,
historical computation, fresh checking and external review remain distinct.

Future work should target the remaining interval or improve the sufficient bound;
raising a finite cutoff does not solve the conjecture. Licence map: [LICENSES.md](LICENSES.md).
