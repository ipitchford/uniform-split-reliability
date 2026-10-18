# Uniformly most split reliable graphs

Anonymous · v1.0.0 · unrefereed mathematical candidate · 18 October 2026

The [paper](paper/split_reliability.md) studies connected simple graphs with two
distinct terminals. Each edge survives independently with common probability p.
Success means exactly two surviving components, one containing each terminal.
A uniformly most split reliable graph maximises this probability for every p in
[0,1], among graphs with the same vertex and edge counts.

## Results and limits

For n >= 7 vertices, a uniform optimiser exists for m = n-1 (the previously known
path case); none exists for n <= m <= binom(n,2)-n-1; and the two-terminal balloon
candidate is optimal for m >= binom(n,2)-n+2. The two remaining edge counts,
binom(n,2)-n and one more, remain conditional at unrestricted orders. Exact finite
cone certificates, small-order base cases and the written reduction complete
the classification through n = 51.

The analytic extension proves the sufficient cone inequality at every cone order
N >= 7 on 1/2 <= p <= 1, and on an explicit positive-length interval near zero in
the target parameter range. The intermediate interval remains unresolved in
general. These interval results do not solve the full uniform boundary problem.

This is an existence classification, not an enumeration of all terminal-preserving
isomorphism classes of optimisers. The model is not an empirical reliability claim
about physical networks.

## Reader and agent entry points

* [Manuscript](paper/split_reliability.md): statements, proofs, worked example,
  literature dependencies and remaining problem.
* [AI_INDEX.md](AI_INDEX.md), [CLAIMS.json](CLAIMS.json): claim-to-evidence map.
* [STATUS.md](STATUS.md), [ASSURANCE.md](ASSURANCE.md), [PROVENANCE.md](PROVENANCE.md).
* [REPLAY.md](REPLAY.md), [ENVIRONMENT.txt](ENVIRONMENT.txt): checking commands.
* [SOURCES.md](SOURCES.md), [CITATION_AUDIT.md](CITATION_AUDIT.md), [LICENSES.md](LICENSES.md).

`cone_star2.py` computes exact cone connected-subgraph coefficients using root-set
recursion. `cone_star.py` is another implementation, not unaffiliated reproduction.
`verification/` contains fail-closed checks and small independent enumerations
within the same AI-assisted production process. Original `logs/` and `results/`
are supplied historical evidence; fresh checks are separately identified.

Rebuild `splitrel.c` from source; prebuilt executables are excluded. Graph6 corpora
can be read without nauty, but regeneration requires upstream graph-generation
software. Some large historical raw outputs were not supplied; this is disclosed
in REPLAY.md rather than treated as fresh successful replay. Code execution does
not establish the universal written arguments or novelty.
