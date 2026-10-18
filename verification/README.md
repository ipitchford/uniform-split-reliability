# Exact replay and regression checks

From the repository root, with Python 3 and a C compiler available:

```sh
bash verification/reproduce.sh /tmp/umsr-full-replay
```

The output directory must not already exist. This compiles `splitrel.c` from
source and runs every check in both ordinary and optimized (`python -O`) mode.
It uses the Python standard library only; it never executes a bundled binary.
For bounded CI smoke testing, append `--quick`. That limits the cone checks to
N=7 through 10 and **does not replay the full finite theorem range**.

## What is checked

- `replay_cones.py`: replay of the supplied exact integer algorithm for every
  required cone triple with N=7 through 50, both epsilon values. The full run
  checks 2,024 triples and 1,354,694 scalar coefficient inequalities. It also
  checks the first excluded parameter and explicitly retains the five
  small exceptions. This is a replay, not an independent implementation.
- `independent_sanity.py`: independently implemented ordinary coefficient-list
  recurrence for the excluded cases at N=7 through 10; integer matrix-tree
  checks of the zero-star formulas for n=7 through 30; direct edge-subset checks
  of the terminal-identification identity on complete graphs, paths and cycles
  of orders 3 through 5. All 57 terminal cases pass when the loop factor is
  retained; 40 deliberate loop omissions are rejected.
- `crosscheck_roots.py`: brute edge-subset checks against the supplied packed
  root-set algorithm for clique orders 1 through 6 (21 root-set pairs and
  33,867 edge subsets), and disconnected zero-neighbour cone controls.
- `base_cases.py`: fresh compiled-C replay of 5,422 supplied graph
  representatives at the two boundary edge counts for n=7 through 10. Checks
  every coefficientwise minimum and the balloon property of every minimizer.
  **Completeness of the graph corpora is inherited from their supplied
  enumeration provenance; this script does not independently enumerate
  unlabeled graphs.** Corpus hashes are recorded in its receipts.
- `test_cli.py`: six regression tests covering valid K3 output; malformed,
  truncated, padded, header-prefixed and empty graph6 input; unknown flags;
  mixed graph orders in maximum-only mode; empty cone ranges; and deliberately
  injected comparison failure and missing coverage. These errors must produce
  nonzero exits, including under optimized Python.
- `check_interval_algebra.py`: exact rational spot checks of the interval
  proof's factorization and margin, cut constants, and tree ratios through
  N=100. These finite checks are supplementary; the universal argument is
  written in the manuscript, not established by the scan.

## Receipts

The new output directory contains `cli-normal.log`, `cli-optimized.log`,
`independent-{normal,optimized}.json`, `rootsets-{normal,optimized}.json`,
`base-{normal,optimized}.json`, and `cones/cone-{normal,optimized}.json`, plus
`interval-algebra-{normal,optimized}.log`,
progress logs, raw base-case inputs/outputs, and the freshly compiled executable.
Cone receipts bind the algorithm source hash and Python version. Progress files
are not completion receipts. A failed run may leave diagnostic partial output;
only a zero exit from the entire script establishes completion of its stated
scope. Output directories should remain outside the source archive.

These checks do not establish an unbounded boundary result, novelty, formal
verification, or external peer review. Small independent checks are distinct
from the larger supplied-algorithm replay. The older top-level exploratory
scripts are not the release's authoritative fail-closed verification entrypoint.
