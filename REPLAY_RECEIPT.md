# Development checks and final replay route

On 18 October 2026 the revised source passed the bounded verification route in
ordinary and optimized Python 3.14.7 on ARM macOS, using Apple clang 21.0.0
(clang-2100.3.34.2) to compile the C source. No supplied binary was executed.

Retained development receipts in `verification/receipts/` establish:

- 5,422 supplied boundary-case graph representatives at n=7 through 10,
  recomputed with the hardened C parser;
- independent integer matrix-tree checks of the zero-star formulas at n=7
  through 30, the five first-excluded exceptions, and 57 loop-aware identities
  with 40 loop-omission negative controls;
- 21 root-set cases, checked by 33,867 brute-force edge subsets;
- six CLI regression tests, including malformed/truncated/NUL/overlong input,
  deliberately failed comparison and missing coverage;
- exact interval-algebra checks on 31,920 rational cases, cut constants and
  target tree ratios through N=100;
- bounded cone replay through N=10 in both modes.

These retained development receipts are **not** the final full-archive replay
receipt. The final publication gate runs `bash verification/reproduce.sh NEW_DIR`
from a fresh extraction, including N=7 through 50, in both modes, and retains
the exact archive-bound gate externally. A fresh run is reproducible from the
included code; the quick mode alone cannot establish the finite n<=51 claim.
The full required range has 2,024 triples and 1,354,694 scalar inequalities.

All receipt claims have the limitations in `verification/README.md` and
`ASSURANCE.md`. In particular, graph-corpus completeness is inherited, finite
algebra tests are not universal proofs, and alternate implementations in the
same production process are not unaffiliated reproduction.
