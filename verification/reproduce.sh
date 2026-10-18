#!/bin/bash
# Usage: bash verification/reproduce.sh NEW_OUTPUT_DIRECTORY [--quick]
# --quick tests cones only through N=10 and cannot certify the n<=51 theorem.
set -euo pipefail
test "$#" -ge 1 && test "$#" -le 2 || { echo 'Expected NEW_OUTPUT_DIRECTORY [--quick]' >&2; exit 2; }
max_n=50
if test "$#" -eq 2; then
  test "$2" = --quick || { echo 'Unknown option' >&2; exit 2; }
  max_n=10
fi
output="$1"
mkdir "$output"
output="$(cd "$output" && pwd)"
cd "$(dirname "$0")/.."
export PYTHONDONTWRITEBYTECODE=1
python_bin="${PYTHON:-python3}"
cc -O2 -Wall -Wextra -o "$output/splitrel" splitrel.c
for mode in normal optimized; do
  flags=(-B)
  if test "$mode" = optimized; then flags=(-O -B); fi
  "$python_bin" "${flags[@]}" -m unittest discover -s verification -p test_cli.py > "$output/cli-$mode.log" 2>&1
  "$python_bin" "${flags[@]}" verification/independent_sanity.py > "$output/independent-$mode.json"
  "$python_bin" "${flags[@]}" verification/crosscheck_roots.py > "$output/rootsets-$mode.json"
  "$python_bin" "${flags[@]}" verification/check_interval_algebra.py > "$output/interval-algebra-$mode.log"
  "$python_bin" "${flags[@]}" verification/base_cases.py --output-dir "$output/base-$mode" --binary "$output/splitrel" > "$output/base-$mode.json"
  "$python_bin" "${flags[@]}" verification/replay_cones.py --output-dir "$output/cones" --max-N "$max_n" > "$output/cone-$mode.log"
done
echo "PASS: normal and optimized checks through cone N=$max_n; receipts in $output"
