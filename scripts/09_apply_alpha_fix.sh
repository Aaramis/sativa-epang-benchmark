#!/usr/bin/env bash
# Part of 18.sativa_speedup. What the numbers mean: ../RESULTS.md
#
# sativa_epang/sativa.cfg now builds the reference tree under GTRGAMMA, so EPA-ng reads a
# fitted alpha instead of the GTRCAT placeholder. Every EPA-ng cell measured before that
# change describes the old behaviour, so this re-measures them.
#
# Step 1 is timed and runs alone. Step 2 records the previous behaviour under its own name
# (*_catinfo), so the before/after comparison in RESULTS.md stays reproducible.
# One-off: delete this script once the results are regenerated.
set -uo pipefail
cd "$(dirname "$0")/.."
source config.sh
BENCH=(--datasets-json data_ascii/datasets.json --runs-dir results/runs --timeout 10800)

rm -rf results/runs/*_gammainfo__rep*

echo "=== 1. EPA-ng bars, re-measured with the fitted alpha ($(date)) ==="
"$PY3" scripts/02_bench.py "${BENCH[@]}" --datasets n400,n800,n1600,n5402 \
    --conditions port_py3_epang_T1,port_py3_epang_T8,port_py3_epang_T1_K25,port_py3_epang_T8_K25 \
    --repeat 1 --force

echo "=== 2. the previous behaviour, kept as a control ($(date)) ==="
"$PY3" scripts/02_bench.py "${BENCH[@]}" --datasets n400,n800,n1600,n5402 \
    --conditions port_py3_epang_T8_catinfo,port_py3_epang_T8_K25_catinfo \
    --repeat 1 --untimed

"$PY3" scripts/05_concordance.py --details
"$PY3" scripts/06_cutoff_sensitivity.py --by-bin
"$PLOT_PY" scripts/03_plot.py
"$PLOT_PY" scripts/08_build_report.py
echo "=== alpha fix measured ($(date)) ==="
