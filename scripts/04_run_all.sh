#!/usr/bin/env bash
# Every measurement in RESULTS.md, in the order it is reported. About 3 h on 32 cores.
#
# Phase 1 and 2 are timed, so they run one job at a time and alone on the machine. Phases 3
# to 5 answer agreement questions, where wall clock does not matter: they pass --untimed,
# which records timing_trusted=false so no figure can pick those seconds up by accident.
set -euo pipefail
cd "$(dirname "$0")/.."
source config.sh

TIMEOUT="${TIMEOUT:-10800}"
DATASETS_JSON="${DATASETS_JSON:-data_ascii/datasets.json}"
RUNS_DIR="${RUNS_DIR:-results/runs}"
BENCH=(--datasets-json "$DATASETS_JSON" --runs-dir "$RUNS_DIR" --timeout "$TIMEOUT")

echo "=== phase 1: small sizes, every condition, 3 replicates ($(date)) ==="
"$PY3" scripts/02_bench.py "${BENCH[@]}" --datasets n100,n200,n400 \
    --conditions upstream_py2_raxml_T1,upstream_py2_raxml_T8,upstream_py3_raxml_T1,upstream_py3_raxml_T8,port_py3_raxml_T1,port_py3_raxml_T8,port_py3_epang_T1,port_py3_epang_T2,port_py3_epang_T4,port_py3_epang_T8,port_py3_epang_T16 \
    --repeat 3

echo "=== phase 2: size gradient, timed ($(date)) ==="
"$PY3" scripts/02_bench.py "${BENCH[@]}" --datasets n800,n1600,n5402 \
    --conditions upstream_py3_raxml_T1,upstream_py3_raxml_T8,port_py3_epang_T1,port_py3_epang_T2,port_py3_epang_T4,port_py3_epang_T8,port_py3_epang_T16,port_py3_epang_T1_K25,port_py3_epang_T8_K25 \
    --repeat 1

echo "=== phase 3: fold sweep ($(date)) ==="
"$PY3" scripts/02_bench.py "${BENCH[@]}" --datasets n100,n200,n400 \
    --conditions port_py3_epang_T8_K10,port_py3_epang_T8_K25,port_py3_epang_T8_K50,port_py3_epang_T8_Kexact \
    --repeat 1 --untimed

echo "=== phase 4: which placement rule matters ($(date)) ==="
"$PY3" scripts/02_bench.py "${BENCH[@]}" --datasets n200,n400 \
    --conditions port_py3_epang_T8_noheur,port_py3_epang_T8_blo,port_py3_epang_T8_acc999,port_py3_epang_T8_raxmlmatch \
    --repeat 1 --untimed

echo "=== phase 5: fitted alpha, and the reference pinned to GTRGAMMA ($(date)) ==="
"$PY3" scripts/02_bench.py "${BENCH[@]}" --datasets n800,n1600,n5402 \
    --conditions port_py3_epang_T8_gammainfo,port_py3_epang_T8_K25_gammainfo \
    --repeat 1 --untimed
"$PY3" scripts/02_bench.py "${BENCH[@]}" --datasets n800,n1600 \
    --conditions upstream_py3_raxml_T8_gamma --repeat 1 --untimed

echo "=== tables and figure ($(date)) ==="
"$PY3" scripts/05_concordance.py --runs-dir "$RUNS_DIR" --details
"$PY3" scripts/06_cutoff_sensitivity.py --runs-dir "$RUNS_DIR"
"$PLOT_PY" scripts/03_plot.py --runs-dir "$RUNS_DIR"
echo "=== done ($(date)) ==="
