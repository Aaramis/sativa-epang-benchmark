#!/usr/bin/env bash
# Shared paths for the SATIVA benchmark. Everything downstream sources this.
# Results and interpretation: RESULTS.md
#
# Nothing site specific lives here. Put your own paths in config.local.sh (untracked),
# or export the four SATIVA_BENCH_* variables before sourcing this file:
#
#   SATIVA_BENCH_PY3_ENV      conda env with python 3, ete3 and epa-ng 0.3.8
#   SATIVA_BENCH_PY2_ENV      bare python 2.7 env, for the 2017 SATIVA commit
#   SATIVA_BENCH_PLOT_PY      python with matplotlib, markdown and pillow
#   SATIVA_BENCH_SOURCE_MSA   alignment to draw the size gradient from, headers
#                             >id|kingdom=...|phylum=...|...|species=...

BENCH_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ -f "${BENCH_ROOT}/config.local.sh" ] && source "${BENCH_ROOT}/config.local.sh"

for _var in SATIVA_BENCH_PY3_ENV SATIVA_BENCH_PY2_ENV SATIVA_BENCH_PLOT_PY SATIVA_BENCH_SOURCE_MSA; do
    if [ -z "${!_var:-}" ]; then
        echo "config.sh: $_var is not set. See the header of this file." >&2
    fi
done
unset _var

# The three SATIVA under test.
# sativa_epang/ is the modified SATIVA this report describes; see its PROVENANCE.md for the
# diff against upstream, the environment variables and the checksums.
SATIVA_PORT="${BENCH_ROOT}/sativa_epang/sativa.py"                          # EPA-ng version
SATIVA_UPSTREAM_PY3="${BENCH_ROOT}/env/sativa_upstream/sativa.py"           # v0.9.3, commit 68284c2
SATIVA_UPSTREAM_PY2="${BENCH_ROOT}/env/sativa_upstream_py2/sativa.py"       # last python-2 commit, 8c31962

PY3_ENV="${SATIVA_BENCH_PY3_ENV}"
PY3="${PY3_ENV}/bin/python"
PY2_ENV="${SATIVA_BENCH_PY2_ENV}"
PY2="${PY2_ENV}/bin/python"
PLOT_PY="${SATIVA_BENCH_PLOT_PY}"
SOURCE_MSA="${SATIVA_BENCH_SOURCE_MSA}"

# SATIVA parameters, identical in every condition
SATIVA_MODE="ultrafast"
SATIVA_CONF="0.4"
SATIVA_TAXCODE="BOT"
SATIVA_SEED="42"

export BENCH_ROOT SATIVA_PORT SATIVA_UPSTREAM_PY2 SATIVA_UPSTREAM_PY3
export PY3_ENV PY3 PY2_ENV PY2 PLOT_PY SOURCE_MSA
export SATIVA_MODE SATIVA_CONF SATIVA_TAXCODE SATIVA_SEED
