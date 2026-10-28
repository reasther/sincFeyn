#!/bin/sh
# scan2.sh -- two-loop convergence scan: h at a fixed, tiny eps.
#
#   sh scripts/scan2.sh            run (or resume; kill and restart at any time)
#   sh scripts/scan2.sh summary    rebuild the tables from what has finished
#
# All 7 two-loop diagrams (5 entries: g2_00..g2_03 and the VP insertion v2_00) on the plain grid
# (no rotations, unit scales).
# eps = 1e-15 leaves a pruning floor of ~2e-10 on g2_02/g2_03; use EPSS=1e-18 to push it to ~1e-11.
# Box |z| <= 40, so h is the only discretisation parameter.
# Results in conv2/: totals.tsv (sum against exact A2), compare.tsv , summary.tsv.
# A few minutes on 8 cores (h = 0.7 already gives the sum to ~1e-9).  Any setting below can be overridden from the
# environment, e.g.  HS2="1.0 0.8" EPSS=1e-18 sh scripts/scan2.sh
export GRAPHS=${GRAPHS:-all2}
export HS2=${HS2:-"1.5 1.3 1.2 1.0 0.9 0.8 0.7 0.6 0.5"}
export EPSS=${EPSS:-1e-15}
export GRIDS=${GRIDS:-plain}
export THREADS=${THREADS:-8}
export OUT=${OUT:-conv2}
exec sh "$(dirname "$0")/convergence.sh" "$@"
