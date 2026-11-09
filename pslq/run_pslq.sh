#!/bin/sh
# run_pslq.sh -- the high-precision Sinc runs behind the PSLQ closed forms of the two-loop graphs
# (Volkov's subtraction scheme, mirror graphs included; see README.md).
#
#   sh run_pslq.sh                 the runs used for the paper (about 23 min on 4 cores here, a few GB of memory)
#   QUICK=1 sh run_pslq.sh         the same chain at coarse steps (a few minutes): checks the build
#   THREADS=24 sh run_pslq.sh      threads (default: all)
#   python3 pslq/identify.py pslq_runs/values.tsv     optional PSLQ step (mpmath)
#
# 1. demo/bin/sincqed demo: the unreduced sums over all five free Schwinger parameters, double precision
#    with a double-double fallback (G2_ABSERR=1e-18), at h = 0.5, eps = 1e-20.  Gives the ladder (g2_03)
#    and corner (g2_00) values that PSLQ was run on (~14 digits), plus the crossed, self-energy and Pi2
#    graphs, whose sum is checked against A_1^{(4)}.
# 2. quad/se_qd: the self-energy (g2_01), reduced 3-D sum in quad-double (~64 digits), h = 0.4,
#    eps = 1e-26 (~17 digits; h = 0.3 gives ~21).
# Settings (environment): H_DEMO, EPS_DEMO, H_SE, EPS_SE, OUT (pslq_runs).  Restartable: a step whose
# output file is complete is skipped (delete it to rerun).
DIR=$(cd "$(dirname "$0")" && pwd)
OUT=${OUT:-pslq_runs}
if [ -n "$QUICK" ]; then
    H_DEMO=${H_DEMO:-0.7}; EPS_DEMO=${EPS_DEMO:-1e-18}; H_SE=${H_SE:-0.7}; EPS_SE=${EPS_SE:-1e-22}
else
    H_DEMO=${H_DEMO:-0.5}; EPS_DEMO=${EPS_DEMO:-1e-20}; H_SE=${H_SE:-0.4}; EPS_SE=${EPS_SE:-1e-26}
fi
[ -n "$THREADS" ] && export OMP_NUM_THREADS=$THREADS
SQ="$DIR/demo/bin/sincqed"; SE="$DIR/quad/se_qd"
[ -x "$SQ" ] && [ -x "$SE" ] || { echo "build first: make" >&2; exit 1; }
mkdir -p "$OUT"
T0=$(date +%s)

# 1. unreduced double / double-double sums
D="$OUT/demo_h${H_DEMO}_eps${EPS_DEMO}.txt"
if ! grep -q '^# SUM L=2' "$D" 2>/dev/null; then
    rm -f "$D"
    echo "[$(( $(date +%s) - T0 ))s] demo h=$H_DEMO eps=$EPS_DEMO (all two-loop graphs)" >&2
    G2_ABSERR=1e-18 SINCQED_DEMO_RESULTS="$D" SINCQED_LOG="$OUT/demo.log" "$SQ" demo "$H_DEMO" "$EPS_DEMO" all2 > "$OUT/demo.out" || exit 1
fi
# 2. self-energy in quad-double
S="$OUT/se_h${H_SE}_eps${EPS_SE}.txt"
if ! grep -q 'value =' "$S" 2>/dev/null; then
    echo "[$(( $(date +%s) - T0 ))s] se_qd h=$H_SE eps=$EPS_SE" >&2
    "$SE" "$H_SE" "$EPS_SE" > "$S.tmp" && mv "$S.tmp" "$S" || exit 1
fi
echo "[$(( $(date +%s) - T0 ))s] done" >&2

# values.tsv: graph, value, estimated error (the h error: the measured convergence of each sum,
# err ~ 1e-6 exp(-18 (1/h - 1)) for the unreduced sums and ~3e-11 exp(-17.6 (1/h - 5/3)) for the
# reduced self-energy), source.  This is what pslq/identify.py reads.
{
printf 'graph\tvalue\terr\tsource\n'
tol=$(awk -v h="$H_DEMO" 'BEGIN { printf "%.1e", 1e-6*exp(-18*(1/h - 1)) }')
awk -v t="$tol" -v src="demo h=$H_DEMO eps=$EPS_DEMO" '
    $1 == "g2_00" { n = "corner" } $1 == "g2_01" { n = "self-energy_demo" } $1 == "g2_02" { n = "crossed" }
    $1 == "g2_03" { n = "ladder" } $1 == "v2_00" { n = "Pi2" }
    $1 ~ /^(g2_0[0-3]|v2_00)$/ { printf "%s\t%s\t%s\t%s\n", n, $7, t, src }
    $1 == "#" && $2 == "SUM" && $3 == "L=2" { printf "sum\t%s\t%s\t%s\n", $6, t, src }' "$D"
tol=$(awk -v h="$H_SE" 'BEGIN { printf "%.1e", 3e-11*exp(-17.6*(1/h - 5/3)) }')
# the quad-double value as its four exact double components (summed exactly by identify.py)
awk -v t="$tol" -v src="se_qd h=$H_SE eps=$EPS_SE" '{ for (i = 1; i <= NF; i++) if ($i == "[parts") { printf "self-energy\t%s,%s,%s,%s\t%s\t%s\n", $(i+1), $(i+2), $(i+3), substr($(i+4), 1, length($(i+4)) - 1), t, src } }' "$S"
} > "$OUT/values.tsv"
cat "$OUT/values.tsv"
