#!/bin/sh
# demo_converge.sh -- eps-converged unreduced one- and two-loop Sinc sums (sincqed demo) on a ladder of h.
#
#   sh scripts/demo_converge.sh
#
# Runs all1 + all2 (g1_00, g2_00..03, v2_00) at every h from HMAX down to HMIN in steps of HSTEP, each
# with eps = 1e-18, lowered at small h so the eps truncation stays negligible next to the h error:
#   eps(h) = min(1e-18, 1e-15 exp(-18 (1/h - 1)))      (1e-18 at h >= 0.76, 4.5e-19 at 0.7, 6e-21 at 0.6)
# The two-loop h error falls about as 1e-6 exp(-18 (1/h - 1)) and the truncation is up to ~1e6 eps
# (at h = 0.6, eps 1e-18 vs 1e-20 moves single graphs by up to 9e-13, against a crossed-graph h error
# of 5e-13), so a fixed 1e-18 is not enough below h ~ 0.75.
# Output in $OUT: demo_conv.txt (one line per graph and h, plus "# SUM" lines), graphs.tsv (value and
# value - exact per graph and h), sums.tsv (one row per h: one- and two-loop totals and their errors),
# demo.log.  Restartable: h values already in
# demo_conv.txt are skipped.  Settings (environment): HMAX 1.6, HMIN 0.6, HSTEP 0.02, THREADS, OUT
# (demo_conv), BIN; EPSFAC (1e-15) and EPSMAX (1e-18) set the eps rule
# (EPSFAC=1 EPSMAX=1e-18 gives a fixed 1e-18); GRAPHS ("all1 all2") picks the graphs, e.g. a fine scan of
# the corner through its sign change:
#   GRAPHS=g2_00 HMAX=1.3 HMIN=1.25 HSTEP=0.001 OUT=corner_fine sh scripts/demo_converge.sh
HMAX=${HMAX:-1.6}; HMIN=${HMIN:-0.6}; HSTEP=${HSTEP:-0.02}
EPSFAC=${EPSFAC:-1e-15}; EPSMAX=${EPSMAX:-1e-18}
OUT=${OUT:-demo_conv}
GRAPHS=${GRAPHS:-"all1 all2"}
BIN=${BIN:-$(cd "$(dirname "$0")/.." && pwd)/bin/sincqed}
[ -n "$THREADS" ] && export OMP_NUM_THREADS=$THREADS
export G2_ABSERR=${G2_ABSERR:-1e-18}
[ -x "$BIN" ] || { echo "no executable $BIN (run make)" >&2; exit 1; }
mkdir -p "$OUT"
R="$OUT/demo_conv.txt"
h=$HMAX
while awk -v h="$h" -v m="$HMIN" 'BEGIN { exit !(h >= m - 1e-9) }'; do
    if ! grep -q "^# SUM L=2 h=$(printf '%.4f' "$h") " "$R" 2>/dev/null; then
        eps=$(awk -v h="$h" -v f="$EPSFAC" -v m="$EPSMAX" 'BEGIN { e = f*exp(-18*(1/h - 1)); if (e > m) e = m; printf "%.1e", e }')
        echo "h=$h eps=$eps" >&2
        SINCQED_DEMO_RESULTS="$R" SINCQED_LOG="$OUT/demo.log" "$BIN" demo "$h" "$eps" $GRAPHS > /dev/null || exit 1
    fi
    h=$(awk -v h="$h" -v s="$HSTEP" 'BEGIN { printf "%.4f", h - s }')
done
{ printf 'h\teps\tone_loop\tone_loop_minus_exact\ttwo_loop\ttwo_loop_minus_exact\n'
  awk '$1 == "#" && $2 == "SUM" { L = substr($3, 3); h = substr($4, 3); e = substr($5, 5)
                                   if (L == 1) { v1[h] = $6; d1[h] = $8 } else { v2[h] = $6; d2[h] = $8; ee[h] = e } }
       END { for (h in v2) printf "%s\t%s\t%s\t%s\t%s\t%s\n", h, ee[h], v1[h], d1[h], v2[h], d2[h] }' "$R" | sort -gr; } > "$OUT/sums.tsv"
# per graph, against the exact values (Volkov's scheme, mirror graphs included)
awk 'BEGIN { x["g1_00"] = 0.5; x["g2_00"] = -0.064020941344682815593; x["g2_01"] = -0.589978022282742145491
             x["g2_02"] = -0.467645446093613651600; x["g2_03"] = 0.777478022282742145491; x["v2_00"] = 0.015687421859102682611
             print "graph\th\teps\tvalue\tminus_exact" }
     $1 in x { printf "%s\t%s\t%s\t%s\t%+.4e\n", $1, $5, $6, $7, $7 - x[$1] }' "$R" > "$OUT/graphs.tsv"
echo "results in $R; per graph in $OUT/graphs.tsv; totals (all1 all2 only) in $OUT/sums.tsv" >&2
