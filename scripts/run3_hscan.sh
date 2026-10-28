#!/bin/sh
# run3_hscan.sh -- every three-loop graph (all3: 41 entries = all 72 diagrams) on a ladder of h
# values, coarsest first, each h with its own eps.  Optionally (BUDGET) it stops before an h level
# predicted to overrun a time budget.
#
#   sh scripts/run3_hscan.sh [OUTDIR]            (default OUTDIR = runs3)
#   sh scripts/sum3.sh [OUTDIR]                  table of totals per h (also run after every level)
#
# Environment (all optional):
#   BIN=bin/sincqed       binary; on the GPU machine e.g. BIN=bin/sincqed-gpu
#   PAR=1                 graphs run at the same time (each one is a separate process)
#   BUDGET=0              hours; 0 = no limit (stop it by hand).  Otherwise no new h level is started
#                         if it is predicted to end after BUDGET hours
#   LADDER="h:eps ..."    the h values and their eps (default below)
#   ZMAX=48               window |z| <= ZMAX at every h (N is chosen from it)
#   GRAPHS="..."          graph names (default: all of all3, most expensive first)
#   GROWTH=3              assumed time ratio between successive levels until two have been timed
#   SINCQED_PREDICT=0     plain grid, as in all earlier three-loop runs
#
# Output: OUTDIR/h<h>/ holds res_<graph>.txt, log_<graph>.log, out_<graph>.txt per graph,
# results.txt and sincqed.log (all graphs of that h); OUTDIR/levels.txt holds one line per
# finished h (h eps N seconds).  Rerunning the same command resumes: graphs with a result are
# skipped, and a graph that was killed restarts from its checkpoint (in OUTDIR/h<h>/).
set -u
OUT=${1:-runs3}
BIN=${BIN:-bin/sincqed}
PAR=${PAR:-1}
BUDGET=${BUDGET:-0}
ZMAX=${ZMAX:-48}
GROWTH=${GROWTH:-3}
LADDER=${LADDER:-"1.5:1e-10 1.4:1e-10 1.3:1e-11 1.2:1e-11 1.1:1e-12 1.0:1e-12 0.9:1e-13 0.8:1e-14 0.7:1e-15"}
GRAPHS=${GRAPHS:-"g3_lbl g3_10 g3_20 g3_00 g3_05 g3_14 g3_22 g3_24 g3_26 g3_13 g3_02 g3_03 g3_15 g3_18 g3_23 g3_27 g3_12
 g3_17 g3_21 g3_19 g3_01 g3_06 g3_11 g3_04 g3_08 g3_16 g3_25 g3_07 g3_09
 v3_p4Vc v3_p4Vm v3_01 v3_05 v3_02 v3_06 v3_07 v3_03 v3_04 v3_p4Sm v3_00 v3_p4Sc"}
export SINCQED_PREDICT=${SINCQED_PREDICT:-0}
[ -x "$BIN" ] || { echo "no binary $BIN (build first, or set BIN)"; exit 1; }
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"
T0=$(date +%s); last=0; before=0
for he in $LADDER; do
    h=${he%%:*}; eps=${he#*:}
    # N = 2K+1 with K*h >= ZMAX and K a multiple of 4 (top = 1)
    N=$(awk -v h="$h" -v z="$ZMAX" 'BEGIN{k = int(z/h); if(k*h < z - 1e-9) k++; k = int((k+3)/4)*4; print 2*k+1}')
    D=$OUT/h$h
    now=$(date +%s); used=$((now - T0))
    if [ "$last" -gt 0 ] && [ "$BUDGET" != 0 ]; then
        est=$(awk -v a="$last" -v b="$before" -v g="$GROWTH" 'BEGIN{r = (b > 0) ? a/b : g; if(r < 1.5) r = 1.5; printf "%d", a*r}')
        if [ $((used + est)) -gt $(awk -v b="$BUDGET" 'BEGIN{printf "%d", b*3600}') ]; then
            echo "# stop before h=$h: predicted $((est/60)) min would end after the ${BUDGET} h budget (used $((used/60)) min)"
            break
        fi
        echo "# h=$h eps=$eps N=$N: predicted $((est/60)) min (used $((used/60)) min of ${BUDGET} h)"
    else
        echo "# h=$h eps=$eps N=$N"
    fi
    mkdir -p "$D"
    todo=""
    for g in $GRAPHS; do grep -q "^$g " "$D/res_$g.txt" 2>/dev/null || todo="$todo $g"; done
    t1=$(date +%s)
    if [ -n "$todo" ]; then
        printf '%s\n' $todo | xargs -P "$PAR" -I{} sh -c \
            'SINCQED_RESULTS="$1/res_$2.txt" SINCQED_LOG="$1/log_$2.log" SINCQED_CKPT_DIR="$1" "$3" run "$4" "$5" "$6" 1 "$2" > "$1/out_$2.txt" 2>&1 || echo "# $2 failed at h=$4 (see $1/out_$2.txt)"' \
            sh "$D" {} "$BIN" "$h" "$N" "$eps"
    fi
    t2=$(date +%s)
    for g in $GRAPHS; do grep -v '^#' "$D/res_$g.txt" 2>/dev/null; done > "$D/results.txt"
    for g in $GRAPHS; do cat "$D/log_$g.log" 2>/dev/null; done > "$D/sincqed.log"
    n=$(grep -c . "$D/results.txt")
    if [ -n "$todo" ] || ! grep -q "^$h " "$OUT/levels.txt" 2>/dev/null; then   # seconds add up over resumed runs
        p=$(awk -v h="$h" '$1 == h {s = $4} END {print s + 0}' "$OUT/levels.txt" 2>/dev/null || echo 0)
        echo "$h $eps $N $((p + t2 - t1)) $n" >> "$OUT/levels.txt"
    fi
    sh "$HERE/sum3.sh" "$OUT" | tail -1
    if [ -n "$todo" ]; then before=$last; last=$((t2 - t1)); fi
done
sh "$HERE/sum3.sh" "$OUT"
