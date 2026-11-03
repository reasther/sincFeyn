#!/bin/sh
# run5_hscan.sh -- the five-loop runs behind Fig. 5: an h scan at matched eps plus one eps step, per graph.
#
#   sh scripts/run5batch.sh data/g5/graphs.txt 2.4 1e-6     # once: generate, build and run the h = 2.4 point
#   sh scripts/run5_hscan.sh                                 # the rest of the scan
#   python3 scripts/fit5.py runs5                            # h extrapolation + eps correction vs Volkov
#   python3 scripts/plot5_cost.py runs5                      # Fig. 5
#
# For every graph (default: data/g5/graphs.txt), top = 0 (the top-down sum) at
#   h = 2.4, 2.3, 2.2, 2.1 with eps = 1e-6 (h/2.4)^D, which keeps the pruned fraction about the same at
#       every h (so the exp(-c/h) fit in fit5.py sees only the h dependence), N = 2 ceil(48/h) + 1;
#   h = 2.4 with eps = 1e-5, for the eps correction delta = 0.24/0.76 (S(1e-6) - S(1e-5)).
# A run is skipped when OUT already holds a top-0 result for that graph and h at the target eps or below,
# so stopping and restarting is safe.  Settings: GRAPHS, HS ["2.4 2.3 2.2 2.1"], OUT [runs5], THREADS.
# Cost on a 24-thread Mac Studio (M2 Ultra): h = 2.4 ~4-10 min per graph, h = 2.1 ~20-80 min; the
# whole scan for eleven graphs ~15 h.  Memory below ~15 GB (D = 13 at h = 2.1).
# Run from the top of the tree after bin/sincqed has been built with the g5 kernels.
set -e
OUT=${OUT:-runs5}; HS=${HS:-"2.4 2.3 2.2 2.1"}
ROOT=$(cd "$(dirname "$0")/.." && pwd); cd "$ROOT"
GRAPHS=${GRAPHS:-$(grep -v '^#' data/g5/graphs.txt | grep .)}
mkdir -p "$OUT"; [ -n "$THREADS" ] && export OMP_NUM_THREADS=$THREADS
bin/sincqed list > "$OUT/.g5list" 2>/dev/null
dof(){ awk -v g="$1" '$1 == g { for(i = 2; i <= NF; i++) if($i ~ /^D=/){ sub(/D=/, "", $i); print $i; exit } }' "$OUT/.g5list"; }
done_at(){   # done_at NAME H EPS: a top-0 result at this h with eps <= EPS (1e-5 runs count only for 1e-5)?
    for f in "$OUT/$1_h$2_eps"*"_top0_N"*.txt; do
        [ -s "$f" ] && grep -q "^$1 " "$f" || continue
        e=$(echo "$f" | sed 's/.*_eps\(.*\)_top0.*/\1/')
        awk -v a="$e" -v b="$3" 'BEGIN{ exit !(a + 0 <= b * 1.0001 && a + 0 >= b * 0.1) }' && return 0
    done; return 1
}
run(){       # run NAME H N EPS
    TAG="$1_h$2_eps$4_top0_N$3"
    if done_at "$1" "$2" "$4"; then echo "-- $1 h$2 eps $4 done, skipped"; return; fi
    echo "-- $TAG  ($(date '+%H:%M'))"
    ( cd "$OUT" && SINCQED_TIMING=1 SINCQED_RESULTS="$TAG.txt" SINCQED_LOG="$TAG.log" \
        "$ROOT/bin/sincqed" run "$2" "$3" "$4" 0 "$1" < /dev/null ) | grep "^$1 \|resuming" || true
}
for h in $HS; do
  N=$(awk -v h=$h 'BEGIN{ m = int(48 / h); if(m * h < 47.99) m++; print 2 * m + 1 }')
  for D in 12 13; do
    for g in $GRAPHS; do
      d=$(dof $g); [ -z "$d" ] && { echo "-- $g not in this build, skipped"; continue; }
      [ "$d" = "$D" ] || continue
      eps=$(awk -v h=$h -v d=$d 'BEGIN{ printf "%.1e", 1e-6 * (h / 2.4)^d }' | sed 's/\.0e/e/; s/e-0/e-/')
      run $g $h $N $eps
    done
  done
done
for g in $GRAPHS; do [ -n "$(dof $g)" ] && run $g 2.4 41 1e-5; done
echo "== done: python3 scripts/fit5.py $OUT"
