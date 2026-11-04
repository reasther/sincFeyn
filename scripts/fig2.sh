#!/bin/sh
# fig2.sh -- the data for Fig. 2 of the paper: fractional residuals of the five two-loop terms and of A1^(4)
# against h, as one table for plotting (Mathematica: Import["fig2.tsv", "TSV"], skipping "#" lines).
#
#   sh scripts/fig2.sh [OUTDIR]          (default fig2/; builds bin/sincqed with g2 if needed)
#
# h runs from HMAX [1.60] down to HMIN [0.55] in steps of DH [0.01], eps = EPS [1e-20], top = 1,
# N = 2 ceil(72/h) + 1 capped at 249 (8-bit lattice indices; below h ~ 0.58 the window narrows to
# |z| <= 124 h).  Each h writes OUTDIR/two_loop_h<h>.txt (skipped if present, so restarts are cheap)
# and the table OUTDIR/fig2.tsv has one row per h:
#   h  N  Delta(Pi2)  Delta(corner)  Delta(self-energy)  Delta(ladder)  Delta(crossed)  Delta(A1^(4))
# with Delta = (S - exact)/exact (signed; mirror weights included).  About 6 minutes on 4 cores.
set -e
OUT=${1:-fig2}; HMAX=${HMAX:-1.60}; HMIN=${HMIN:-0.55}; DH=${DH:-0.01}; EPS=${EPS:-1e-20}
ROOT=$(cd "$(dirname "$0")/.." && pwd); BIN=$ROOT/bin/sincqed
if [ ! -x "$BIN" ] || ! "$BIN" list 2>/dev/null | grep -q '^  g2_00 '; then
  echo "-- building bin/sincqed with g2 (set CC for macOS, e.g. CC=gcc-14)"
  ( cd "$ROOT" && make -j4 CLASSES="g2" ${CC:+CC=$CC} >/dev/null )
fi
mkdir -p "$OUT"
for h in $(awk -v a=$HMAX -v b=$HMIN -v d=$DH 'BEGIN{ n = int((a - b)/d + 0.5); for(i = 0; i <= n; i++) printf "%.2f\n", a - i*d }'); do
  f="$OUT/two_loop_h$h.txt"
  [ -s "$f" ] && grep -q '^# SUM' "$f" && continue
  N=$(awk -v h=$h 'BEGIN{ n = 72/h; m = (n == int(n)) ? n : int(n) + 1; if (m > 124) m = 124; print 2*m + 1 }')
  echo "-- h $h N $N"; rm -f "$f"
  SINCQED_RESULTS="$f" SINCQED_LOG=/dev/null "$BIN" run "$h" "$N" "$EPS" 1 all2 > /dev/null
done
{ echo "# h	N	Pi2	corner	selfenergy	ladder	crossed	A1sum   (Delta = (S - exact)/exact; sincqed run h N $EPS 1 all2)"
  for f in "$OUT"/two_loop_h*.txt; do
    awk 'BEGIN{ A1 = -0.328478965579193 }
      $1 ~ /^(g2|v2)_/ { n = NF; d[$1] = $n / $(n-1); tot += $(n-10); h = $(n-14); N = $(n-13) }
      END { printf "%s\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", h, N, d["v2_00"], d["g2_00"], d["g2_01"], d["g2_03"], d["g2_02"], (tot - A1)/A1 }' "$f"
  done | sort -g -r; } > "$OUT/fig2.tsv"
echo "== $OUT/fig2.tsv ($(grep -vc '^#' "$OUT/fig2.tsv") h values)"
