#!/bin/sh
# table1.sh -- Table I of the paper: the five two-loop terms (Volkov's scheme) at h = 0.6 against their exact
# values, as fractional residuals x 1e12, and the residual of their sum against A1^(4).
#
#   sh scripts/table1.sh            # builds bin/sincqed with g2 if needed; ~15 s on 4 cores, any machine
#
# Settings: H [0.6], EPS [1e-20], WINDOW [72] (N = 2 ceil(WINDOW/h) + 1; with |z| <= 48 the crossed and
# ladder graphs stop near 1e-10 below h = 0.7; capped at N = 249, see below), OUT [table1_h$H.txt] for the raw results file.
# Values include the mirror weights.  The residual is (S - exact)/exact; for the sum, against
# A1^(4) = 197/144 + pi^2/12 - pi^2 ln2/2 + 3 zeta(3)/4 = -0.328478965579193.
# There are no timings in the table, so any machine will do (Linux 4 cores, 2026-11-04: Pi2 +1046.06,
# corner +88.66, self-energy -51.49, ladder -3.10, crossed -8.74, sum -130.27); another machine or
# thread count can move the last digit shown.
set -e
H=${H:-0.6}; EPS=${EPS:-1e-20}; WINDOW=${WINDOW:-72}; OUT=${OUT:-table1_h$H.txt}
ROOT=$(cd "$(dirname "$0")/.." && pwd); BIN=$ROOT/bin/sincqed
if [ ! -x "$BIN" ] || ! "$BIN" list 2>/dev/null | grep -q '^  g2_00 '; then
  echo "-- building bin/sincqed with g2 (set CC for macOS, e.g. CC=gcc-14)"
  ( cd "$ROOT" && make -j4 CLASSES="g2" ${CC:+CC=$CC} >/dev/null )
fi
# N <= 249: the driver stores lattice indices in 8 bits and rounds K = (N-1)/2 up to a multiple of 4 at
# top = 1, so K <= 124.  Below h ~ 0.58 this narrows the window (h = 0.5: |z| <= 62).
N=$(awk -v h=$H -v w=$WINDOW 'BEGIN{ n = w / h; m = (n == int(n)) ? n : int(n) + 1; if (m > 124) m = 124; print 2 * m + 1 }')
echo "-- sincqed run $H $N $EPS 1 all2"; rm -f "$OUT"     # the driver appends to an existing results file
SINCQED_RESULTS="$OUT" SINCQED_LOG=/dev/null "$BIN" run "$H" "$N" "$EPS" 1 all2 > /dev/null
awk -v h=$H -v N=$N -v eps=$EPS '
BEGIN { name["v2_00"] = "Pi2 insertion"; name["g2_00"] = "Corner"; name["g2_01"] = "Self-energy"
        name["g2_02"] = "Crossed"; name["g2_03"] = "Ladder"
        ord[1] = "v2_00"; ord[2] = "g2_00"; ord[3] = "g2_01"; ord[4] = "g2_03"; ord[5] = "g2_02"; A1 = -0.328478965579193 }
$1 in name { n = NF; S[$1] = $(n - 10); ref[$1] = $(n - 1); d[$1] = $n; tot += $(n - 10) }
END {
  printf "Table I: h = %s, N = %d, eps = %s\n\n", h, N, eps
  printf "%-14s %-8s %20s %20s %12s\n", "graph", "kernel", "Sinc sum", "exact", "1e12 Delta"
  for (i = 1; i <= 5; i++) { k = ord[i]; if (!(k in S)) { print "missing " k; exit 1 }
    printf "%-14s %-8s %+20.15f %+20.15f %+12.2f\n", name[k], k, S[k], ref[k], d[k] / ref[k] * 1e12 }
  printf "%-14s %-8s %+20.15f %+20.15f %+12.2f\n", "A1^(4) (sum)", "", tot, A1, (tot - A1) / A1 * 1e12
}' "$OUT"
