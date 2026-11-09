#!/bin/sh
# demo_scan.sh -- h and eps scan of the unreduced one- and two-loop Sinc sums (sincqed demo).
#
#   sh scripts/demo_scan.sh
#
# For each eps in EPSS the demo sum over all free Schwinger parameters (2 at one loop, 5 at two
# loops) is run at every h in HS; at each eps the error against the exact coefficients falls with
# h until it reaches that eps's truncation floor.  Results are appended to $OUT/demo_<eps>.txt (one
# line per graph and h, plus "# SUM" lines with the one- and two-loop totals); the summary
# $OUT/demo_sums.tsv has one row per (eps, h).  Memory grows as about h^-4: 1.1 GB at h = 0.6,
# eps = 1e-16.  Override any setting from the environment, e.g.
#   HS="1.0,0.8,0.6,0.5" EPSS="1e-14 1e-16" THREADS=24 sh scripts/demo_scan.sh
HS=${HS:-"1.5,1.25,1.0,0.9,0.8,0.7"}
EPSS=${EPSS:-"1e-8 1e-10 1e-12 1e-14 1e-16"}
OUT=${OUT:-demo}
BIN=${BIN:-$(cd "$(dirname "$0")/.." && pwd)/bin/sincqed}
[ -n "$THREADS" ] && export OMP_NUM_THREADS=$THREADS
mkdir -p "$OUT"
for e in $EPSS; do
    SINCQED_DEMO_RESULTS="$OUT/demo_$e.txt" SINCQED_LOG="$OUT/demo.log" "$BIN" demo "$HS" "$e" all1 all2 || exit 1
done
{ printf 'L\teps\th\tsum\tsum_minus_exact\n'
  for e in $EPSS; do
      awk -v e="$e" '$1 == "#" && $2 == "SUM" { sub("L=", "", $3); sub("h=", "", $4); printf "%s\t%s\t%s\t%s\t%s\n", $3, e, $4, $6, $8 }' "$OUT/demo_$e.txt"
  done; } > "$OUT/demo_sums.tsv"
echo "summary in $OUT/demo_sums.tsv"
