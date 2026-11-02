#!/bin/sh
# sum3.sh [OUTDIR] -- totals per h from a run3_hscan.sh output directory.
# "shadow" is the sum over graphs of the pruning loss left at eps (from the "eps extrapolation" log
# lines); "corrected" adds it to the total.  diff columns are against the exact
# A3 = 1.18124145658720 (Laporta-Remiddi), shown only when all 41 graphs are present.
# Also writes OUTDIR/pergraph.tsv: graph h eps value shadow value+shadow time_s, for every graph and h.
OUT=${1:-runs3}
[ -f "$OUT/levels.txt" ] || { echo "no $OUT/levels.txt"; exit 1; }
printf '%-5s %-7s %-4s %-3s %-17s %-10s %-17s %-10s %-10s %-8s\n' h eps N n total shadow corrected diff diff_corr wall_min
: > "$OUT/pergraph.tsv.tmp"
awk '{l[$1] = $0} END{for(h in l) print l[h]}' "$OUT/levels.txt" | sort -k1,1gr | while read h eps N secs n; do
    D=$OUT/h$h
    awk -v h="$h" -v eps="$eps" -v N="$N" -v secs="$secs" -v pg="$OUT/pergraph.tsv.tmp" '
        FILENAME == ARGV[1] && / eps (shadow|extrapolation): / { for(i = 1; i <= NF; i++) if($i == "~") sh[$2] = $(i+1) }
        FILENAME == ARGV[2] && !/^#/ {
            for(i = 4; i <= NF; i++) if($i ~ /^[+-][0-9]\.[0-9]+e[+-][0-9]+$/){ v = $i; break }
            n++; tot += v; s += sh[$1]
            printf "%s\t%s\t%s\t%.15e\t%.3e\t%.15e\t%s\n", $1, h, eps, v, sh[$1], v + sh[$1], $(NF-2) >> pg
        }
        END {
            A3 = 1.18124145658720
            d = (n == 41) ? sprintf("%+.2e", tot - A3) : "-"; dc = (n == 41) ? sprintf("%+.2e", tot + s - A3) : "-"
            printf "%-5s %-7s %-4s %-3d %+.14f %+.3e %+.14f %-10s %-10s %.1f\n", h, eps, N, n, tot, s, tot + s, d, dc, secs/60
        }' "$D/sincqed.log" "$D/results.txt"
done
sort -k1,1 -k2,2gr "$OUT/pergraph.tsv.tmp" > "$OUT/pergraph.tsv"; rm -f "$OUT/pergraph.tsv.tmp"
