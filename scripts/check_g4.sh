#!/bin/sh
# Checks before a four-loop group V run (all 269 g4 graphs, 518 diagrams).
#   sh scripts/check_g4.sh                  kernels only: kernels-ext/g4 (KD=DIR to override) holds the 269 FORM kernels, numbered as
#                                           the generator numbers them (scripts/g4_kernels.tsv), nothing else
#   sh scripts/check_g4.sh bin/sincqed-gpu  also runs four quick graphs (D = 7, 8, 9, 10) at h=1.5 N=57
#                                           eps=1e-10 with the predictor and compares them with the Mac Studio
#                                           group V run (same settings, FORM kernels): minutes on an A100
# Each kernel's Volkov label, dimension and weight must match the table and its header must say
# "Dirac algebra: FORM".  The kernels may be the 2026-10-10 FORM release (sincqed_2026-10-10_form.tar.xz)
# or the 2026-10-21 ones (kernels-g4_2026-10-21.tar.xz: same-momentum series classes and the
# mirror-invariant Cheng-Wu photon, 67 graphs differ, 65 of them with a smaller D), or the 2026-10-23
# ones (kernels-g4_2026-10-23.tar.xz or generator/regen.sh: the 2026-10-21 kernels with the one-point
# u rule and the packed-point wrapper, every file differs, same D as 2026-10-21), or the 2026-10-25
# ones (kernels-g4_2026-10-25.tar.xz: line identifications, every file differs, same D), or the v0.9.0
# release (kernels-g4.tar.xz: generated with GEN_HBASIS=2, the generator default since 2026-10-31); a
# mix of the two is reported.  A kernel whose sha256 differs from both is reported but passes: kernels
# regenerated with FORM differ in the order of subexpressions only.  Exit status 0 if everything passes.
set -u
TAB=$(dirname "$0")/g4_kernels.tsv
KD=${KD:-$( [ -d kernels-ext/g4 ] || [ ! -d kernels/g4 ] && echo kernels-ext/g4 || echo kernels/g4 )}
BIN=${1:-}
if command -v sha256sum >/dev/null 2>&1; then SHA="sha256sum"; else SHA="shasum -a 256"; fi
fail=0; nid=0; nre=0; n=0; nold=0; nnew=0; nh10=0; nh21=0; nh23=0; nh25=0; nh11=0
[ -d "$KD" ] || { echo "FAIL  no directory $KD (unpack the g4 release archive into kernels-ext/)"; exit 1; }
# anything in $KD that is not one of the 269 (e.g. Finder copies "g4_00 2.c", or old Mac-numbered files)
for f in "$KD"/*; do
    b=$(basename "$f"); g=${b%.c}
    if [ "$b" = "$g" ] || ! grep -q "^$g	" "$TAB"; then echo "FAIL  unexpected file $f"; fail=$((fail+1)); fi
done
while IFS='	' read -r g label D w sha Dm sham sh23 sh25 sh11; do
    case $g in \#*) continue;; esac
    n=$((n+1)); f="$KD/$g.c"
    if [ ! -f "$f" ]; then echo "FAIL  $g: missing $f"; fail=$((fail+1)); continue; fi
    line=$(grep -h "^const graph_t G_$g=" "$f")
    want="{\"$g\",\"$label\",4,$D,$w,"; wantm="{\"$g\",\"$label\",4,$Dm,$w,"
    case $line in
        *"$want"*) [ "$D" = "$Dm" ] || nold=$((nold+1));;
        *"$wantm"*) nnew=$((nnew+1));;
        *) echo "FAIL  $g: expected '$label' D=$D (or $Dm in the 2026-10-21 kernels) weight $w, file has: $(echo "$line" | cut -c1-90)"; fail=$((fail+1)); continue;;
    esac
    head -3 "$f" | grep -q "Dirac algebra: FORM" || { echo "FAIL  $g: not a FORM kernel (header lacks 'Dirac algebra: FORM')"; fail=$((fail+1)); }
    h=$($SHA "$f" | cut -d' ' -f1)
    if [ "$h" = "$sha" ]; then nid=$((nid+1)); [ "$sha" = "$sham" ] || nh10=$((nh10+1))
    elif [ "$h" = "$sham" ]; then nid=$((nid+1)); nh21=$((nh21+1))
    elif [ "$h" = "$sh23" ]; then nid=$((nid+1)); nh23=$((nh23+1)); [ "$sha" = "$sham" ] || nh21=$((nh21+1))
    elif [ "$h" = "$sh25" ]; then nid=$((nid+1)); nh25=$((nh25+1)); [ "$sha" = "$sham" ] || nh21=$((nh21+1))
    elif [ "$h" = "$sh11" ]; then nid=$((nid+1)); nh11=$((nh11+1)); [ "$sha" = "$sham" ] || nh21=$((nh21+1)); else nre=$((nre+1)); fi
done < "$TAB"
echo "kernels: $n checked, $fail failed; of the files present, $nid identical to a FORM release ($nh23 of them the 2026-10-23 one-point packed kernels, $nh25 the 2026-10-25 line-identification ones, $nh11 the v0.9.0 release), $nre differ (labels match)"
echo "of the 67 graphs that changed on 2026-10-21: $nh21 are the new kernels, $nh10 the 2026-10-10 ones (by sha256; D: $nnew new, $nold old)"
[ $nh21 -eq 0 ] || [ $nh10 -eq 0 ] || echo "NOTE  $KD mixes 2026-10-10 and 2026-10-21 kernels (both are exact)"
[ $fail -eq 0 ] || { echo "$KD is NOT the FORM group V set: fix it before running"; exit 1; }
[ -n "$BIN" ] || exit 0

# quick graphs against the Mac Studio values (1e-7 relative: x86 and Apple silicon round differently,
# so a point near the pruning threshold can go either way; g4_184 differs by 5.5e-8) (h=1.5 N=57 eps=1e-10, SINCQED_PREDICT=1, SINCQED_DDSLACK=100)
[ -x "$BIN" ] || { echo "FAIL  no binary $BIN"; exit 1; }
OUT=${TMPDIR:-/tmp}/check_g4.$$; mkdir -p "$OUT"
export SINCQED_PREDICT=1 SINCQED_DDSLACK=100 SINCQED_CKPT=0 SINCQED_PROGRESS=0 SINCQED_LOG="$OUT/check.log"
SINCQED_RESULTS="$OUT/check.txt" "$BIN" run 1.5 57 1e-10 g4_73 g4_72 g4_66 g4_184 > "$OUT/out.txt" 2>&1
rfail=0
while read -r g ref ev; do
    v=$(awk -v g="$g" '$1 == g {print $(NF-10), $(NF-3)}' "$OUT/check.txt")
    set -- $v
    if [ $# -lt 2 ]; then echo "FAIL  $g: no result (see $OUT/out.txt)"; rfail=$((rfail+1)); continue; fi
    ok=$(awk -v v="$1" -v r="$ref" 'BEGIN{d=(v-r)/r; if(d<0)d=-d; print (d<=1e-7) ? "PASS" : "FAIL"}')
    [ "$ok" = PASS ] || rfail=$((rfail+1))
    printf "%-5s %-7s %+.15e  Mac %+.15e  evals %s (Mac %s)\n" "$ok" "$g" "$1" "$ref" "$2" "$ev"
done <<'REF'
g4_73 -5.732451123610411e-01 48754542
g4_72 -5.601438731696874e-01 307212358
g4_66 -1.303697451634703e+00 1187913294
g4_184 -1.566357003768172e+00 3340183117
REF
echo "runs: $((4-rfail))/4 match the Mac Studio values to 1e-7 (outputs in $OUT)"
[ $rfail -eq 0 ]
