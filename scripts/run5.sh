#!/bin/sh
# run5.sh -- one five-loop graph without lepton loops, from nothing to a value compared with Volkov.
#
#   sh scripts/run5.sh "6; 1-8, 2-9, 3-10, 4-11, 5-7" [H] [EPS]
#   sh scripts/run5.sh g5_3086 2.4 1e-8
#
# The graph is a Volkov string (either orientation; "6;1*8;2*9;..." as in his supplement also works)
# or a kernel name.  Steps, each skipped when its output already exists:
#   1. generate the kernel (gen.py, FORM O1: O3 does not finish five-loop nodes) into $KEXT/g5/
#   2. vectorise (vecg.py), optimise (optimize.sh) and split the FORM blocks (splitform.py), without
#      which gcc needs more than 14 GB per five-loop kernel
#   3. build bin/sincqed with $CLASSES (every kernel in $KEXT/g5 is linked)
#   4. run at step H [2.4] and threshold EPS [1e-7], N = 2*ceil(WINDOW/h)+1 with WINDOW [48]
#   5. print the value, the driver's prune estimate (and its eps-shadow estimate) and the difference from
#      Volkov's 2019 value when VOLKOV_SUPP points at his class_*_aver.txt files.
#
# Settings: JOBS [nproc] generator processes; MAKEJ [1] compile jobs (a split five-loop kernel takes
# ~4 GB and ~4 minutes to compile); THREADS [all] for the run; OUT [runs5] for results and logs;
# KEXT [kernels-ext]; CLASSES ["g2 g5"]; VOLKOV_SUPP [unset].
# Needs python3 with numpy and sympy >= 1.13, and FORM on the PATH or in $FORM (the script stops without it: the
# Python Dirac fallback takes hours per five-loop graph).  Run from the top of the tree.
# Cost on 4 cores (h = 2.4): D = 12 graphs 5-100 min at eps = 1e-6, about 10x more at 1e-8; the D = 13
# graph tried took 3 h at 1e-6.  Memory during the run stays below 100 MB.
set -e
[ $# -ge 1 ] || { sed -n '2,24p' "$0"; exit 1; }
G=$1 H=${2:-2.4} EPS=${3:-1e-7}
KEXT=${KEXT:-kernels-ext} OUT=${OUT:-runs5} CLASSES=${CLASSES:-"g2 g5"} WINDOW=${WINDOW:-48}
JOBS=${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)} MAKEJ=${MAKEJ:-1}
ROOT=$(cd "$(dirname "$0")/.." && pwd); cd "$ROOT"
info=$(python3 scripts/g5info.py "$G" $VOLKOV_SUPP)
NAME=$(echo "$info" | cut -f1) VOLK=$(echo "$info" | cut -f2) D=$(echo "$info" | cut -f3)
REF=$(echo "$info" | cut -f6) SU=$(echo "$info" | cut -f7) SD=$(echo "$info" | cut -f8)
echo "== $NAME  '$VOLK'  D=$D  w=$(echo "$info" | cut -f4)  terms=$(echo "$info" | cut -f5)  Volkov 2019: $REF (+$SU -$SD)"
[ "$D" -le 13 ] || { echo "D=$D exceeds MAXD-1 (driver/sq.h)"; exit 1; }
mkdir -p "$KEXT/g5" "$OUT"
K="$KEXT/g5/$NAME.c"
if [ ! -s "$K" ]; then
    python3 -c "from sympy.polys.matrices import DomainMatrix as M; M.adj_det" 2> /dev/null || { echo "SymPy >= 1.13 needed (DomainMatrix.adj_det): pip install -U 'sympy>=1.13'"; exit 1; }
    command -v "${FORM:-form}" > /dev/null || { echo "FORM executable '${FORM:-form}' not found (set FORM=/path/to/form): without it the five-loop Dirac tables take many hours"; exit 1; }
    T=$(mktemp -d "${TMPDIR:-/tmp}/run5.XXXXXX"); mkdir -p "$T/g5"
    echo "-- generating $NAME (FORM O1)"
    ( cd generator && GEN_FORM_LEVEL=O1 PYTHONHASHSEED=0 python3 gen.py 5 "$T/g5" --jobs "$JOBS" --only "$VOLK" ) | grep -v "table_5.h" || true
    python3 tools/vecg.py -i "$T/g5/$NAME.c"
    sh tools/optimize.sh "$T" > /dev/null
    python3 tools/splitform.py "$T/g5/$NAME.c"
    mv "$T/g5/$NAME.c" "$K"; rm -rf "$T"
fi
python3 tools/splitform.py "$KEXT"/g5/*.c > /dev/null      # idempotent: kernels made by hand earlier
if ! bin/sincqed list 2>/dev/null | grep -q "^ *$NAME "; then
    echo "-- building ($CLASSES)"
    make -j"$MAKEJ" KEXT="$KEXT" CLASSES="$CLASSES" QOPT=O0 > "$OUT/build.log" 2>&1 || { tail -5 "$OUT/build.log"; exit 1; }
fi
N=$(awk -v h="$H" -v w="$WINDOW" 'BEGIN{ n = w/h; if(n > int(n)) n = int(n) + 1; print 2*n + 1 }')
TAG="${NAME}_h${H}_eps${EPS}"
echo "-- running $NAME at h=$H N=$N eps=$EPS (results $OUT/$TAG.txt)"
[ -n "$THREADS" ] && export OMP_NUM_THREADS=$THREADS
( cd "$OUT" && SINCQED_RESULTS="$TAG.txt" SINCQED_LOG="$TAG.log" "$ROOT/bin/sincqed" run "$H" "$N" "$EPS" "$NAME" ) | grep "^$NAME \|resuming" || true
V=$(awk -v g="$NAME" '$1 == g { for(i = 2; i <= NF; i++) if($i ~ /^[-+][0-9]\.[0-9]+e[-+][0-9]+$/ && length($i) > 18){ v = $i; break } } END{ print v }' "$OUT/$TAG.txt")
# pruning: the driver's prune estimate (weight dropped below eps) and the eps-shadow extrapolation.  At five
# loops the shadow assumes a geometric fall-off in eps that the g5_3086 eps ladder does not show (it
# understated the loss 10x at eps = 1e-7), so the prune estimate is the one to compare with.
P=$(sed -n "s/^$NAME .* prune \([-+0-9.e]*\) .*/\1/p" "$OUT/$TAG.log" | tail -1)
SH=$(sed -n "s/.*$NAME eps [a-z]*.*pruning loss left at eps ~ \([-+0-9.e]*\).*/\1/p" "$OUT/$TAG.log" | tail -1)
# the same pass also sums what eps*10 and eps*100 would have kept: the values one and two decades up
D12=$(sed -n "s/.*$NAME eps [a-z]*: eps\*10 would change the value by \([-+0-9.e]*\), eps\*100 by a further \([-+0-9.e]*\).*/\1 \2/p" "$OUT/$TAG.log" | tail -1)
awk -v v="$V" -v p="${P:-0}" -v sh="${SH:--}" -v dd="$D12" -v r="$REF" -v su="$SU" -v sd="$SD" -v g="$NAME" -v h="$H" -v e="$EPS" 'BEGIN{
    printf "%s h=%s eps=%s  value %+.6f  prune estimate %.1e  (eps-shadow estimate %s)", g, h, e, v, p, sh
    if(r != "-") printf "  Volkov %+.6f(+%s/-%s)  diff %+.2e (%.2f%%) = %.2f x prune estimate", r, su, sd, v - r, 100*(v - r)/r, (p > 0 ? (v - r)/p : 0)
    if(split(dd, d, " ") == 2) printf "\n    same run at eps x10: %+.6f   x100: %+.6f", v + d[1], v + d[1] + d[2]
    printf "\n" }' | tee -a "$OUT/summary.txt"
