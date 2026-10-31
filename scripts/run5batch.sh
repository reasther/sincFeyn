#!/bin/sh
# run5batch.sh -- a batch of five-loop graphs: make every missing kernel first, build once, then run each.
#
#   sh scripts/run5batch.sh LIST [H] [EPS]
#
# LIST is a file with one graph per line (a Volkov string or a kernel name, as for run5.sh; blank lines
# and lines starting with # are skipped).  Steps:
#   1. look up every graph and list those with no kernel in $KEXT/g5/;
#   2. generate all of the missing kernels in ONE gen.py call (FORM O1), so that the Dirac tables and the
#      per-graph emission run in parallel over JOBS processes and share one table cache (separate gen.py
#      processes would each rewrite generator/cache_dirac.pkl and lose the others' tables); then
#      vectorise, optimise and split each kernel, one at a time;
#   3. build bin/sincqed once, serially (MAKEJ [1]: a split five-loop kernel needs ~4 GB to compile);
#   4. run the graphs one after another with scripts/run5.sh, which finds the kernels and the build in
#      place and skips straight to the run.
# Settings as for run5.sh: JOBS [nproc], MAKEJ [1], THREADS, OUT [runs5], KEXT [kernels-ext],
# CLASSES ["g2 g5"], WINDOW [48], VOLKOV_SUPP.  Run from the top of the tree.
set -e
[ $# -ge 1 ] || { sed -n '2,19p' "$0"; exit 1; }
LIST=$1 H=${2:-2.4} EPS=${3:-1e-7}
KEXT=${KEXT:-kernels-ext} OUT=${OUT:-runs5} CLASSES=${CLASSES:-"g2 g5"}
JOBS=${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)} MAKEJ=${MAKEJ:-1}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
case $LIST in /*) ;; *) LIST=$(pwd)/$LIST ;; esac
cd "$ROOT"; [ -f "$LIST" ] || { echo "no list file $LIST"; exit 1; }
mkdir -p "$KEXT/g5" "$OUT"

# 1. which kernels are missing
names= missing= nmiss=0
while IFS= read -r G || [ -n "$G" ]; do
    case $G in ''|'#'*) continue ;; esac
    info=$(python3 scripts/g5info.py "$G")
    NAME=$(echo "$info" | cut -f1) VOLK=$(echo "$info" | cut -f2) D=$(echo "$info" | cut -f3)
    [ "$D" -le 13 ] || { echo "$NAME: D=$D exceeds MAXD-1 (driver/sq.h), skipped"; continue; }
    case " $names " in *" $NAME "*) continue ;; esac
    names="$names $NAME"
    if [ ! -s "$KEXT/g5/$NAME.c" ]; then missing="$missing
$VOLK"; nmiss=$((nmiss + 1)); fi
done < "$LIST"
echo "== $(echo $names | wc -w) graphs:$names"
echo "== $nmiss kernels to generate"

# 2. generate the missing kernels together, then post-process them one by one
if [ "$nmiss" -gt 0 ]; then
    python3 -c "from sympy.polys.matrices import DomainMatrix as M; M.adj_det" 2> /dev/null || { echo "SymPy >= 1.13 needed (DomainMatrix.adj_det): pip install -U 'sympy>=1.13'"; exit 1; }
    command -v "${FORM:-form}" > /dev/null || { echo "FORM executable '${FORM:-form}' not found (set FORM=/path/to/form)"; exit 1; }
    T=$(mktemp -d "${TMPDIR:-/tmp}/run5batch.XXXXXX"); mkdir -p "$T/g5"
    echo "-- generating $nmiss kernels (FORM O1, $JOBS jobs)"
    echo "$missing" | sed '/^$/d' > "$T/only"
    ( cd generator && tr '\n' '\0' < "$T/only" | GEN_FORM_LEVEL=O1 PYTHONHASHSEED=0 xargs -0 python3 gen.py 5 "$T/g5" --jobs "$JOBS" --only ) | grep -v "table_5.h" || true
    for K in "$T"/g5/g5_*.c; do
        N=$(basename "$K" .c); echo "-- post-processing $N"
        mkdir -p "$T/one/g5"; mv "$K" "$T/one/g5/"
        python3 tools/vecg.py -i "$T/one/g5/$N.c"
        sh tools/optimize.sh "$T/one" > /dev/null
        python3 tools/splitform.py "$T/one/g5/$N.c"
        mv "$T/one/g5/$N.c" "$KEXT/g5/$N.c"
    done
    rm -rf "$T"
fi
python3 tools/splitform.py "$KEXT"/g5/*.c > /dev/null      # idempotent: kernels made by hand earlier

# 3. one serial build if any graph is not in the driver yet
build=0
for N in $names; do bin/sincqed list 2>/dev/null | grep -q "^ *$N " || build=1; done
if [ $build = 1 ]; then
    echo "-- building ($CLASSES, make -j$MAKEJ; log $OUT/build.log)"
    make -j"$MAKEJ" KEXT="$KEXT" CLASSES="$CLASSES" QOPT=O0 > "$OUT/build.log" 2>&1 || { tail -5 "$OUT/build.log"; exit 1; }
fi

# 4. run them in turn
for N in $names; do sh scripts/run5.sh "$N" "$H" "$EPS"; done
