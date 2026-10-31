#!/bin/sh
# regen_update.sh -- regenerate every kernel with the current generator and install the result in the tree.
#
#   sh scripts/regen_update.sh OUTDIR [JOBS]
#
# For a generator change that alters the kernels (e.g. the GEN_HBASIS=2 default, 2026-10-31).  Steps:
#   1. scripts/regen_all.sh OUTDIR JOBS: g2 g3 v3 v4 w4 pi6 lbl4 g4 into OUTDIR/<class>/, OUTDIR/SHA256SUMS
#      and OUTDIR.tar.xz (resumable for g4; ~30 core-hours of Dirac tables from an empty cache, plus
#      emission, which GEN_HBASIS=2 makes ~6x slower in FORM);
#   2. every five-loop kernel in $KEXT/g5 that scripts/g5info.py recognises (g5_<n>.c) again, in ONE gen.py
#      call (FORM O1, JOBS processes, shared Dirac cache), then vectorised, optimised and split one by one,
#      into OUTDIR/g5/ (skip with G5=0);
#   3. installs: the classes in git go to kernels/ (and kernels/SHA256SUMS is rewritten for them), g4 lbl4 g5
#      to $KEXT/.  The previous kernels are moved to OUTDIR/previous/ first, never deleted;
#   4. make clean; make -j$MAKEJ (MAKEJ [1]: one compile at a time, five-loop kernels need ~4 GB each).
# Afterwards: sh tests/check_regen.sh (must now match), sh scripts/check_g4.sh, then commit kernels/ and
# publish OUTDIR.tar.xz (g4 lbl4) as the new release archives.
# Settings: KEXT [kernels-ext], MAKEJ [1], G5 [1], INSTALL [1] (0: steps 1-2 only), BUILD [1].
set -e
OUT=${1:?usage: sh scripts/regen_update.sh OUTDIR [JOBS]}; JOBS=${2:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}
KEXT=${KEXT:-kernels-ext} MAKEJ=${MAKEJ:-1} G5=${G5:-1} INSTALL=${INSTALL:-1} BUILD=${BUILD:-1}
ROOT=$(cd "$(dirname "$0")/.." && pwd); cd "$ROOT"
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
echo "== GEN_HBASIS=${GEN_HBASIS:-2 (default)}"
t0=$(date +%s)

# 1. two-, three- and four-loop classes
sh scripts/regen_all.sh "$OUT" "$JOBS"

# 2. the five-loop kernels already in use
if [ "$G5" = 1 ] && ls "$KEXT"/g5/g5_*.c > /dev/null 2>&1; then
    T=$OUT/.g5raw; rm -rf "$T"; mkdir -p "$T/g5" "$OUT/g5"; : > "$T/only"
    for K in "$KEXT"/g5/g5_*.c; do
        N=$(basename "$K" .c)
        info=$(python3 scripts/g5info.py "$N" 2> /dev/null) || { echo "-- $N: not a generator graph number, left as it is"; continue; }
        [ "$(echo "$info" | cut -f1)" = "$N" ] || { echo "-- $N: name mismatch, left as it is"; continue; }
        echo "$info" | cut -f2 >> "$T/only"
    done
    n=$(wc -l < "$T/only" | tr -d ' ')
    echo "== regenerating $n five-loop kernels (FORM O1, $JOBS jobs)"
    if [ "$n" -gt 0 ]; then
        ( cd generator && tr '\n' '\0' < "$T/only" | GEN_FORM_LEVEL=O1 PYTHONHASHSEED=0 xargs -0 python3 gen.py 5 "$T/g5" --jobs "$JOBS" --only ) | grep -v "table_5.h" || true
        for K in "$T"/g5/g5_*.c; do
            N=$(basename "$K" .c); echo "-- post-processing $N"
            rm -rf "$T/one"; mkdir -p "$T/one/g5"; mv "$K" "$T/one/g5/"
            python3 tools/vecg.py -i "$T/one/g5/$N.c"
            sh tools/optimize.sh "$T/one" > /dev/null
            python3 tools/splitform.py "$T/one/g5/$N.c"
            mv "$T/one/g5/$N.c" "$OUT/g5/"
        done
    fi
    rm -rf "$T"
fi

# 3. install, keeping the old kernels
if [ "$INSTALL" = 1 ]; then
    mkdir -p "$OUT/previous"
    for d in "$OUT"/*/; do
        c=$(basename "$d"); case $c in previous) continue ;; esac
        ls "$d"*.c > /dev/null 2>&1 || continue
        dst=kernels; case $c in g4|lbl4|g5) dst=$KEXT ;; esac
        mkdir -p "$dst/$c" "$OUT/previous/$c"
        if [ "$c" = g5 ]; then           # only the regenerated g5 files are replaced; hand-made variants stay
            for f in "$d"*.c; do [ -f "$dst/$c/$(basename "$f")" ] && mv "$dst/$c/$(basename "$f")" "$OUT/previous/$c/"; cp "$f" "$dst/$c/"; done
        else
            [ -n "$(ls "$dst/$c" 2>/dev/null)" ] && mv "$dst/$c"/* "$OUT/previous/$c/"
            cp "$d"*.c "$dst/$c/"
        fi
        echo "-- installed $(ls "$d"*.c | wc -l | tr -d ' ') $c kernels in $dst/$c (old ones in $OUT/previous/$c)"
    done
    if command -v sha256sum > /dev/null 2>&1; then SHA=sha256sum; else SHA="shasum -a 256"; fi
    ( cd kernels && $SHA ./*/*.c | sort -k2 > SHA256SUMS )
    echo "-- rewrote kernels/SHA256SUMS ($(wc -l < kernels/SHA256SUMS | tr -d ' ') files)"
fi

# 4. rebuild
if [ "$INSTALL" = 1 ] && [ "$BUILD" = 1 ]; then
    echo "== make clean; make -j$MAKEJ (log $OUT/build.log)"
    make clean > /dev/null 2>&1 || true
    make -j"$MAKEJ" KEXT="$KEXT" > "$OUT/build.log" 2>&1 || { tail -5 "$OUT/build.log"; exit 1; }
    tail -1 "$OUT/build.log"
fi
echo "== done in $(( ($(date +%s) - t0)/60 )) min"
