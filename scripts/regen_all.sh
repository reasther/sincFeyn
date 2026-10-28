#!/bin/sh
# Regenerate every kernel class (g2 g3 v3 v4 w4 pi6 lbl4 g4) with the current generator, on a
# workstation:
#   sh scripts/regen_all.sh OUTDIR [JOBS]
# Needs FORM (brew install form / apt install form), python3 with numpy and sympy, and
# generator/cache_lbl.pkl (the light-by-light Dirac tables: without it lbl4 recomputes them in
# Python, tens of core-hours).  Writes OUTDIR/<class>/*.c, OUTDIR/SHA256SUMS (sorted, paths
# <class>/<file>) and OUTDIR.tar.xz.  Nothing in kernels/ is touched; to use the result (g4 and lbl4,
# which are not in git, go to kernels-ext/):
#   for c in OUTDIR/*/; do c=$(basename $c); d=kernels; case $c in g4|lbl4) d=kernels-ext;; esac
#     mkdir -p $d; rm -rf $d/$c; cp -r OUTDIR/$c $d/; done; make clean; make
# Resumable for g4 (finished tables and kernels are kept in OUTDIR/.raw); the other classes take
# minutes and are always redone.
set -e
OUT=${1:?usage: sh scripts/regen_all.sh OUTDIR [JOBS]}; JOBS=${2:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
command -v "${FORM:-form}" >/dev/null || { echo "FORM not found (brew install form / apt install form, or set FORM=...)"; exit 1; }
python3 -c 'import numpy, sympy' 2>/dev/null || { echo "python3 needs numpy and sympy (pip3 install numpy sympy)"; exit 1; }
[ -f "$ROOT/generator/cache_lbl.pkl" ] || echo "WARNING: generator/cache_lbl.pkl missing: lbl4 Dirac tables will be recomputed in Python (slow)"
echo "FORM: $(${FORM:-form} -v 2>&1 | head -1)"; echo "python: $(python3 --version), sympy $(python3 -c 'import sympy; print(sympy.__version__)')"
t0=$(date +%s)
sh "$ROOT/generator/regen.sh" "$OUT" "$JOBS" g2 g3 v3 v4 w4 pi6 lbl4 g4
OUT=$(cd "$OUT" && pwd)
( cd "$OUT" && if command -v sha256sum >/dev/null; then sha256sum */*.c; else shasum -a 256 */*.c; fi | sort -k2 > SHA256SUMS )
( cd "$(dirname "$OUT")" && tar -cJf "$(basename "$OUT").tar.xz" --exclude .raw "$(basename "$OUT")" )
echo "== $(ls "$OUT"/*/*.c | wc -l) kernels, $(du -sh "$OUT" | cut -f1) (with .raw), $(( ($(date +%s) - t0)/60 )) min"
echo "== archive $OUT.tar.xz"
