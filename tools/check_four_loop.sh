#!/bin/sh
# Regenerate the four-loop group-V (g4) and light-by-light (lbl4) kernels from scratch (with the
# default speed-ups, as in the production tree), and compare them file by file with an
# existing kernel tree (for example the one a production run was built from).
#
#   sh tools/check_four_loop.sh EXISTING_KERNELS_DIR [JOBS] [OUTDIR]
#
# EXISTING_KERNELS_DIR is the directory that holds g4/ and lbl4/ (e.g. ~/sinc/sincqed/kernels).
# Cost: the Dirac tables are about 30 core-hours for g4 and of the same order for lbl4, then
# C emission; run it from a fresh checkout so that no cached tables are reused.
# Writes OUTDIR/SHA256SUMS.four (checksums of the regenerated g4 and lbl4 files) and a summary
# OUTDIR/check_four_loop.txt.  Needs python3 with numpy and sympy.
set -e
OLD=${1:?usage: sh tools/check_four_loop.sh EXISTING_KERNELS_DIR [JOBS] [OUTDIR]}
[ -d "$OLD" ] || { echo "no directory $OLD"; exit 1; }
OLD=$(cd "$OLD" && pwd)    # absolute: the comparison below runs from inside OUTDIR
JOBS=${2:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
OUT=${3:-./regen4}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
sh "$ROOT/generator/regen.sh" "$OUT" "$JOBS" g4 lbl4
cd "$OUT"
if command -v sha256sum >/dev/null 2>&1; then SH=sha256sum; else SH="shasum -a 256"; fi
find g4 lbl4 -name '*.c' | sort | xargs $SH > SHA256SUMS.four
{
  echo "# check_four_loop $(date)  existing: $OLD"
  for c in g4 lbl4; do
    same=0; diff=0; miss=0
    for f in $c/*.c; do
      b=$(basename "$f")
      if [ ! -e "$OLD/$c/$b" ]; then miss=$((miss+1)); echo "missing in existing tree: $c/$b"
      elif cmp -s "$f" "$OLD/$c/$b"; then same=$((same+1))
      else diff=$((diff+1)); echo "differs: $c/$b"; fi
    done
    extra=$(ls "$OLD/$c" 2>/dev/null | grep -c '\.c$' || true)
    echo "$c: $same identical, $diff differ, $miss missing  (existing tree has $extra .c files)"
  done
} | tee check_four_loop.txt
