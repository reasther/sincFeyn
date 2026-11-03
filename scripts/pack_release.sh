#!/bin/sh
# pack_release.sh -- the kernel archives attached to a GitHub release, from kernels-ext/.
#
#   sh scripts/pack_release.sh [OUTDIR]          (default release/; run from the top of the tree)
#
# Writes OUTDIR/kernels-<class>.tar.xz for each of g4, lbl4 and g5 present in $KEXT [kernels-ext], with
# paths kernels-ext/<class>/... so that `tar xJf kernels-<class>.tar.xz` at the top of a checkout puts
# them where the Makefile looks, plus OUTDIR/SHA256SUMS over the archives and, inside each archive,
# kernels-ext/<class>/SHA256SUMS over its .c files.  Check g4 first with scripts/check_g4.sh.
# Then:  gh release create vX.Y.Z OUTDIR/kernels-*.tar.xz OUTDIR/SHA256SUMS --title ... --notes ...
set -e
OUT=${1:-release}; KEXT=${KEXT:-kernels-ext}
ROOT=$(cd "$(dirname "$0")/.." && pwd); cd "$ROOT"
mkdir -p "$OUT"
sha(){ if command -v sha256sum >/dev/null; then sha256sum "$@"; else shasum -a 256 "$@"; fi; }
n=0
for c in g4 lbl4 g5; do
  d="$KEXT/$c"; ls "$d"/*.c >/dev/null 2>&1 || { echo "-- no kernels in $d, skipped"; continue; }
  ( cd "$d" && sha *.c > SHA256SUMS )
  # COPYFILE_DISABLE: no macOS ._ files.  The archive path is always kernels-ext/<class>.
  ( cd "$KEXT/.." && COPYFILE_DISABLE=1 tar cf - --transform "s,^$(basename "$KEXT")/,kernels-ext/," \
      "$(basename "$KEXT")/$c"/*.c "$(basename "$KEXT")/$c/SHA256SUMS" 2>/dev/null \
    || COPYFILE_DISABLE=1 tar cf - -s ",^$(basename "$KEXT")/,kernels-ext/," "$(basename "$KEXT")/$c"/*.c "$(basename "$KEXT")/$c/SHA256SUMS" ) \
    | xz -9 -T0 > "$OUT/kernels-$c.tar.xz"
  echo "-- $OUT/kernels-$c.tar.xz  $(ls "$d"/*.c | wc -l | tr -d ' ') files, $(du -h "$OUT/kernels-$c.tar.xz" | cut -f1)"
  n=$((n + 1))
done
[ $n -gt 0 ] || { echo "nothing packed"; exit 1; }
( cd "$OUT" && sha kernels-*.tar.xz > SHA256SUMS ); cat "$OUT/SHA256SUMS"
