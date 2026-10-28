#!/bin/sh
# Regenerate the committed kernels from the Python generator and check that they are byte-identical
# to kernels/ (via kernels/SHA256SUMS).  Needs python3 with numpy and sympy, and FORM 4.3 (g2, g3,
# v3, v4 are FORM-optimised).  About 6 minutes on
# 4 cores from an empty cache.
#   sh tests/check_regen.sh [OUTDIR] [JOBS]
set -e
OUT=${1:-${TMPDIR:-/tmp}/sincfeyn_regen.$$}; JOBS=${2:-4}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
sh "$ROOT/generator/regen.sh" "$OUT" "$JOBS"
cd "$OUT"
if command -v sha256sum >/dev/null 2>&1; then SUM="sha256sum -c --quiet"; else SUM="shasum -a 256 -c --quiet"; fi
if $SUM "$ROOT/kernels/SHA256SUMS"; then echo "regenerated kernels match kernels/SHA256SUMS ($OUT)"; else echo "MISMATCH (see above)"; exit 1; fi
