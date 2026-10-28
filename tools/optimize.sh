#!/bin/sh
# Value-preserving speed-ups of the kernel sources, applied in place.  generator/regen.sh runs this
# by default (REGEN_PLAIN=1 skips it); run it by hand on a tree generated with REGEN_PLAIN=1 or by
# an older version of the generator:
#   sh tools/optimize.sh [KERNELDIR]        (default: kernels)
#  1. compact.py on every kernel: each sub-sum of the factored Dirac polynomials is computed once,
#     as straight-line code in functions of about 1500 terms, and the |.| companion sums go in a
#     block the driver can skip (1.3-4x faster, 2-15x faster to compile);
#  2. glpts.py on every kernel: the rule over the external-vertex split u cut to one point at the
#     weighted mean of u, which is exact because the point function is linear in u (2x);
#  3. packpts.py on the vectorised kernels (g3, g4): a wrapper fv_ that puts different lattice
#     points in the SIMD lanes, which the one-point rule left idle.
# Values agree with the plain kernels to rounding.  The scripts are idempotent, so running this
# twice changes nothing.  Rebuild afterwards (make clean; make ...).
# Files are independent: OPT_JOBS=N (default 1) runs N at a time; the output does not depend on it.
set -e
K=${1:-kernels}; T=$(cd "$(dirname "$0")" && pwd); J=${OPT_JOBS:-1}
for d in "$K"/*/; do
    ls "$d"*.c >/dev/null 2>&1 || continue
    echo "compact $d"; ls "$d"*.c | xargs -n 1 -P "$J" python3 "$T/compact.py" -i > /dev/null
done
echo "glpts $K/*"; ls "$K"/*/*.c | xargs -n 8 -P "$J" python3 "$T/glpts.py" > /dev/null
echo "packpts $K/*"; ls "$K"/*/*.c | xargs -n 8 -P "$J" python3 "$T/packpts.py" > /dev/null
