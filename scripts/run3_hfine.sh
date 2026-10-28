#!/bin/sh
# run3_hfine.sh -- fine h scan of all three-loop graphs (all3, 41 entries = 72 diagrams) on the A100.
# Run from the top of the tree:
#   BIN=bin/sincqed-gpu PAR=4 sh scripts/run3_hfine.sh [OUTDIR]          (default OUTDIR = runs3_fine)
# Ladder: h = 1.50, 1.49, ..., 1.00 (51 values), then 0.95, 0.90, 0.85, 0.80.  Coarsest first.
# eps per band is the 2026-10-09 A100 scan's ladder divided by EPSDIV (default 10, one decade tighter):
#   h >= 1.35: 1e-10   1.15-1.34: 1e-11   0.95-1.14: 1e-12   0.85-0.94: 1e-13   0.80: 1e-14   (before /EPSDIV)
# At EPSDIV=1 the 2026-10-09 scan had pruning loss <= 6e-7, far below the h error at every h; /10 makes
# that ~5x smaller still.  Estimated A100 time (PAR=4, from the 2026-10-09 timings): ~10 h at EPSDIV=1,
# ~13 h at EPSDIV=10.  Window |z| <= 48 (ZMAX); N is rounded so the window is 48..48+4h, a truncation
# change of ~1e-9 between neighbouring h, below everything else.
# Everything else (resume, PAR, BUDGET, GRAPHS, output layout) is as in scripts/run3_hscan.sh.
set -u
EPSDIV=${EPSDIV:-10}
LADDER=$(awk -v d="$EPSDIV" 'BEGIN{
  for (i = 150; i >= 100; i--) h[++n] = i/100; h[++n] = 0.95; h[++n] = 0.90; h[++n] = 0.85; h[++n] = 0.80;
  for (j = 1; j <= n; j++) { x = h[j];
    e = (x >= 1.345) ? 1e-10 : (x >= 1.145) ? 1e-11 : (x >= 0.945) ? 1e-12 : (x >= 0.845) ? 1e-13 : 1e-14;
    printf "%.2f:%.0e ", x, e/d } }')
export LADDER
exec sh "$(dirname "$0")/run3_hscan.sh" "${1:-runs3_fine}"
