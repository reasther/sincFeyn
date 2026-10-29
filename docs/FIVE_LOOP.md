# Five loops: first tests

Notes from a first attempt at tenth-order (five-loop) graphs without lepton loops, run on a 4-core,
16 GB machine. Nothing here is production quality yet.

## Quick start

```
VOLKOV_SUPP=/path/to/volkov_supplement sh scripts/run5.sh "6; 1-8, 2-9, 3-10, 4-11, 5-7" 2.4 1e-7
```

`scripts/run5.sh` takes one graph (a Volkov string in either orientation, his supplement's
`6;1*8;...` form, or a kernel name such as `g5_3086`) and does everything below: it generates the
kernel, vectorises, optimises and splits it, builds `bin/sincqed` with every kernel in
`kernels-ext/g5`, runs at the given h and eps, and prints the value, the driver's prune estimate (with the eps-shadow
estimate alongside), and the difference from Volkov's 2019 value. Steps whose output exists are skipped,
so a second h or eps for the same graph goes straight to the run. `python3 scripts/g5info.py GRAPH
[SUPPDIR]` prints the kernel name, D, mirror weight, number of subtraction terms and Volkov's value
without running anything. VOLKOV_SUPP is the directory of Volkov's `class_*_aver.txt` files (the
supplement of PRD 100, 096004); the tables are not in this repository. Settings are listed at the
top of the script (JOBS, MAKEJ, THREADS, OUT, WINDOW).

## Generator changes

- `generator/symb.py`: the determinant and adjugate of the loop matrix are now computed over the
  polynomial ring ZZ[s] (`DomainMatrix.adj_det`). SymPy's `Matrix.det()`/`adjugate()` (Bareiss
  with `dotprodsimp`) did not finish a five-loop node in 30 minutes; the ring version takes seconds.
- `generator/gen.py`, `generator/build.py`: the graded arrays carry Wick counts 0..NW-1 with NW
  fixed at 4, enough up to four loops. Five-loop nodes reach a Wick count of 4, which overran the
  arrays (stack smashing in the top F2 node). gen.py now sets `NW = max(4, L)`.
- `tests/check_regen.sh` still regenerates every committed kernel byte for byte.

`enumerate_graphs(5)` gives 3213 integrals (6354 diagrams), matching Aoyama et al.,
arXiv:2412.06473, and Volkov, arXiv:1909.08015.

## Building five-loop kernels

```
cd generator
GEN_FORM_LEVEL=O1 PYTHONHASHSEED=0 python3 gen.py 5 OUTDIR --jobs 2 --only "6; 1-8, 2-9, 3-10, 4-11, 5-7"
cd ..; mkdir -p kernels-ext/g5; cp OUTDIR/g5_*.c kernels-ext/g5/
python3 tools/vecg.py -i kernels-ext/g5/*.c; sh tools/optimize.sh kernels-ext
make -j2 CLASSES="g2 g5" QOPT=O0
```

- FORM `O3` (the default) did not finish the 12-13 MB node expressions in 30+ minutes; `O1`
  takes about 3 minutes per graph. The kernels are 12-22 MB of C.
- New Dirac tables come first and are cached (`cache_dirac.pkl`). With FORM the large ones (~8e5
  structures) took 4-6 minutes each here; without FORM the Python fallback is several hundred times
  slower, so run5.sh refuses to generate unless FORM is found.
- Compiling two kernels takes about 23 minutes and several GB per compile job. Running the
  generator at the same time got the compiler OOM-killed on 16 GB.
- Do not run two gen.py processes that both compute new Dirac tables: they race on
  `cache_dirac.pkl` (use `--cache-dir`, or one process).
- D at five loops (from build.classes): 66 integrals at D = 9, 257 at 10, 638 at 11, 1125 at 12 and
  1127 at 13. MAXD (driver/sq.h) stays 12 unless g5 is in CLASSES, when the Makefile sets 14. The packed 64-bit lattice key is used only for the one
  dilation at the coarsest level (step 2h). At D = 13 it fits for h >= 2; below that the driver uses
  the old sort-based dilate, which is cheap at the ~1e5 cells these runs keep.
- A D = 13 node (g5_1221) is one 8.8 MB function, and gcc needed more than 14 GB to compile it even at
  -O1. `python3 tools/splitform.py FILE` cuts the FORM blocks into chunk functions of 4000
  statements. Two kernels then build in 8 minutes with 4.3 GB peak. The values are bit-identical:
  make test on split g2/g3 kernels, and g5_3086 again gives -0.360975824517. g5 is compiled with
  OPT_K4 (-O2).
- Building recipe, after optimize.sh: `python3 tools/splitform.py kernels-ext/g5/*.c`, then
  `make -j2 CLASSES="g2 g5" QOPT=O0`.

## The test class

The self-energy topology "abcdeeabcd" (one of the 98 classes in arXiv:2412.06473 with a single
second-order self-energy) is mirror symmetric: its 9 vertex graphs are 5 integrals. Reference values
are Volkov's 2019 per-graph values (supplement of PRD 100, 096004, `class_*_aver.txt`). They use
the same subtraction scheme (on-shell U and L, M^2 = m^2) and include the mirror weight.

| kernel | graph | w | Volkov 2019 |
|---|---|---|---|
| g5_3086 | 6; 1-8, 2-9, 3-10, 4-11, 5-7 | 1 | -0.387249(566) |
| g5_1078 | 3; 1-8, 2-9, 4-10, 5-11, 6-7 | 2 | -0.921393(634) |
| g5_372  | 2; 1-8, 3-9, 4-10, 5-11, 6-7 | 2 | +0.148167(611) |
| g5_1784 | 4; 1-8, 2-9, 3-10, 5-11, 6-7 | 2 | +0.148707(628) |
| g5_2490 | 5; 1-8, 2-9, 3-10, 4-11, 6-7 | 2 | +0.771280(641) |
| sum | | | -0.240488 |

## Results (4 threads, double-double fallback, N chosen so that |z| <= 48)

| graph | h | eps | value | evals | time | peak RSS |
|---|---|---|---|---|---|---|
| g5_3086 | 2.4 | 1e-6 | -0.360976 | 9.4e7 | 318 s | 73 MB |
| g5_3086 | 2.2 | 1e-6 | -0.351277 | 1.6e8 | 575 s | 58 MB |
| g5_3086 | 2.0 | 1e-6 | -0.349492 | 2.7e8 | ~1650 s | 52 MB |
| g5_3086 | 2.4 | 1e-7 | -0.370940 | 3.1e8 | 1087 s | 75 MB |
| g5_3086 | 2.4 | 1e-8 | -0.381621 | 8.6e8 | 2959 s | 73 MB |
| g5_3086 | 2.4 | 1e-9 | -0.383597 | | (Mac, 24 thr) | |
| g5_372 | 2.4 | 1e-7 | +0.100746 | | (Mac, 24 thr) | |
| g5_1078 | 2.4 | 1e-6 | -0.948473 | 2.2e8 | 5797 s | 87 MB |
| g5_1221 (D = 13) | 2.4 | 1e-6 | -0.610670 | 4.2e8 | ~3 h | 62 MB |

g5_1221 = 3; 1-10, 2-5, 4-7, 6-9, 8-11 is a primitive D = 13 graph (one term, no subdivergences),
Volkov -0.730246(533). Its eps shadow estimates the pruning loss left at eps = 1e-6 at -0.128,
which gives -0.739, 1.2% from Volkov. |S_2h - S_h| is only 1.2e-2, so at D = 13 the pruning loss is
also much larger than the discretisation error.

The g5_3086 eps ladder at h = 2.4 is -0.370940, -0.381621, -0.383597 for eps = 1e-7, 1e-8, 1e-9:
per-decade changes -0.0107 and -0.0020, so the pruning loss is now falling (ratio 5.4) and extrapolates
to about -0.3840. The remaining 0.8% to Volkov's -0.3872(6) is presumably the h = 2.4 discretisation;
h = 2.2 and 2.0 at eps = 1e-9 are the next runs. g5_372 at eps = 1e-7 (+0.1007 against +0.1482) is a
cancelling 4-term graph and is far from converged in eps.

An independent build on a 24-thread Mac (Homebrew gcc, FORM built from source) gave g5_3086 at h = 2.4,
eps = 1e-7 as -0.370940243117 with the same 306728281 evaluations, identical to the Linux run, in 321 s.

The eps shadow does not extrapolate the pruning loss at five loops: at eps = 1e-7 it predicts -9.5e-4
left, but the next decade alone moves the value by -1.07e-2 (the per-decade moves are -0.0099 and
-0.0107, not falling geometrically). The prune column is the better guide (2.6e-2 at 1e-7, against
1.6e-2 actually left), and run5.sh reports the difference from Volkov in units of it.

(The h = 2.0 run was interrupted and resumed from its checkpoint, so its time is approximate.)

- **Pruning dominates at h ≈ 2.** Each tenfold cut in eps moves g5_3086 by about +0.01 towards
  Volkov's value (-0.3610, -0.3709, -0.3816). The prune column (8.1e-2, 2.6e-2, 6.9e-3) is the
  right size: value + prune at eps = 1e-8 is -0.3885, against -0.3872(6). With eps = 1e-6, the
  h-ladder only follows the pruning loss, so h-dependence should be studied at eps <= 1e-8.
- **Memory is not a limit.** The hierarchical sum keeps 5e4-4e5 cells and the peak RSS stays
  under 90 MB. What sets the cost is time: evaluations grow about 3x per tenfold cut in eps and
  about 1.7x per 10% cut in h.
- **g5_1078 costs about 8x more per evaluation than g5_3086** (26 vs 3.4 µs per evaluation per
  thread). It has twice as many Dirac structures and no symmetry fold; the rest is not yet
  understood.
- The `disc` column is far too large at these h (as at two loops) and should be ignored.
- `sincqed scale` on a five-loop graph means 3^14 directions evaluated in double-double. That did
  not finish in 10 minutes; a sampled version of the test is needed.

## Next steps

1. g5_3086 at eps = 1e-9 (h = 2.4), and at h = 2.0 with eps = 1e-8, to separate the
   discretisation and pruning errors (a few hours on 4 cores; quicker on 24+ threads).
2. Compile g5_372, g5_1784 and g5_2490 and run the whole class against -0.240488.
3. Profile g5_1078, and try FORM O2 on the largest nodes for faster kernels.
4. A sampled scale test for D = 12.
