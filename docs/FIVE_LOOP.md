# Five loops

The paper evaluates a sample of eleven tenth-order (five-loop) vertex graphs without lepton loops and compares
them with Volkov's 2019 Monte Carlo values (Fig. 5). This page describes those runs, the h extrapolation and
pruning correction applied to them, and how the five-loop kernels are made.

## The graphs

`data/g5/graphs.txt` lists them. Reference values, errors and GPU times are Volkov's per-graph values
(supplement of PRD 100, 096004, `class_*_aver.txt`), copied for these eleven graphs into
`data/g5/volkov2019.tsv`. They use the same subtraction scheme as sincFeyn (on-shell U and L, M^2 = m^2) and
include the mirror weight.

| kernel | Volkov graph | D | w | Volkov 2019 | set |
|---|---|---|---|---|---|
| g5_3086 | 6; 1-8, 2-9, 3-10, 4-11, 5-7 | 12 | 1 | -0.387249(566) | class abcdeeabcd |
| g5_1784 | 4; 1-8, 2-9, 3-10, 5-11, 6-7 | 12 | 2 | +0.148707(628) | class abcdeeabcd |
| g5_1078 | 3; 1-8, 2-9, 4-10, 5-11, 6-7 | 12 | 2 | -0.921393(634) | class abcdeeabcd |
| g5_372  | 2; 1-8, 3-9, 4-10, 5-11, 6-7 | 12 | 2 | +0.148167(611) | class abcdeeabcd |
| g5_2490 | 5; 1-8, 2-9, 3-10, 4-11, 6-7 | 12 | 2 | +0.771280(641) | class abcdeeabcd |
| g5_1663 | 4; 1-7, 2-5, 3-11, 6-10, 8-9 | 12 | 2 | -4.441411(1321) | AHKN-Volkov tension classes |
| g5_3024 | 6; 1-7, 2-4, 3-11, 5-10, 8-9 | 12 | 2 | -6.061594(1369) | AHKN-Volkov tension classes |
| g5_717  | 3; 1-4, 2-7, 5-8, 6-11, 9-10 | 12 | 2 | -5.416099(1108) | AHKN-Volkov tension classes |
| g5_2567 | 5; 1-9, 2-7, 3-10, 4-8, 6-11 | 13 | 2 | -0.625916(201) | tension classes, no subtractions |
| g5_1795 | 4; 1-8, 2-10, 3-7, 5-9, 6-11 | 13 | 2 | +0.521242(242) | tension classes, no subtractions |
| g5_1221 | 3; 1-10, 2-5, 4-7, 6-9, 8-11 | 13 | 2 | -0.730246(533) | no subtractions |

The class abcdeeabcd (one of the 98 self-energy classes of arXiv:2412.06473 with a single second-order
self-energy) is all five integrals of one gauge-invariant self-energy class; its Volkov sum is -0.240488.
The tension classes are those where Aoyama et al. (2024) and Volkov (2019) differ most (X382, X227, X352, X350).

## The recipe

Every run is top = 0 (the top-down sum, `sincqed run H N EPS 0 GRAPH`) with N = 2 ceil(48/h) + 1.

1. **h scan at matched eps.** h = 2.4, 2.3, 2.2, 2.1 (and 2.0 for the class graphs) with
   eps = 1e-6 (h/2.4)^D. Scaling eps with the number of points keeps the pruned fraction about the same at
   every h, so the scan sees only the h dependence. Mixing eps across h spoils the fit.
2. **h extrapolation.** Least squares S(h) = S_inf + A exp(-c/h) over all the matched-eps points (c on a grid,
   S_inf and A linear). Equation 7 of the paper motivates the form; the fitted c is 13-23. An exact
   three-point fit is ill-conditioned (a 3e-4 wobble in one point moves S_inf by 0.02), so four or five h
   are used.
3. **Pruning correction.** One more run at h = 2.4, eps = 1e-5. The loss from pruning falls geometrically
   with eps, by a ratio r = 0.24 per decade at D = 12 (Aitken on eps 1e-5, 1e-6, 1e-7 at h = 2.4 agrees
   to ~5e-4 for g5_717, 1663, 3024), so delta = r/(1-r) (S(1e-6) - S(1e-5)) is added to S_inf.
   At D = 13 the measured ratio is larger (0.28 for g5_2567 from 1e-6 to 1e-7, still rising), so the
   rule under-corrects there; `K13=...` in fit5.py changes the coefficient for D = 13.

Every run that enters a result is charged in the cost of Fig. 5.

## Reproducing Fig. 5 and the table below

The results files of the runs used in the paper are in `data/g5/` (one `bin/sincqed` results line each).
From them, with no runs:

```
python3 scripts/fit5.py                 # per-graph extrapolated values against Volkov (the table below)
python3 scripts/plot5_cost.py           # Fig. 5 (g5_cost_mean.png) and the per-graph version (g5_cost.png)
```

These runs are heterogeneous: they were made between 31 October and 3 November 2026 on several machines
(a 24-thread Mac Studio M2 Ultra, a 12-thread laptop and a desktop) with successive versions of the driver,
all of them top = 0 at the settings above. The figure in the preprint is drawn from them, and its costs
are scaled from their evaluation counts (see `scripts/plot5_cost.py`), not from wall time. The final
version of the paper will use runs made with the tagged release, which the commands below reproduce.
`data/g5/README.md` lists what is in the directory.

To run everything again (FORM and python3 with numpy and sympy needed for the kernels):

```
sh scripts/run5batch.sh data/g5/graphs.txt 2.4 1e-6    # generate, split and build the 11 kernels, run h = 2.4
sh scripts/run5_hscan.sh                                # h = 2.3, 2.2, 2.1 at matched eps, and h = 2.4 at 1e-5
python3 scripts/fit5.py runs5
python3 scripts/plot5_cost.py runs5
```

The kernels take a few hours to generate and build; the runs take about 20 hours on a 24-thread machine.

## Results

`python3 scripts/fit5.py` on `data/g5/` (eps rule r = 0.24 for every D; sigma is Volkov's, on the side facing
our value):

| graph | D | h used | c | S_inf | eps corr | result | Volkov | diff | sigma |
|---|---|---|---|---|---|---|---|---|---|
| g5_3086 | 12 | 2.4-2.0 | 12.9 | -0.386595 | -1.0e-3 | -0.387594 | -0.387249 | -0.09% | -0.9 |
| g5_1784 | 12 | 2.4-2.0 | 22.7 | +0.148866 | +6.6e-4 | +0.149522 | +0.148707 | +0.55% | +1.3 |
| g5_1078 | 12 | 2.4-2.0 | 18.1 | -0.918099 | -3.1e-3 | -0.921171 | -0.921393 | +0.02% | +0.3 |
| g5_372  | 12 | 2.4-2.0 | 13.8 | +0.149514 | +2.4e-4 | +0.149758 | +0.148167 | +1.07% | +2.6 |
| g5_2490 | 12 | 2.4-2.0 | 17.9 | +0.771342 | +9.3e-4 | +0.772267 | +0.771280 | +0.13% | +1.5 |
| g5_1663 | 12 | 2.4-2.1 | 19.6 | -4.429827 | -1.0e-2 | -4.440009 | -4.441411 | +0.03% | +1.1 |
| g5_3024 | 12 | 2.4-2.1 | 15.8 | -6.042537 | -1.5e-2 | -6.057191 | -6.061594 | +0.07% | +3.2 |
| g5_717  | 12 | 2.4-2.1 | 17.2 | -5.401530 | -1.1e-2 | -5.412350 | -5.416099 | +0.07% | +3.4 |
| g5_2567 | 13 | 2.4-2.1 | 22.9 | -0.621272 | -3.1e-3 | -0.624332 | -0.625916 | +0.25% | +7.9 |
| g5_1795 | 13 | 2.4-2.1 | 20.3 | +0.516881 | +2.8e-3 | +0.519654 | +0.521242 | -0.31% | -10.4 |
| g5_1221 | 13 | 2.4-2.1 | 16.5 | -0.722909 | -4.6e-3 | -0.727514 | -0.730246 | +0.37% | +5.1 |

The mean relative difference is 2.7e-3 (Volkov's mean relative sigma is 1.0e-3). The three D = 13
residuals are each 0.52-0.59 of their eps correction and on the side of under-correction, which is what
the larger D = 13 eps ratio predicts. Per-diagram mean cost of the projected values is 1.04 TFLOP-hours at
nominal peak, against 32 for Volkov's Monte Carlo (Fig. 5).

## Running one graph

```
VOLKOV_SUPP=/path/to/volkov_supplement sh scripts/run5.sh "6; 1-8, 2-9, 3-10, 4-11, 5-7" 2.4 1e-6
```

`scripts/run5.sh` takes one graph (a Volkov string in either orientation, his supplement's `6;1*8;...` form,
or a kernel name such as `g5_3086`). It generates the kernel, vectorises, optimises and splits it, builds
`bin/sincqed` with every kernel in `kernels-ext/g5`, runs at the given h and eps, and prints the value and,
with VOLKOV_SUPP set, the difference from Volkov's 2019 value. Steps whose output exists are skipped.
`python3 scripts/g5info.py GRAPH [SUPPDIR]` prints the kernel name, D, mirror weight, number of subtraction
terms and Volkov's value without running anything. For several graphs, `sh scripts/run5batch.sh LIST H EPS`
generates all the missing kernels in one `gen.py` call (one shared Dirac-table cache), builds once and runs
them one after another. Settings are listed at the top of each script.

## Generator changes for five loops

- `generator/symb.py`: the determinant and adjugate of the loop matrix are now computed over the
  polynomial ring ZZ[s] (`DomainMatrix.adj_det`). SymPy's `Matrix.det()`/`adjugate()` (Bareiss
  with `dotprodsimp`) did not finish a five-loop node in 30 minutes; the ring version takes seconds.
- `generator/gen.py`, `generator/build.py`: the graded arrays carry Wick counts 0..NW-1 with NW
  fixed at 4, enough up to four loops. Five-loop nodes reach a Wick count of 4, which overran the
  arrays (stack smashing in the top F2 node). gen.py now sets `NW = max(4, L)`.
- `tests/check_regen.sh` still regenerates every committed kernel byte for byte.

`enumerate_graphs(5)` gives 3213 integrals (6354 diagrams), matching Aoyama et al.,
arXiv:2412.06473, and Volkov, arXiv:1909.08015.

## Building five-loop kernels by hand

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
  1127 at 13. MAXD (driver/sq.h) stays 12 unless g5 is in CLASSES, when the Makefile sets 14. The packed lattice key (mixed radix, W = N at the
  key's step) is used by the coarsest-level grow and dilation. It is 64 bits while W^D <= 2^63 and
  128 bits beyond (driver/sumkey.h, included twice by sum.c); the 64-bit path gives the same points and
  sums as before, and SINCQED_KEY128=1 forces the 128-bit one for testing. Before the 128-bit keys a
  top = 0 run (no coarse level, so the grow runs at step h) needed N <= 37 for D = 12 and N <= 28 for
  D = 13, and silently fell back to a full box scan otherwise.
- top = 0 (no coarse level; `sincqed run H N EPS 0 GRAPH`): the grow then runs at step h and expands
  down to eps*1e-3, evaluating only the kept points and a one-point boundary shell, with no refinement.
  At five loops this removes most of the pruning loss: g5_3086 at h = 2.4, top = 0, eps 1e-6 gave
  -0.383953 (173M evaluations, 504 s on 24 threads) where top = 1 needed eps ~1e-9 (top = 1 at 1e-7:
  -0.370940 from 306M evaluations; Volkov -0.387249).  The refinement's 2^D children per kept cell,
  mostly below threshold, are what top = 1 spends its evaluations on.  The grow evaluates each frontier
  through eval_batch (packed SIMD kernel calls, same values; SINCQED_GROWBATCH=0 for one point at a
  time) and, at top = 0 only, skips the |.| sums for points reached from a checked, non-near point
  within SINCQED_GROWSKIP generations (default 2; 0 = check all).  g5_3086 h = 2.4 N = 37 eps 1e-4,
  4 cores: 825 s one at a time, 190 s batched, 137 s batched + skip (57% skipped; 114 s at GROWSKIP=4),
  identical values; SINCQED_ABCHECK=1 found 0 skipped points that needed the fallback (2 at GROWSKIP=4).
- A D = 13 node (g5_1221) is one 8.8 MB function, and gcc needed more than 14 GB to compile it even at
  -O1. `python3 tools/splitform.py FILE` cuts the FORM blocks into chunk functions of 4000
  statements. Two kernels then build in 8 minutes with 4.3 GB peak. The values are bit-identical:
  make test on split g2/g3 kernels, and g5_3086 again gives -0.360975824517. g5 is compiled with
  OPT_K4 (-O2).
- Building recipe, after optimize.sh: `python3 tools/splitform.py kernels-ext/g5/*.c`, then
  `make -j2 CLASSES="g2 g5" QOPT=O0`.
