# Using sincFeyn

This guide covers building and running the `sincqed` program. For how it works see
[IMPLEMENTATION.md](IMPLEMENTATION.md); for the commands behind the paper's numbers see
[REPRODUCE.md](REPRODUCE.md).

## Build

```
make -j8                                  # all classes present in kernels/ (and kernels-ext/)
make -j8 CLASSES="g2 g3 v3"               # everything up to three loops
make -j8 CLASSES="g2 g3 v3 lbl4"          # only some classes: smaller, faster link
make -j8 CC=gcc-14                        # macOS: a specific Homebrew GCC (default: the newest on the PATH)
make test                                 # checks against known values (a few minutes)
```

No Python is needed. Each class becomes a static library `build/libsq_<class>.a`, and
`bin/sincqed` links the driver with the libraries you chose, so a change to one kernel recompiles
only that file.

| class | contents | files |
|---|---|---|
| g2   | two loops, and the Pi2 insertion in the one-loop vertex | 5 |
| g3   | three loops without lepton loops, and three-loop light-by-light (g3_lbl) | 29 |
| v3   | three-loop vacuum-polarisation insertions (w3_p4* are a Pi4 cross-check) | 18 |
| g4   | four loops without lepton loops (release archive, or `generator/regen.sh`) | 269 |
| v4   | four-loop Pi2 insertions | 91 |
| w4   | four-loop Pi4 insertions | 27 |
| pi6  | the proper Pi6 insertion | 27 |
| lbl4 | four-loop light-by-light (release archive) | 19 |

Build options:
- `QOPT` (default `O1`): optimisation of the high-precision half of each kernel, which is rarely
  executed. `QOPT=O0` compiles the four-loop kernels fastest.
- `DD` (default 1): the high-precision fallback is double-double (32 digits from FP64
  arithmetic). `DD=0` uses `__float128` (x86-64).
- The largest lbl4 kernels need about 13 GB of memory per compile job; use a lower `-j`.
- `GPU=1`: the CUDA build, see [GPU.md](GPU.md).

## Run

```
bin/sincqed list
bin/sincqed run h N [eps] [top] GRAPH...     # GRAPH: name, Volkov string or group
bin/sincqed compare RESULTS_h RESULTS_1.2h   # the realistic error: |S(h)-S(1.2h)|/3 + round
bin/sincqed scale GRAPH                      # finiteness (scaling) test of the subtracted integrand
```

- **Groups:** all2, all3, all4, q3, vp3, vp4, lbl, lbl4, pi6.
- **N** = 2K+1 points per axis, so the window is |z| ≤ Kh; K is rounded up to a multiple of
  2^(top+1). Points are stored in 8 bits, so N ≤ 255.
- **top:** the coarsest level is 2^top h. The default (1 for h > 0.75, 2 below) keeps it at most
  3 in z (`SINCQED_TOPSTEP`). A coarser start misses whole patches of the integrand that no eps
  recovers.
- **Output** goes to the screen, `sincqed.log`, and one line per graph in
  `sincqed_results.txt` (columns in [IMPLEMENTATION.md](IMPLEMENTATION.md), Sec. 4).
- **Checkpoints.** Every 10 minutes (`SINCQED_CKPT` seconds, 0 = off) a run writes
  `sincqed_<graph>.ckpt` (in `SINCQED_CKPT_DIR`, default the current directory). If the run is
  killed, the same command resumes from it. The file is used only if h, N, eps and top all match,
  and it is removed when the graph finishes.

## Choosing h, N and eps

- **Window.** Use |z| ≤ 40–48, i.e. N = 2 ceil(48/h) + 1 (h = 1.5: N = 65; 1.0: 97). For two-loop
  values below about 1e-10 widen it to |z| ≤ 72. The sum visits only the support, so a wide window
  costs little; a narrow one (|z| ≤ 24) truncates seven-dimensional graphs at the 1e-4 level.
- **Step.** The error falls like exp(−c/h). Typical accuracy of the sum of all graphs at an
  order: two loops 2e-6 at h = 1.0 and 1e-9 at h = 0.7; three loops 3e-4 at h = 1.2, 4e-6 at
  h = 1.0, 2.5e-7 at h = 0.8, 2e-8 at h = 0.7.
- **Threshold.** eps is relative to the reference scale of each graph. The pruning loss is
  one-signed and adds up across graphs, and it grows as h shrinks, so lower eps with h: two loops
  1e-15 to 1e-20, three loops 1e-11 to 1e-15, four loops 1e-10 to 1e-11 at h = 1.5 with
  `SINCQED_PREDICT=1`.
- **Judging convergence.** Compare runs at two step sizes. The eps shadow line in the log gives
  the eps dependence of each graph from a single run, and an estimate of the pruning loss still
  left.

## Settings (environment variables)

| variable | default | effect |
|---|---|---|
| `SINCQED_PREDICT=alpha` | off | skip a final-level point when the mean \|f\| at the corners of its 2h cube is below alpha × the keep threshold; use with eps ten times smaller. Large gains at four loops. |
| `SINCQED_SYM` | on | fold over axis permutations that leave the integrand unchanged (found automatically; `0` = off, or name the involutions as `p0,p1,...`) |
| `SINCQED_PACK` | on | the vectorised kernels evaluate several lattice points per call (`0` = one at a time) |
| `SINCQED_DDSLACK=s` | 1 | send a point to the double-double fallback only when the bound exceeds s × tolerance; 100 removes most fallbacks at four loops with negligible change |
| `SINCQED_ABREL`, `G2_ABSERR` | 1e-3, 1e-15 | the per-point fallback tolerance is max(G2_ABSERR, ABREL·eps·Ref) |
| `SINCQED_ABSKIP` | 1 | skip the \|.\| sums at points whose parent was far from needing the fallback |
| `SINCQED_DROPUNRES` | 1 | count points unresolved even in double-double as zero (reported in the log) |
| `SINCQED_GROW=0` | | scan the whole coarse window instead of following the support |
| `SINCQED_HIST=prefix` | | per-level histograms of the term sizes |
| `SINCQED_DUMP=prefix` | | write the evaluated points (coarse levels, and a 1-in-`SINCQED_DUMPSAMPLE` sample of the final level) |
| `SINCQED_LOG`, `SINCQED_RESULTS` | sincqed.log, sincqed_results.txt | output files (`-` = none) |
| `SINCQED_PROGRESS` | | seconds between progress lines |

## Run scripts

| script | what it runs |
|---|---|
| `scripts/scan2.sh` | the two-loop graphs over a range of h (wraps `scripts/convergence.sh`; exact values in `scripts/refs.tsv`) |
| `scripts/run3_hscan.sh [OUTDIR]` | every three-loop graph on a ladder of h (`LADDER="h:eps ..."`), resumable; `scripts/sum3.sh` tabulates the totals |
| `scripts/run3_hfine.sh [OUTDIR]` | the fine three-loop ladder, h = 1.50 ... 1.00 in steps of 0.01, then 0.95 ... 0.80 |
| `scripts/check_g4.sh` | checks a g4 kernel tree (labels, dimensions, weights) and, with a binary, four reference values |
| `scripts/run_groupV_h1.3_gpu.sh` | group V at h = 1.3 on a GPU |
| `scripts/regen_all.sh OUTDIR` | regenerate every kernel class |
| `scripts/regen_update.sh OUTDIR` | regenerate every kernel class and the five-loop kernels in use, install them (old ones kept in OUTDIR/previous), rewrite kernels/SHA256SUMS, rebuild |
| `scripts/run5.sh GRAPH [H] [EPS]` | one five-loop graph, from nothing to a value compared with Volkov |
| `scripts/run5batch.sh LIST [H] [EPS]` | a list of five-loop graphs: missing kernels generated together at the start, one serial build, then each run |

The header of each script lists its environment variables.

## Tried and dropped

Grid rotations and per-axis step scales, a tail map z = L sinh(t/L), sparse grids, refinement
cells pointing outward or inward, edge-cell trimming, and log-linear, quadratic, weighted and
inner-first versions of the point predictor were all implemented and removed: each lost to the
plain grid with a smaller eps at equal cost.
