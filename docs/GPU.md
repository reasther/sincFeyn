# CUDA back end (`make GPU=1`)

`make GPU=1` builds `bin/sincqed-gpu` next to the CPU build `bin/sincqed`.

The algorithm is unchanged. The same hierarchical Sinc sum, pruning, predictor and checkpoints evaluate the
same points, and sum them in the same order. What changes is where the work runs:
- the integrand is evaluated on the GPU;
- the host work around it overlaps with the GPU;
- the double-double fallback runs on the CPU, vectorised (8 points at a time with AVX-512).

GPU kernels are built for the classes `g2 g3 v3 g4 g5`. The other classes (`v4 w4 pi6 lbl4`) run on
the CPU in `bin/sincqed-gpu`.

Developed and tested on an NVIDIA A100 80GB (sm_80) with a 28-core Xeon (Cascade Lake, AVX-512) and 93 GB
of RAM, using CUDA 13.3 and GCC 13.3.

## Build and run

```
make GPU=1 -j16 CLASSES="g2 g3 v3"       # bin/sincqed-gpu (the CPU kernel libraries are built and linked too)
sh tests/run_tests.sh bin/sincqed-gpu    # the regression checks with the GPU build

BIN=bin/sincqed-gpu PAR=4 MODE=matched sh scripts/run3_hscan.sh    # the three-loop h-scan, four graphs at a time
```

Four loops, group V (all 269 g4 graphs): unpack the g4 release archive into
`kernels-ext/` first, then

```
sh scripts/check_g4.sh                         # kernels-ext/g4 is the 269 FORM kernels (labels, D, weights), nothing else
make GPU=1 -j16 CLASSES="g2 g3 g4"             # g2 and g3 are for the regression tests
sh tests/run_tests.sh bin/sincqed-gpu
sh scripts/check_g4.sh bin/sincqed-gpu         # g4_73, g4_72, g4_66, g4_184 against the Mac Studio values to 1e-7
BIN=bin/sincqed-gpu HS=1.3 sh scripts/run_groupV.sh   # rerun the same command to resume
```

- **Requirements.** The build needs the CUDA toolkit and Python 3 (for `tools/gpugen.py` and
  `tools/ddvgen.py`) in addition to GCC. The CPU build is unchanged and needs neither.
- **Build time.** About a minute for `g2 g3 v3` here.

`make GPU=1` options:
- `GPUARCH`: the target GPU. The default `native` is the GPU in the machine; `sm_80` is an A100, `sm_90` an H100.
- `CUDA`: install prefix (default `/usr/local/cuda`).
- `NVOPT`: nvcc flags (default `-O3 -arch=$(GPUARCH) -fmad=false -diag-suppress 177`, the last silencing the unused-temporary warnings of the FORM kernels). Keep `-fmad=false`: with nvcc's default
  fused multiply-adds, 11 g4 kernels lose exact cancellations between terms (g4_191, 149, 90, 258, 54, 122,
  03, 39, 241, 174, 01; g4_191 gave -1.46e7 instead of -0.0516). Without fusing, all 269 g4 graphs agree with the
  CPU to 7.9e-9 at h=3 N=33 eps=1e-8. It costs 5-10%.
- `GPUCLASSES`: classes that get GPU kernels (default: those of `g2 g3 v3 g4 g5` in `CLASSES`).
- `DDVGRAPHS`, `DDVOPT`: the graphs that get the SIMD double-double fallback (default: all of them), and its
  compiler flags (default `-O1 -march=native`).

Run-time environment variables:
- `SINCQED_GPU=0` runs `bin/sincqed-gpu` on the CPU only (identical to `bin/sincqed`). `SINCQED_GPU=d` picks
  device d-1.
- `SINCQED_DDV=0` uses the scalar double-double fallback.
- `SINCQED_GPUCELLS=0` makes the last level's children on the host (the node path), as a check.
- `SINCQED_TIMING=1` prints a per-graph breakdown of host time.
- `SINCQED_DDSLACK` and the other `SINCQED_*` settings work as in the CPU build.

**Several graphs at once.** With the GPU, `PAR=4` (four graphs at a time) is fastest for the three-loop
scan, because the processes keep the GPU busy while the others are in host-only phases. The processes
must share the cores: `run3_hscan.sh` gives each `OMP_NUM_THREADS` =
cores/`PAR` unless it is already set. With every process on all 28 cores, the h = 1.5 scan took 160 s at
`PAR=4`. Split 4 × 7 it takes 20 s, and `PAR=1` takes 61 s (2026-11-05, same totals). Mind the memory,
though. At h=0.7 each of the 17 large D=7 graphs of all3 needs about 35 GB, with about 55 GB peaks in
`dilate()`. On a 93 GB machine those must run one at a time (`GRAPHS=... PAR=1`).

## Validation

At h = 1.5, eps = 1e-10, all 41 entries of all3, comparing `bin/sincqed` with `bin/sincqed-gpu`:
- **Same work.** All 41 graphs have identical evaluation counts and identical fallback counts.
- **Values.** They agree to within 3.8e-12 relative.
- **Time.** The CPU build (28 threads) took 239 s; the GPU build with `PAR=4` took 55 s.

`tests/run_tests.sh` gives identical values on both builds, and the SIMD double-double fallback
is bit-identical to the scalar `dd.h` (checked with the fallback forced on up to 47% of the
points). With `-fmad=false` all 269 g4 graphs agree with the CPU to 7.9e-9 at h = 3, N = 33,
eps = 1e-8.

**Small differences at very small eps.** At about eps ≤ 1e-13 the keep threshold reaches the
per-point round-off tolerance. The points near it are the heavily cancelling ones, where the CPU
and GPU doubles disagree in the last bits, so the two builds can keep slightly different point
sets (for g3_00 at h = 1.5, eps = 1e-14: evaluation counts differ by 2 in 10^6, values agree to
12 digits). Two CPUs or compilers would differ the same way.

## How it works

### The integrand on the GPU (`tools/gpugen.py`, `driver/gpu_kernel.h`, `driver/gpu.cu`)

The high-precision (QREAL) half of each kernel is a complete scalar copy of the integrand, written against
`REAL`. `gpugen.py` compiles it with `REAL = double` into `build/gpu/<graph>.cu`:
- **Chunk functions** become `__device__ __noinline__`. Inlining them gave 92 KB stack frames.
- **Every other function** takes the global `SQ_AB` as an extra parameter and passes it on at every call.
- **The Gauss-Legendre tables** go into `__constant__` memory.

A batch kernel evaluates one lattice node per thread. It computes `w = prod(s h)` in the order of
`lattice_node()` and returns `f*w` and `sum|terms|*w`. On the four-loop kernels the A100 sustains about
2.4e7 evaluations/s, a quarter of its FP64 peak. The limit is local-memory traffic: the per-thread `T_[]`
sub-sum array, and register spills.

### Host and GPU in step (`driver/eval.c`, `driver/sum.c`)

`eval_gpu_submit` / `eval_gpu_finish` replace `eval_batch` when a graph has a GPU kernel:
- **submit** packs the nodes as int16, copies them, launches the kernel and returns.
- **finish** waits, then applies exactly the tests of `leafx()`: the near flag, and the fallback test with
  `DDSLACK`.

The refinement runs in blocks of cells, two in flight on two CUDA streams. Block j+1 is prepared and
submitted before block j is finished and summed:

```
prep(0); submit(0)
loop:  prep(j+1); submit(j+1);  wait(j); check + fallback + sum(j)
```

All host work stays in one thread in cell order. The sums, the kept set and the checkpoints are therefore
those of a serial pass, whatever the batch size. The frontiers of the coarse grow are GPU batches too.

### The grow on the GPU (`grow()` in `driver/sumkey.h`)

The grow (the coarse level, and the whole sum at top = 0) sends each frontier to the GPU in chunks of
`SQ_GROWCHUNK` (4M) points, two in flight: chunk c+1 is unpacked and submitted before chunk c is
finished, so the host's tests and double-double fallback for one chunk overlap the GPU's evaluation of
the next. At top = 0 the |.|-sum skip (`SINCQED_GROWSKIP`) is on, as in the CPU build: the skipped points
are evaluated without the |.| sums on the GPU too. Values, evaluation counts and skip counts are the
same as the CPU build's (checked with a CPU stand-in for `gpu.cu`, also with 64k-point chunks so that
every layer takes several). The hash bookkeeping between frontiers stays on the CPU and does not overlap
the GPU.

## Five loops (`g5`)

The g5 kernels go through `gpugen.py` like the others. `tools/splitform.py` cuts their large FORM
blocks into pieces (`<node>_cB_K`); `gpugen.py` makes these real calls (`__noinline__`) rather than
inlining them, which for g5_3086 (sm_80, CUDA 13.4) takes the kernel from a 51 KB stack frame and
119 KB of spill stores to 31 KB and 18 KB (plus 1-3 KB per piece). One g5 kernel compiles in about
2.5 minutes with nvcc using 2.8 GB, so `make -j` is limited by the CPU kernels (about 4 GB each), not by
nvcc.

```
make GPU=1 -j2 CLASSES="g2 g5" QOPT=O0          # g5 kernels from kernels-ext/g5 (MAXD = 14)
SINCQED_GPU=1 bin/sincqed-gpu run 2.4 41 1e-5 0 g5_3086
```


### Children made on the GPU at the last level (`prepcells()` in `driver/sum.c`)

On the plain grid (no predictor, dump or histogram), a cell's children at the last level are its 2^D corners
`k - st b`. Only the cells go to the GPU (`gpu_submit_cells`), and each thread computes its child's
coordinates. The host checks the results and sums them in prep()'s order: cell, then b ascending, skipping
children below the box. The cell weight is 1 at the last level, so the sums are exact.

This path is bit-identical to the node path in the same binary. It replaces about 1e10 host-made 48-byte
nodes per three-loop graph at h=0.9. On g3_00 at h=1.1 it cut the wall time from 80 s to 54 s, and the CPU
time by half.

### The double-double fallback, vectorised (`tools/ddvgen.py`, `driver/ddv.h`)

`ddvgen.py` makes a SIMD version of each kernel's double-double half: `REAL = ddv`, a pair of 8-wide
AVX-512 vectors (4 with AVX2). One call evaluates 8 fallback points. `ddv.h` performs the operations of
`dd.h` lane by lane: the same two-sum and fma two-prod, the sloppy add unless `DD_IEEE_ADD`, and dd op
double without promoting the double. So each lane gives what `dd.h` gives.

This works because the half has no branches on values. The only branches are on `SQ_AB`, a constant, and
`RABS`, which becomes a lane-wise select. The exception is the Pi4 parts (`v3_p4*`, `w3_*`), which branch
on computed coefficients. They get an empty stub, and their weak registry entry sends them to the scalar
fallback. On the four-loop g4_181 (1.2M fallback points at h=2) the SIMD fallback cut the graph from 117 s
to 34 s.

Making the double-double fallback a GPU kernel was tried and dropped. nvcc's ptxas runs out of memory on
the double-double four-loop kernels: 94 GB, also at `-O0`.

### Identically zero integrands

`box_zero()` probes z = 0 and 4096 lattice nodes before the grow and reports a graph such as
`v3_p4Vc` as zero without scanning the box. It works unchanged in both builds; the probes run on the CPU.
