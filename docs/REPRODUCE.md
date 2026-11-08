# Reproducing the results

The commands below produce the numbers in the paper. Times are wall-clock. A run with the same
settings and thread count gives bitwise-identical sums; a different thread count or machine
changes only the last digits.

**Preprint figures.** The figures in the preprint are drawn from heterogeneous runs: they were made over
several weeks on different machines (Apple M2 Max and M2 Ultra, x86 desktops and laptops, an NVIDIA A100)
with successive versions of the driver, as the code was developed. Rerunning the commands here with this
release gives the same values to well within the accuracies plotted, but not bit for bit, and different
times. The five-loop runs behind Fig. 5 are included in `data/g5/` (see [FIVE_LOOP.md](FIVE_LOOP.md)).
The figures in the final version of the paper will be redrawn from runs of the tagged release, which
these commands then reproduce exactly.

| paper | section here |
|---|---|
| Table I, Fig. 2 | [2. Two loops](#2-two-loops-table-i-fig-2) |
| Fig. 4 top | [3. Three loops](#3-three-loops-fig-4-top) |
| Fig. 4 bottom | [4. Four loops, group V](#4-four-loops-group-v-fig-4-bottom) |
| Fig. 5 | [7. Five loops](#7-five-loops-fig-5) |

Build first (Linux; on macOS add `CC=gcc-14` or your Homebrew GCC):

```
make -j8 CLASSES="g2 g3 v3"
make test
```

`make test` runs four checks against analytic values: the two-loop sum, the Pi2 insertion,
three-loop graph 16 and three-loop light-by-light. All four should PASS in about two minutes on
4 cores.

## 1. Kernels from the generator

```
python3 -m pip install numpy sympy
sh tests/check_regen.sh
```

This regenerates every committed kernel (two loops, three loops, and the four-loop
vacuum-polarisation classes v4, w4, pi6) from an empty Dirac-table cache and compares them with
`kernels/SHA256SUMS`. It needs FORM 4.3 on the PATH for the classes made with FORM. Expected:
"regenerated kernels match".

## 2. Two loops (Table I, Fig. 2)

Table I (h = 0.55, N = 249), with the residuals already divided out (about 20 seconds):

```
sh scripts/table1.sh
```

Every two-loop graph is compared with its exact value in Volkov's subtraction scheme (the
closed forms are in Table I and in `driver/graphs.c`). Fig. 2 is the same comparison over a range of h;
`sh scripts/fig2.sh` runs h = 1.60 ... 0.55 in steps of 0.01 (about 6 minutes on 4 cores) and writes the
signed residuals of every term and of A1^(4) to `fig2/fig2.tsv` for plotting. The self-energy residual
changes sign between h = 1.21 and 1.20, and the A1^(4) residual between 1.10 and 1.09. By hand:

```
for h in $(seq 1.60 -0.02 0.60); do
  N=$(awk -v h=$h 'BEGIN{n = 72/h; print 2*(n == int(n) ? n : int(n)+1) + 1}')
  SINCQED_RESULTS=two_loop_h$h.txt bin/sincqed run $h $N 1e-20 1 all2
done
```

The window is |z| ≤ 72: with |z| ≤ 48 the crossed and ladder graphs stop near 1e-10 below
h = 0.7. The last column of each results file is the difference from the exact value; the sum
of the seven diagrams differs from A1^(4) by +9e-7 at h = 1.0, +2.3e-9 at h = 0.7 and +4.3e-11
at h = 0.6. The whole ladder takes a few minutes on 4 cores.

## 3. Three loops (Fig. 4)

All 41 three-loop entries (72 diagrams) at h = 2.10, 2.05, ..., 1.00, once at fixed eps = 1e-15 and once
with eps matched to the h error, eps = 1e-11 (h/1.5)^10:

```
MODE=fixed sh scripts/run3_hscan.sh          # runs3_fixed/:   about 3.7 hours on a 24-thread M2 Ultra
MODE=matched sh scripts/run3_hscan.sh        # runs3_matched/: about 1.2 hours
sh scripts/sum3.sh runs3_matched             # totals per h against A1^(6)
```

The runs in the figure are in `data/g3/` and `python3 scripts/fig4.py` tabulates them (with group V, below).
At h = 1.00 the matched total is 1.1812457980 (eps extrapolation included), +4.3e-6 from A1^(6).
Smaller h works the same way with an explicit ladder, e.g. `LADDER="0.70:1e-15" sh scripts/run3_hscan.sh runs3`;
the D = 7 graphs need up to about 55 GB of memory each at h = 0.7, so run them one at a time on a machine
with less than about 128 GB (`GRAPHS=... PAR=1`). Each level writes
`OUTDIR/h<h>/results.txt`; the gauge-invariant sets are sums over its entries:

| set | entries |
|---|---|
| no closed lepton loop (50 diagrams) | g3_00 ... g3_27 |
| Pi2 in the fourth-order vertex (12) | v3_01 ... v3_07 |
| light-by-light (6) | g3_lbl |
| vacuum polarisation in the second-order vertex (4) | v3_00, v3_p4Sc, v3_p4Sm, v3_p4Vc, v3_p4Vm |

Exact values: A1^(6) = 1.181241456587 (Laporta and Remiddi), light-by-light 0.3710052921, and
the vacuum-polarisation set 0.0554291774123 (the Pi4 and Pi2Pi2 insertions). The other two sets
are known exactly only as their sum, 0.7548069870749.

On a GPU: `BIN=bin/sincqed-gpu PAR=4 MODE=matched sh scripts/run3_hscan.sh` ([GPU.md](GPU.md)).

## 4. Four loops, group V (Fig. 4)

### Kernels

Fetch the g4 archive from the release and unpack it into `kernels-ext/` (not tracked by git; the
Makefile reads g4 there). Alternatively, regenerate it (FORM needed; a few hours on 16 cores):

```
sh generator/regen.sh kernels-ext "$(nproc)" g4
```

Then check the tree and build:

```
sh scripts/check_g4.sh
make clean; make -j16 CLASSES="g4" QOPT=O0
```

A full g4 build takes a few hours. `tools/check_four_loop.sh EXISTING_KERNELS` regenerates g4
and lbl4 from scratch and compares them file by file with an existing tree.

### Runs

```
sh scripts/run_groupV.sh                                  # h = 1.5, N = 57: about 14 hours on a 24-thread M2 Ultra
HS="2.2 2.1 1.9 1.7 1.3" sh scripts/run_groupV.sh         # the other points of Fig. 4 (about 7 + 50 hours)
python3 scripts/sum_groupV.py sincqed_results_groupV_h1.5.txt sincqed_groupV_h1.5.log
```

eps is 1e-10 for 199 graphs and 1e-11 for the other 70, both scaled by (h/1.5)^10 (`P=0` keeps them fixed).
The total is the sum over the 269 graphs of the value plus the pruning loss left at eps, which each graph's
`eps extrapolation` log line reports. At h = 1.5 it is -2.175054, +0.0018 from Laporta's -2.176866; at
h = 1.3, -2.174897. The runs in the figure are in `data/g4/`; `python3 scripts/fig4.py` writes the
Fig. 4 table (`fig4.tsv`: total, relative error and TFLOP-h per set and h) from `data/g3/` and `data/g4/`.

## 5. Four-loop lepton-loop sets

These classes are coded and checked against Laporta's values ([LEPTON_LOOPS.md](LEPTON_LOOPS.md)),
but the paper does not report production values for them:

```
bin/sincqed run 1.0 97 1e-10 v4_00                       # set 17, Pi2^3
bin/sincqed run 1.0 97 1e-10 w4_00Vm w4_00Sm w4_00Sc     # set 18
bin/sincqed run 1.0 97 1e-10 w4_0iVm w4_0iSm w4_0iSc     # set 19
bin/sincqed run 1.0 97 1e-10 pi6                         # set 20, proper Pi6
bin/sincqed run 1.2 81 1e-10 lbl4                        # sets 21-25, light-by-light
```

## 6. Finiteness (scaling) tests

```
bin/sincqed scale g3_15              # every direction in {-1,0,1}^(n-1)
SCALE_E=2 bin/sincqed scale g3_27    # {-2..2}
```

A subtracted integrand passes when there is no direction of slow decay. These tests were run for
all two- and three-loop graphs and for the four-loop light-by-light classes.

## 7. Five loops (Fig. 5)

Eleven graphs without lepton loops, each an h scan at matched eps (h = 2.4 to 2.1, top = 0) plus one
eps step, extrapolated to h -> 0 and eps -> 0 and compared with Volkov's 2019 values. From the
runs included in `data/g5/` (the preprint's figure):

```
python3 scripts/fit5.py                 # per-graph values against Volkov; mean relative difference 2.7e-3
python3 scripts/plot5_cost.py           # g5_cost_mean.png (Fig. 5) and g5_cost.png (per graph)
```

From scratch (FORM needed for the kernels; about a day on a 24-thread machine):

```
sh scripts/run5batch.sh data/g5/graphs.txt 2.4 1e-6
sh scripts/run5_hscan.sh
python3 scripts/fit5.py runs5
python3 scripts/plot5_cost.py runs5
```

[FIVE_LOOP.md](FIVE_LOOP.md) gives the graphs, the extrapolation recipe and the results table.

