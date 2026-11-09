# sincFeyn

[![DOI](https://zenodo.org/badge/1436884324.svg)](https://doi.org/10.5281/zenodo.23975727)

The code implements Sinc-function based evaluations of Feynman diagrams. It currently focuses on 
the QED vertex graphs that contribute to the electron anomalous magnetic moment, at two, three and
four loops, and a sample of five-loop graphs but the overall summation code can be extended to 
other scenarios. 

This code should be regarded as "beta" software and, as currently constituted, is primarily intended
to be a "proof of concept" implementation of Sinc function integration for nontrivial diagrams. 

The author welcomes feedback and suggestions! 

Each graph is written as an integral over Schwinger parameters. Its ultraviolet and infrared
divergences are removed point by point with Volkov's forest formula. The integral is then
evaluated as a pruned, hierarchical trapezoid ("Sinc") sum in z = log s, whose error falls
exponentially in 1/h. 

The method and the results are described in the accompanying paper. This repository holds
everything needed to reproduce its numbers:
- the integration driver (C, OpenMP; optional CUDA back end);
- the generated graph kernels (C);
- the Python generator that writes those kernels, including the Dirac algebra (FORM when it is
  installed);
- the source transformations that make the kernels fast;
- the run scripts.


## What it computes

| order | diagrams | integrals | check |
|---|---|---|---|
| two loops | 7 | 5 | every graph against its exact value (Volkov's subtraction scheme) |
| three loops | 72 | 41 | the sum against A1^(6) (Laporta–Remiddi): -1.7e-8 at h = 0.7 |
| four loops | 891 | 423 | group V (518 diagrams without lepton loops) and every lepton-loop set against Laporta ([LEPTON_LOOPS.md](docs/LEPTON_LOOPS.md)) |
| five loops | 12672 (6354 without lepton loops) | 3213 without lepton loops | generator covers all; 11 graphs run against Volkov's values ([FIVE_LOOP.md](docs/FIVE_LOOP.md)) |

Every diagram is computed from scratch; the literature is used only for comparison.

## Quick start

You need a C/C++ compiler with OpenMP (GCC on Linux; Homebrew GCC on macOS). No Python is
needed to build or run.

```
make -j8 CLASSES="g2 g3 v3"        # everything up to three loops (a few minutes)
make test                          # four checks against analytic values (about 2 minutes)
bin/sincqed list                   # the registered graphs
bin/sincqed run 1.0 97 1e-12 all2  # the seven two-loop diagrams at h = 1.0
```

The program is called `sincqed`. On macOS use `make CC=gcc-14` (or whichever Homebrew GCC you
have). The usage guide is [docs/USAGE.md](docs/USAGE.md). With an NVIDIA GPU, `make GPU=1`
builds `bin/sincqed-gpu`, which evaluates the integrand on the GPU; see [docs/GPU.md](docs/GPU.md).

## Four- and five-loop kernels

The four-loop group-V kernels (g4, 269 files), the four-loop light-by-light kernels (lbl4) and the
five-loop kernels of the paper (g5, 11 files) are too large for the repository. They come as
archives attached to the [v0.9.0 release](https://github.com/reasther/sincFeyn/releases/tag/v0.9.0)
(with `SHA256SUMS`), which unpack into `kernels-ext/` (not tracked by git; the Makefile reads them
there, `KEXT=DIR` to keep them elsewhere).

They can also be regenerated from scratch (FORM needed). The four-loop sets take a few hours on 16
cores:

```
sh generator/regen.sh kernels-ext 16 g4 lbl4
```

A five-loop kernel, with full optimisation, takes around 20 minutes of single-threaded FORM time on a
desktop; [docs/FIVE_LOOP.md](docs/FIVE_LOOP.md) gives the commands.


## Documentation

| file | contents |
|---|---|
| [docs/IMPLEMENTATION.md](docs/IMPLEMENTATION.md) | how a graph becomes a Sinc sum: parametrisation, subtractions, the hierarchical sum, precision |
| [docs/USAGE.md](docs/USAGE.md) | building, running, settings, choosing h, N and eps |
| [docs/REPRODUCE.md](docs/REPRODUCE.md) | the commands behind the numbers in the paper |
| [docs/FIVE_LOOP.md](docs/FIVE_LOOP.md) | the five-loop graphs, their h extrapolation, and how five-loop kernels are made |
| [docs/GENERATOR.md](docs/GENERATOR.md) | the kernel generator and its checks |
| [docs/LEPTON_LOOPS.md](docs/LEPTON_LOOPS.md) | the four-loop lepton-loop classes and their exact values |
| [docs/GPU.md](docs/GPU.md) | the CUDA back end |
| [mathematica/](mathematica/README.md) | the two-loop integrands and Sinc sums as a Mathematica notebook, for checking the C by hand |

## Layout

```
driver/      the integration driver: hierarchical Sinc sum, precision fallback, checkpoints, GPU
kernels/     generated C, one directory per class (g2 g3 v3 v4 w4 pi6; g4 lbl4 from the release)
generator/   Python: graphs and forests, Dirac tables, kernel emission, lepton-loop generators
tools/       kernel transformations (shared sub-sums, one-point rule, SIMD), GPU/SIMD code generators, checks
scripts/     run scripts for the paper's scans, Table I and Figs. 2, 4, 5; release packing
data/        the runs behind Figs. 4 and 5 (g3, g4, g5) and Volkov's five-loop reference values
tests/       regression runs (make test) and the regeneration checks
docs/        documentation
mathematica/ the two-loop Sinc sums in Mathematica, translated from the C kernels
pslq/        the high-precision two-loop runs (quad-double) behind the PSLQ closed forms of Table I
```

## How this code was written

The code in this repository was written almost entirely by an AI coding assistant (Claude, by
Anthropic), working under the direction of Richard Easther, who specified the method, chose the
tests and decided what to keep. The evidence that it is right is the set of independent checks
described in the paper: agreement with the analytic two- and three-loop coefficients and with
exactly known gauge-invariant subsets, pointwise finiteness (scaling) tests of every subtracted
integrand, and external checks of the summations for low dimension graphs. The information needed
to reproduce the results can be found in this repository.

## Method

The overall method builds on work by Easther, Guralnik, Hahn, Petrov and Wang, from 1999-2001
([hep-ph/9903255](https://arxiv.org/abs/hep-ph/9903255),
[hep-ph/9912255](https://arxiv.org/abs/hep-ph/9912255),
[hep-ph/0010143](https://arxiv.org/abs/hep-ph/0010143)).


## Citing

If you use this code, please cite the paper and the code. The code's Zenodo DOI,
[10.5281/zenodo.23975727](https://doi.org/10.5281/zenodo.23975727), always resolves to the latest
release; `CITATION.cff` has the details (GitHub's "Cite this repository" button reads it).

```
@misc{sincFeyn,
  author       = {Easther, Richard},
  title        = {{sincFeyn}},
  year         = {2026},
  howpublished = {\url{https://github.com/reasther/sincFeyn}},
  note         = {Zenodo, doi:10.5281/zenodo.23975727}
}
```

## Licence

MIT, see [LICENSE](LICENSE).
