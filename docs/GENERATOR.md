# The kernel generator

`generator/` holds the Python that writes every kernel in `kernels/`. You only need it to
regenerate or extend the kernels; building and running need no Python. It requires Python 3
(3.11 to 3.13 tested) with NumPy and SymPy, and FORM 4.3 for the classes whose Dirac algebra
FORM computes.

## Pipeline

`generator/regen.sh OUTDIR [JOBS] [CLASS ...]` runs the whole pipeline and leaves
`OUTDIR/<class>/*.c`. With no CLASS it builds the committed classes g2 g3 v3 v4 w4 pi6.

1. **Graphs, forests, terms** (`graphs.py`).
   - Enumerates the vertex graphs without lepton loops as pairings of the 2L non-external
     positions on the electron line. It keeps the 1PI graphs and one representative per mirror
     pair (weight 2) or self-mirror graph (weight 1). The order, and so the file numbering
     g<L>_NN, is deterministic.
   - Finds the UV-divergent subgraphs (contiguous segments closed under the photon arcs), the
     forests, and the subtraction terms of Volkov's formula (operators A, U, L).
2. **Dirac tables.** For each node (subgraph, children, operator), the integer F2 coefficient,
   or the rational vertex and self-energy coefficients, of every numerator structure.
   - `dirac.py` obtains them numerically: explicit 4×4 gamma matrices at sample kinematic
     points, a polynomial fit in t, and rounding to exact integers or rationals with assertions
     ([IMPLEMENTATION.md](IMPLEMENTATION.md), Sec. 2).
   - `formdirac.py` computes the same tables algebraically with FORM (symbolic traces) and writes
     each node polynomial as FORM `Format O3` code (Horner scheme plus common subexpressions).
     Every expression is checked against the exact polynomial at a random rational point before
     it is emitted. The two routes agree exactly on every two- and three-loop table and on the
     largest four-loop table checked.
   - `regen.sh` uses FORM for the classes made by `gen.py` and `vp.py` (g2, g3, g4, v3, v4);
     their kernels carry `Dirac algebra: FORM O3` in the header. `GEN_FORM=0` selects the
     numerical route (gen.py warns at four loops, where it takes hours and gives slower kernels).
     The lepton-loop generators take their tables from `dirac.py`.
   - Tables are cached in `generator/cache_dirac.pkl` and `generator/cache_lbl.pkl`; delete the
     caches to recompute from scratch.
3. **Symbolic polynomials and emission** (`symb.py`, `build.py`).
   - Computes the Symanzik polynomial U, the exponent V, the shifted-momentum coefficients and
     the contraction coefficients of each node with SymPy.
   - **Line identifications.** A node's coefficients depend only on its momentum routing, so
     lines with equal loop and external rows share one symbol, and equal monomials are summed
     before emission.
   - **Same-momentum classes.** Lines with the same momentum at q = 0 share one Sinc variable σ:
     each of the K lines gets σ/K, with weight σ^(K−1)/(K−1)!. This is exact because U, V and the
     contraction coefficients see only the sum of their parameters. It includes lines of one row
     on opposite sides of the external vertex.
   - **Cheng–Wu photon.** The photon of arc 0 has its parameter set to 1, except in a self-mirror
     graph whose arc 0 is not its own mirror; there the first arc that is its own mirror is used,
     so the mirror permutes the remaining Sinc axes and the driver's symmetry fold applies.
   - Writes the term algebra and the wrapper over the Sinc variables.
4. **Lepton loops.**

   | generator | kernels |
   |---|---|
   | `vp.py` | Pi2 insertions: v2_00, v3_00..07, and v4_* with argument 4 |
   | `vp4.py` | Pi4 in the one-loop vertex: v3_p4* |
   | `vp4x.py` | Pi4 (with Pi2) in any vertex graph: w4_* |
   | `pi6.py` | the proper Pi6: w4_p6*; with argument `pi4`, the Pi4 cross-check w3_p4* |
   | `lbl.py` | three-loop light-by-light: g3_lbl |
   | `lbl4b.py`, `lbl4c.py`, `lbl4d.py`, `lblvp.py` | four-loop light-by-light IV(b), IV(c), IV(d), IV(a) |

   In the light-by-light kernels the Cheng–Wu parameter is one that the open-line mirror keeps in
   place, so the driver folds over the mirror as well as the loop reversal.
5. **Vectorisation** (`tools/vec4.py`, `tools/vecg.py`): a SIMD copy of the double half.
6. **Speed-ups** (`tools/optimize.sh`), below.

Steps 1–6 are deterministic. With `PYTHONHASHSEED=0` and node code inline (`GEN_SPLIT=0`,
except pi6.py), the output is byte-identical to the committed kernels; `tests/check_regen.sh`
verifies this (about 6 minutes on 4 cores from an empty cache). For g4 and lbl4, `regen.sh` also
sets `GEN_CHUNK=0`, the layout used for the four-loop runs. These two classes take tens of
core-hours, mostly Dirac tables (819 node tables for g4: about 30 core-hours numerically, under
an hour with FORM). `tools/check_four_loop.sh` regenerates both and compares them with an
existing tree.

## Speed-ups (on by default)

The last step of `regen.sh` runs `tools/optimize.sh`, which applies value-preserving rewrites to
every kernel. They are idempotent, and the committed kernels have them.
- `compact.py`: every sub-sum of the factored trees is computed once, as straight-line code in
  functions of about 1500 terms. The absolute-value sums go in a block the driver can skip.
- `glpts.py`: the rule in u (the split at the external vertex) becomes one point at the
  weighted mean of u, which is exact because the point function is linear in u.
  `tests/check_onept.sh` compares it with the two-point Gauss–Legendre rule on every kernel.
- `packpts.py` (the vectorised g3 and g4 kernels): adds the wrapper `fv_`, which puts a different
  lattice point in each SIMD lane.

`REGEN_PLAIN=1 sh generator/regen.sh OUTDIR` stops before this step and leaves the plain
factored form, which is easier to read and to compare with the generator. `sh tools/optimize.sh
OUTDIR` converts it later.

## Generator options

| variable or option | effect |
|---|---|
| `GEN_FORM=0` | numerical Dirac tables and factored trees instead of FORM |
| `GEN_SPLIT=1` (the default in build.py) | one C function per node: lower compile memory, about 20% slower kernels |
| `GEN_SPLIT=0` | node code inline, as in the committed kernels |
| `GEN_CHUNK=n` | split factored expressions longer than n characters into functions (default 20000) |
| `GEN_HBASIS=0` | F2 (external-vertex) numerators in the coefficient symbols cp and cpp, as before 2026-10-31: quick FORM generation. The default `2` writes them per node in whichever of the two bases FORM makes shorter (cp/cpp, or cH = B0 - cp and cP; `1`: always cH/cP). Same values; g4_24, g4_59, g4_133 1.5-2.2x faster, g5_3086 1.7x and g5_1078 5.2x faster at FORM O1, and O3 then finishes at five loops; FORM generation ~6x slower (`1`: ~3x), with a note on stderr |
| `gen.py L OUTDIR --only "VOLKOV STRING"` | generate selected graphs |
| `gen.py L OUTDIR --range A:B` | generate a slice of the graph list (names stay global) |
| `gen.py L OUTDIR --dirac-only` | tables only |
| `gen.py L OUTDIR --missing` | skip existing files |

The driver's graph registry is built by `make` from the `G_<name>` descriptors in the kernel
files.

## Drawing the diagrams

`python3 generator/diagrams.py diagrams.tex` writes TikZ drawings of every two- and
three-loop diagram, each labelled with its kernel name, Volkov string and weight. Compile with
`pdflatex`.

## Adding a class

Write the kernels into `kernels/<class>/` with a `const graph_t G_<name>` descriptor per file.
Add the class to `ALLCLASSES` in the Makefile, and give it a group in `driver/graphs.c` if useful.
Check finiteness with `bin/sincqed scale <graph>` before trusting any value.
