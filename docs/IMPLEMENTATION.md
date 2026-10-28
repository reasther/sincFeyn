# How sincqed computes a g−2 graph

This note follows one graph from its integral to the number in the results line. The kernel
files it refers to are generated; [GENERATOR.md](GENERATOR.md) describes how.

## 1. The integral

Each graph contributes a coefficient of (α/π)^L. After the loop momenta are integrated out
(they are Gaussian in Schwinger form), it is an integral over the Schwinger parameters
s_i ∈ (0, ∞) of its internal lines. Three exact reductions lower the number of variables
before any sum is taken.

- **Overall scale (Cheng–Wu).** The integrand is projective, so one parameter can be fixed to 1.
  The generator fixes a photon line: the photon at the first arc, or, in a graph that is its
  own mirror image, a photon that is its own mirror, so that the mirror becomes a permutation of
  the remaining axes.
- **Lines with the same momentum.** At zero momentum transfer, lines that carry the same loop
  momentum (a series chain, including the two lines at the external vertex) enter the
  Symanzik polynomials only through the sum of their parameters. Each such class of K lines
  becomes one variable σ, with weight σ^(K−1)/(K−1)!.
- **The split at the external vertex.** The integrand is linear in how σ is split between the
  two lines at the external vertex, u = s_b/σ, times a simplex weight (1−u)^p u^q. One point at
  the weighted mean u = (q+1)/(p+q+2) integrates it exactly (`f_<graph>` in the kernel).

What is left is D summation variables: at most 3L − 2 at L loops (D = 10 for the largest
four-loop graphs).

**Example (one loop).** With electron parameters s_A, s_B and the photon set to 1,

    a1 = ∫∫ ds_A ds_B  1 / ((s_A+s_B)(1+s_A+s_B)^3) = ∫ dσ 1/(1+σ)^3 = 1/2.

## 2. Subtractions: Volkov's forest formula, pointwise

The ultraviolet and infrared divergences are removed before integration, point by point in
parameter space, with Volkov's forest formula (arXiv:1507.06435, 1705.05800). Nothing divergent
is ever integrated, and no difference of large numbers is taken after the sum.

- `generator/graphs.py` finds the UV-divergent subgraphs (vertex and self-energy subgraphs),
  the forests, and the terms of the formula, each a product of Volkov's operators A
  (anomalous-magnetic-moment projection), U (UV subtraction) and L (the vertex
  renormalisation constant) applied to nested subgraphs.
- Each forest element becomes a **node**: the U, V and W polynomials of a subgraph, with the
  Dirac numerator reduced to a polynomial in the node's momentum and contraction
  coefficients. The integrand (`pt_<graph>`) is the signed sum over the forest of products of
  W/V^n factors, evaluated at the same Schwinger parameters.
- The price is cancellation between the terms near the edges of the region, which the
  precision fallback (Sec. 4) handles.

**How the Dirac coefficients of the nodes are obtained.** A physicist would reduce each
node's Dirac string algebraically. The generator was built differently, and this is the step
where the code departs most from a hand calculation:

1. For each node it expands the numerator into structures (which slot carries p, p′, the mass,
   or a loop momentum paired with another slot).
2. For each structure it evaluates the Dirac string numerically with explicit 4×4 gamma
   matrices at a set of kinematic points, then projects onto the quantity the operator needs:
   the F2 projector at zero momentum transfer (A), the vertex coefficients (U, L), or the
   self-energy coefficients.
3. It fits the result as a polynomial in t (the momentum transfer), and rounds the coefficients
   to the exact integers or small rationals they must be, with assertions that the fit is exact
   and the rounding is within 1e-6 (`generator/dirac.py`).

The resulting integer tables are exact, but they are not derived symbolically. As an
independent check, `generator/formdirac.py` computes the same tables algebraically with FORM
(symbolic traces). FORM agrees exactly with every two- and three-loop table and with the
largest four-loop table checked, and is the default for the classes without lepton loops (g2,
g3, g4, v3, v4) when FORM is installed. The lepton-loop classes (g3_lbl, pi6, v3_p4, w3, w4,
lbl4) still take their tables from the numerical route.

The other independent check of the subtracted integrand is the **finiteness (scaling) test**,
`sincqed scale GRAPH`: it scales the parameters along every direction in {−1, 0, 1}^(n−1) and
checks that the integrand decays fast enough in all of them. Every two- and three-loop graph
and the four-loop light-by-light classes pass. A missing subtraction shows up as a direction of
slow decay.

## 3. From an integral to a Sinc sum

Set s = e^z. Every variable then runs over all of ℝ, and ds = s dz. The trapezoid ("Sinc")
rule on a lattice of step h is

    I ≈ Σ_k  f(e^{z(k)}) · Π_i s_i(k) · h^D .

For an integrand that is analytic in a strip |Im z| < d around each real axis and decays at
both ends, the error is about A·exp(−2π d/h): exponential in 1/h, not a power of h.

**Example.** For the one-loop integrand in σ, the poles of 1/(1+e^z)^3 are at z = ±iπ, so
d = π. The Sinc sum minus 1/2:

| h    | error     | exp(−2π²/h) |
|------|-----------|-------------|
| 2.0  | 1.0e-3    | 5.2e-5      |
| 1.5  | 5.1e-5    | 1.9e-6      |
| 1.0  | 1.1e-7    | 2.7e-9      |
| 0.75 | 2.0e-10   | 3.7e-12     |
| 0.5  | 5.7e-16   | 7.2e-18     |

The error falls at the predicted rate, with a prefactor of about 20–80 from the triple pole.

## 4. The hierarchical sum (`driver/sum.c`)

Lattice node k ∈ ℤ^D, with |k_i| ≤ K, sits at z_i = k_i h; the window is |z| ≤ Kh, about 40–48.
A full box needs (2K+1)^D nodes, which is about 10^18 at D = 10, N = 65. Almost all of them are negligible,
and the sum visits only the cells that matter, level by level.

1. **Coarsest level** (step 2^top h; by default top = 1 for h > 0.75 and 2 below, so this step
   stays at most 3 in z). A support-following scan starts at z = 0 and walks through axis
   neighbours while the terms are large enough. It records the reference scale
   Ref = Σ|terms|.
2. **Keep and dilate.** Cells above the keep threshold (a fixed fraction of eps·Ref) are kept, then dilated to their axis
   neighbours, so the support is not cut at its edge.
3. **Each finer level** (step halves) evaluates the 2^D children of every kept cell and keeps
   those above the threshold.
4. **Results.** S_h and S_2h are the sums at the last two levels; the value reported is S_h.

Sums are accumulated with compensated (Kahan–Neumaier) summation in blocks, so that the many
small terms are not lost to rounding.

**Symmetry fold.** If axis permutations leave the integrand unchanged (the mirror of a
self-mirror graph, the reversal of a lepton loop, equal insertions), the driver finds them,
evaluates only one point per orbit and weights it by the orbit size.

**Precision (`driver/eval.c`).** Every node is evaluated in double first. The kernel also
returns Σ|terms|. If that bound shows the terms may cancel beyond double precision, the node is
recomputed in double-double (32 digits; `make DD=0` uses `__float128`). Typically 0.1–1% of nodes
need it. A node that still cancels beyond double-double is counted as zero and reported in the
log.

**Error columns** of the results file (already multiplied by the mirror weight w):

| column | meaning |
|--------|---------|
| disc  | 2·\|S_2h − S_h\|² / Σ\|terms at h\|; overestimates at coarse h. The realistic error is \|S(h) − S(1.2h)\|/3 (`sincqed compare`). |
| prune | weight dropped below the keep threshold, summed over every level; understates the true loss at four loops, so use the eps shadow |
| trunc | outer shell of the window, extrapolated geometrically from the last two shells |
| round | points that cancel beyond the high-precision type, in quadrature |
| err   | disc + prune + trunc + round |

**The eps shadow.** Each kept point carries a tier saying whether it would still have been kept
with eps raised by 10 or 100, and the final level is summed per tier. From the two differences
the log estimates the pruning loss still left at eps, at no extra cost:

    # g2_02 eps shadow: eps*10 would change the value by +3.968e-08, eps*100 by a further +3.368e-07 (ratio 8.49); pruning loss left at eps ~ -5.30e-09 (add to the value)

**Example (reading a results line).** g2_02 (the crossed ladder) at h = 1.5, N = 65,
eps = 1e-8:

    g2_02    3; 1-4, 2-5   2 4 1 1.500   65  1.0e-08 1 -4.675953149180136e-01 1.438e-03 1.44e-03 9.80e-08 3.67e-13 0.00e+00 -4.47483885e-01  174738  0.5 -4.676454460936136e-01 +5.013e-05

- **Columns:** name, Volkov string, L = 2, D = 4, w = 1, h, N, eps, top, the value, the error
  columns, S_2h, evaluations, time (s), the exact value and the difference.
- **Size of the sum:** 174,738 integrand evaluations instead of 65^4 = 1.8e7.
- **Accuracy:** 5.0e-5 from the exact value, about the size of a single-propagator Sinc error
  at h = 1.5. disc (1.4e-3) is the coarse-h overestimate mentioned above.

## 5. Kernels and the build

- **Graph code.** `kernels/<class>/<graph>.c` is generated C. Each file has a double and a
  high-precision version of `pt_` (the subtracted integrand at given s) and `f_` (the wrapper
  over the summation variables), and a `G_<graph>` descriptor that registers it with the driver.
- **Classes.** g2 (two loops), g3 (three loops without lepton loops, one file per mirror pair,
  plus g3_lbl for light-by-light), v3 (three-loop vacuum polarisation), and the four-loop classes
  g4, v4, w4, pi6, lbl4 ([USAGE.md](USAGE.md) has the table).
- **Build.** `make CLASSES="g2 g3 v3"` builds a three-loop executable; `make test` checks it.
