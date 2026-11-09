# High-precision two-loop Sinc runs for the PSLQ closed forms

This package reruns the Sinc sums that gave the exact per-graph values of the two-loop g-2 graphs in
Volkov's subtraction scheme (mirror graphs included): the analytic forms in Table I of the paper.

| graph | closed form | value | Sinc run used for PSLQ |
|---|---|---|---|
| ladder (g2_03) | 11/48 + π²/18 | 0.77747802228274215 | unreduced 5-D sum, h = 0.5, eps = 1e-20 (residual 4e-15) |
| corner (g2_00) | −55/24 + π²/18 + π² log 2/3 − ζ(3)/2 | −0.06402094134468282 | unreduced 5-D sum, h = 0.5, eps = 1e-20 (residual 8e-15) |
| self-energy (g2_01) | −1/24 − π²/18 | −0.58997802228274215 | reduced 3-D sum in ~64-digit arithmetic, h = 0.4, eps = 1e-26 (residual 1e-17) |
| crossed (g2_02) | 1/6 + 13π²/36 − 5π² log 2/6 + 5ζ(3)/4 | −0.46764544609361365 | A_1^{(4)} minus the others; matches the literature −0.467645446094 |
| Π2 insertion (v2_00) | 119/36 − π²/3 | 0.01568742185910 | known |

A plain build needs no Python. The PSLQ step itself is optional Python (mpmath).

    make -j8                    # build (g++ with OpenMP; on macOS the newest Homebrew gcc, e.g. brew install gcc)
    make test                   # the whole chain at coarse steps, checked against the closed forms (a few minutes)
    THREADS=24 sh run_pslq.sh   # the runs used for the paper -> pslq_runs/values.tsv
    python3 pslq/identify.py pslq_runs/values.tsv    # optional: PSLQ on the values (mpmath 1.3.0)

## What `run_pslq.sh` runs

1. `demo/bin/sincqed demo 0.5 1e-20 all2`, with `G2_ABSERR=1e-18`. This is the unreduced Sinc sum over
   all five free Schwinger parameters of each two-loop graph: ladder, corner, crossed, self-energy and Π2.
   It runs in double precision, and a point is recomputed in double-double when its rounding bound
   exceeds 1e-18. It is the `sincqed demo` mode of an earlier sincFeyn driver, kept here because the current driver has
   no such mode; see `demo/README.md`.
   The ladder and corner values come from here (about 14 digits), and the sum of all graphs is checked
   against A_1^{(4)}.
2. `quad/se_qd 0.4 1e-26`. This is the self-energy as a reduced 3-D sum in quad-double arithmetic,
   about 64 digits. Its two electron lines in series are merged into one class, and the external-vertex
   split is done with the exact one-point rule (below). It gives about 17 digits. h = 0.3 gives about 21, at
   roughly 2.5 times the cost.

Extended precision is needed for the self-energy because the kernel cancels digits at large Schwinger
parameters, and its |terms| rounding bound does not see that. `__float128` alone stops at about 1e-16.

Settings are environment variables: `H_DEMO`, `EPS_DEMO`, `H_SE`, `EPS_SE`, `OUT` (default `pslq_runs`) and
`THREADS`. The script can be restarted: a step whose output is complete is skipped.

`reference/` holds the outputs of the full run here, including `identify.txt`.

`values.tsv` has one row per graph. Its columns are the graph, the Sinc value, the estimated error and
the run. The error estimate is the measured h convergence of each sum: about 1e-6 e^{−18(1/h−1)} for the
unreduced sums, and about 3e-11 e^{−17.6(1/h−5/3)} for the reduced self-energy. The quad-double value is
written as its four exact double components, which `identify.py` sums exactly. `identify.py` runs PSLQ
on the basis {1, π², π² log 2, ζ(3)} of A_1^{(4)}, with `--basis large` adding log 2, log² 2 and π⁴. It
uses a tolerance of 10 times the estimated error. It flags any relation whose residual is larger than
that, as happens for the coarse `QUICK=1` values.

## Measured here (Linux, 4 cores)

| run | value | − closed form | evals | time |
|---|---|---|---|---|
| ladder, demo h=0.5 eps=1e-20 | +0.77747802228273855 | −3.6e-15 | 1.9e9 | 488 s |
| corner, demo h=0.5 eps=1e-20 | −0.06402094134469101 | −8.2e-15 | 2.0e9 | 277 s |
| crossed, demo h=0.5 eps=1e-20 | −0.46764544609359909 | +1.5e-14 | 2.6e9 | 226 s |
| self-energy, demo h=0.5 eps=1e-20 | −0.58997802228274199 | +1.6e-16 | 8.5e8 | 164 s |
| Π2, demo h=0.5 eps=1e-20 | +0.015687421859138045 | +3.5e-14 | 4.5e8 | 62 s |
| sum of the five vs A_1^{(4)} | −0.32847896557915546 | +3.8e-14 | | |
| self-energy, se_qd h=0.4 eps=1e-26 | −0.58997802228274213478 | +1.1e-17 | 2.7e6 | 138 s |

The whole `run_pslq.sh` takes about 23 minutes, with a peak memory of a few GB. With the small basis, PSLQ
(`identify.py`) returns the closed forms above for every graph. The values are the same as the
original runs to every printed digit.

## Other programs in `quad/`

- `corner_qd`, `crossed_qd` (`h eps [zmax]`): the 4-D reduced corner and crossed sums in quad-double.
  They need h ≈ 0.4 for 16 digits, which is about 1e9 points, so they were not the route to the closed
  forms. They are kept for cross-checks.
- `make quadmath` (x86-64 with libquadmath only): the original programs. `se_qq` uses double-quad (two
  `__float128`), and `corner_q`, `crossed_q` use `__float128`. `se_qd` agrees with `se_qq` to 1e-37 at
  h = 1 and runs 2.3 times faster.

`quad/qd.h` is the quad-double type. Each operation is an exact floating-point expansion (two_sum, with
products split by fma) that is compressed and cut to four doubles, so it is correct to about 2^-208. It
is built with `-ffp-contract=off`. The kernels (`g2_00p.c`, `g2_01.c`, `g2_02p.c`) are the generated
sincqed kernels.

**The external-vertex split.** At zero momentum transfer each kernel's point function is linear in the
split u of the electron line at the external vertex, times a weight (1−u)^p u^q. One point at the
weighted mean, u = (q+1)/(p+q+2) with weight p! q!/((p+q+1)! (1−u)^p u^q), integrates it exactly (the
rule of the main sincFeyn kernels, `tools/glpts.py`). The kernels here evaluate that single point
directly, with u and the weight written as exact ratios in the working precision: u = 1/2 and weight 1
for p = q = 0 (corner, crossed, ladder, Π2, one loop) and u = 1/3, weight 3/4 for the self-energy
(p = 1). There is no Gauss-Legendre rule anywhere in this package; an earlier version used the 2-point
rule, which gives the same values to about 1e-65 at two to three times the cost. `sincqed demo` does
not use the split at all, since it sums over every Schwinger parameter.
