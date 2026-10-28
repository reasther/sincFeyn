# Two-loop Sinc sums in Mathematica

An independent, readable version of what `sincqed` sums at two loops, for checking the C by hand.

- `TwoLoopSinc.wl`: the subtracted integrands of g2_00 ... g2_03 and v2_00 as Wolfram Language
  functions, with the general term of the Sinc sum in three forms: "All" (every Schwinger parameter
  summed, D = 6), "Delta" (Cheng–Wu, D = 5) and "DeltaGL" (Cheng–Wu, the one-point u rule and the
  same-momentum merges, D = 3 or 4: what `sincqed` sums). Exact values in Volkov's scheme are
  included (`exactForm`, `a2Exact` = A_1^{(4)}).
- `TwoLoopSinc.nb`: the notebook. It loads `TwoLoopSinc.wl` from its own directory, checks that the
  three forms agree (homogeneity, the scale integral, the u integral), shows the general terms and
  runs small Sinc sums.
- `check/c2wl.py`: translates `kernels/g2/g2_0N.c` into the `ptg20N` functions of `TwoLoopSinc.wl`
  (the v2_00 kernel is translated by hand). The committed `.wl` is its output for the committed
  kernels:

      python3 mathematica/check/c2wl.py kernels/g2/g2_00.c g2_00 ptg200 kernels/g2/g2_01.c g2_01 ptg201 \
          kernels/g2/g2_02.c g2_02 ptg202 kernels/g2/g2_03.c g2_03 ptg203

- `check/boxsum_gl.c`: full-box DeltaGL sums with the C kernels themselves (no pruning), for the
  reference values in the notebook. From the repository top:

      gcc -O2 -Idriver -I. -DQREAL_DOUBLE -o boxsum_gl mathematica/check/boxsum_gl.c -lm
      ./boxsum_gl 1.0 35

- `check/boxsum.py`, `check/v2_boxsum.py` (with `check/gen_pw.py`, `check/k/`): numpy versions of
  the "All" and "Delta" box sums, which are too long for the notebook
  (`python3 boxsum.py g2_02 delta 1.0 22`).
- `diagrams/`: TikZ drawings of the five two-loop graphs with the line labels x0 ... x5 used in the
  functions (`pdflatex twoloop_diagrams.tex`).
