# Four-loop lepton-loop classes

Exact values: Laporta, arXiv:1704.06996, Table 2 (sets 1-25), matched to Kinoshita's groups.

| group | diagrams | Laporta | generator | entries |
|---|---|---|---|---|
| V no lepton loops | 518 | sets 1-6: -2.176866027739 (sum) | gen.py / gen4all.sh | g4_* (269) |
| I(a) Pi2^3 | 1 | 0.000876865858890 | vp.py 4 | v4_00 |
| I(b) Pi2+Pi4 | 6 | 0.015325282902013 | vp4x.py | w4_00* |
| I(c) Pi4 with Pi2 inside | 3 | 0.011130913987517 | vp4x.py | w4_0i* |
| I(d) Pi6 | 15 | 0.049513202559526 | pi6.py | w4_p6* (27 files, 3 zero) |
| II(a) Pi4 in 2-loop vertices | 36 | -0.420482262508962 | vp4x.py | w4_01..07 |
| II(b)+(c), III Pi2 insertions | 168 | -0.027675, -0.073452, 1.417687 | vp.py 4 | v4_* |
| IV(a) LbL + Pi2 | 18 | 0.598842072031422 | lblvp.py | l4_vp |
| IV(b) photon on the loop | 60 | 0.822284485811034 | lbl4b.py | l4_b (+ parts l4_bp*) |
| IV(c) photon on the open line | 48 | -1.138822876459975 | lbl4c.py | l4_c02..13 |
| IV(d) loop with 4 photons, mu on line | 18 | -0.990607260864552 | lbl4d.py | l4_d1, l4_d2 |

Diagram counts: legs excluded, both loop orientations counted (as the 6 three-loop LbL graphs).
IV(c) = 8 placements x 6, IV(b) = 5!/2 loop sequences, IV(d) = 3 mu positions x 6.  The IV(b)/(c)/(d)
assignment to Laporta's sets 21, 23, 24+25 is by value match to Kinoshita's numbers.

## Method common to the LbL classes
- Each ordering of the photons on the loop is evaluated minus the same loop at zero photon momenta
  (lbl.py `zero=True`), the orderings share one parameter labelling and are summed pointwise.
- Terms follow Volkov's formula with G' in I[G] (vertex-like subgraphs containing mu).
  IV(c) placement (0,4) has the three-loop LbL vertex H as such a subgraph:
  A(G) - A(G/H) U(H) - L(G/H) A(H) + U(G/H) A(H).
- Dirac tables of the loop nodes live in generator/cache_lbl.pkl (md5 of chain+loop tokens).
  8-slot F2 tables take tens of minutes each.

## Light-by-light vertex subgraphs
A light-by-light loop with three photons on one contiguous stretch of the open line forms a
three-loop LbL VERTEX subgraph around the fourth photon (log divergent, needs U).  IV(d) mu at 2
has one (scale test: 2 slow directions, ratio 1).  Check every new class with `bin/sincqed scale`.

## Pi6
pi6.py: general single-loop VP with n internal photons in the one-loop vertex (n=1 reproduces vp4.py
pointwise and integrated).  Forests of U over cyclic intervals; inner nodes are build.Node on the loop
cut at the external vertex outside the interval.  Families: C0 main (dp n-2L-w), C1 contact (n-2L-1-w),
C2 contact with -d/dt of the parent (n-2L-2-w; all C_m scale as lambda^(-2L-w)).
w4_p63d, p65d, p67d are identically zero (the driver detects a zero integrand and returns 0 at once).
p63 and p65 give identical values.
