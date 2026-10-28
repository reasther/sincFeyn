# generator

Python (NumPy, SymPy) that writes the C kernels in `kernels/`. Not needed to build or run.
Start with `sh generator/regen.sh OUTDIR`; the pipeline and options are described in
[../docs/GENERATOR.md](../docs/GENERATOR.md).

| file | role |
|---|---|
| graphs.py | graph enumeration, divergent subgraphs, forests, Volkov's subtraction terms |
| dirac.py | Dirac-structure tables (integer F2 coefficients; vertex and self-energy coefficients) |
| symb.py | Symanzik polynomials and momentum/contraction coefficients (SymPy) |
| build.py | nodes, factored C emission, term algebra, wrapper over the Sinc variables |
| gen.py, gen4all.sh | all graphs without lepton loops at a given order |
| vp.py, vp4.py, vp4x.py, pi6.py | vacuum-polarization insertions (Pi2, Pi4, Pi6) |
| lbl.py, lbl4b.py, lbl4c.py, lbl4d.py, lblvp.py, lblcommon.py | light-by-light classes |
| diagrams.py | TikZ drawings of the two- and three-loop diagrams |
| regen.sh | the whole pipeline, sorted into kernel classes |
