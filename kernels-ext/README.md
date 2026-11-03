# kernels-ext

The four- and five-loop kernels that are too large for git go here, one directory per class:

| directory | files | size | contents |
|---|---|---|---|
| `g4/` | 269 (`g4_00.c` … `g4_268.c`) | ~310 MB | group V: the 518 four-loop diagrams without lepton loops |
| `lbl4/` | 19 (`l4_b.c`, `l4_bp*.c`, `l4_c*.c`, `l4_d1.c`, `l4_d2.c`, `l4_vp.c`) | ~125 MB | four-loop light-by-light |
| `g5/` | 11 (`data/g5/graphs.txt`) | 12-22 MB each | the five-loop graphs of the paper, vectorised, optimised and split |

Everything else in this directory is ignored by git.

Get them from the release archives, which unpack here when run from the top of the repository:

```
tar xJf kernels-g4.tar.xz
tar xJf kernels-lbl4.tar.xz
tar xJf kernels-g5.tar.xz
sh scripts/check_g4.sh          # checks the g4 tree
```

or regenerate them (FORM needed; see [docs/GENERATOR.md](../docs/GENERATOR.md)):

```
sh generator/regen.sh kernels-ext "$(nproc)" g4 lbl4
sh scripts/run5batch.sh data/g5/graphs.txt 2.4 1e-6     # g5: generates, builds and runs h = 2.4
```

`sh scripts/pack_release.sh` makes the archives from this directory.

The Makefile reads `kernels-ext/<class>/`; `make KEXT=DIR` reads them from elsewhere. Without them,
`make` builds the classes that are in git (two and three loops, and the four-loop lepton-loop
classes v4, w4, pi6).
