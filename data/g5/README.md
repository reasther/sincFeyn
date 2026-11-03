# data/g5

The five-loop runs behind Fig. 5 of the preprint, and the reference values they are compared with.

| file | contents |
|---|---|
| `graphs.txt` | the eleven graphs |
| `volkov2019.tsv` | Volkov's values, sigma and V100 GPU-seconds for them (PRD 100, 096004, supplement) |
| `g5_<graph>_h<h>_eps<eps>_top<t>_N<N>.txt` | one `bin/sincqed` results file per run |

`scripts/fit5.py` and `scripts/plot5_cost.py` read this directory by default. They use the top = 0 runs at
the matched eps of each h and at h = 2.4, eps = 1e-5; the other runs (h = 1.8, extra eps values, one
top = 1 run) are earlier tests kept for completeness and do not enter the results.

The runs are heterogeneous: made from 31 October to 3 November 2026 on three machines with successive
driver versions (each file's second line, where present, gives the date and the command). The driver
changes over that period were speed and memory work meant to leave the values unchanged, or changed only
far below the h and eps errors; the times depend strongly on the machine and version, which is why the
cost in Fig. 5 is computed from the evaluation counts. Three h = 2.1 files (g5_1078, g5_372, g5_1221) were rebuilt from
their run logs because the results file was not kept, so they have no header line (nor do two older
files); the time of g5_372 is approximate (the run stalled for about 77 minutes while the machine slept).

The preprint's figure is drawn from these runs; the final version will be redrawn from runs of the
tagged release.
