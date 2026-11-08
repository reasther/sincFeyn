# Data behind the paper's figures

| directory | figure | contents |
|---|---|---|
| `g3/` | Fig. 4, three loops | `runs3_fixed.tsv`, `runs3_matched.tsv`: every three-loop entry (all3, 41 integrals = 72 diagrams) at h = 2.10, 2.05, ..., 1.00, with eps = 1e-15 (fixed) or 1e-11 (h/1.5)^10 (matched) |
| `g4/` | Fig. 4, four loops | `groupV_matched.tsv`: the 269 group V integrals (518 diagrams) at h = 2.2, 2.1, 1.9, 1.7, 1.5, 1.3 with eps = 1e-10 or 1e-11 times (h/1.5)^10 |
| `g5/` | Fig. 5 | the eleven five-loop graphs; see [g5/README.md](g5/README.md) |

One row per graph and h. `value` is the Sinc sum (mirror weight included), `eps_corr` the pruning loss left at
eps that the driver's `eps extrapolation` log line reports, `corrected` their sum, `seconds` the wall time on
the desktop (24-thread Apple M2 Ultra). `python3 scripts/fig4.py` turns `g3/` and `g4/` into the Fig. 4 table
(totals against A1^(6) and Laporta's group V value, and cost in TFLOP-h); `python3 scripts/fit5.py` and
`scripts/plot5_cost.py` do the same for `g5/`.

These runs were made with the code as it stood between 3 and 7 November 2026 (see
[docs/REPRODUCE.md](../docs/REPRODUCE.md) for the commands); the preprint's figures are drawn from them.
