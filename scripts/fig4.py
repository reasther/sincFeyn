#!/usr/bin/env python3
"""Fig. 4 of the paper: relative error of the three-loop total (all 72 diagrams, fixed and matched eps) and of
four-loop group V (518 diagrams, matched eps) against cost, from the runs in data/g3 and data/g4.

    python3 scripts/fig4.py [OUT.tsv]          (default fig4.tsv; from the top of the tree)

For each h the total is the sum over graphs of value + eps_corr (the pruning loss left at eps), compared with
A1^(6) = 1.18124145658720 (Laporta-Remiddi) and group V = -2.176866027739 (Laporta, sets 1-6).  Cost is the
summed wall time times 1.05 TFLOP/s, the nominal FP64 peak taken for the desktop (24-thread Apple M2 Ultra)
the runs were made on.  Columns: set, h, eps (for group V the 1e-10 graphs; the other 70 use a tenth of it),
graphs, corrected total, relative error, hours, TFLOP-h.
"""
import os, sys
from collections import defaultdict

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
PEAK = 1.05
EXACT = {'g3': 1.18124145658720, 'g4': -2.176866027739}

def rows(path):
    head = None
    for l in open(path):
        if l.startswith('#') or not l.strip(): continue
        t = l.rstrip('\n').split('\t')
        if head is None: head = t; continue
        yield dict(zip(head, t))

def totals(path, label, exact, out):
    tot = defaultdict(lambda: [0, 0.0, 0.0, None])        # h -> graphs, corrected, seconds, eps
    for r in rows(path):
        a = tot[float(r['h'])]
        a[0] += 1; a[1] += float(r['corrected']); a[2] += float(r['seconds'])
        e = float(r['eps']); a[3] = e if a[3] is None else max(a[3], e)
    for h in sorted(tot, reverse=True):
        n, s, sec, e = tot[h]
        out.append((label, h, e, n, s, abs(s - exact) / abs(exact), sec / 3600, sec / 3600 * PEAK))

def main():
    out_path = sys.argv[1] if len(sys.argv) > 1 else 'fig4.tsv'
    out = []
    totals(os.path.join(ROOT, 'data', 'g3', 'runs3_fixed.tsv'), '3-loop fixed eps', EXACT['g3'], out)
    totals(os.path.join(ROOT, 'data', 'g3', 'runs3_matched.tsv'), '3-loop matched eps', EXACT['g3'], out)
    totals(os.path.join(ROOT, 'data', 'g4', 'groupV_matched.tsv'), '4-loop group V matched eps', EXACT['g4'], out)
    with open(out_path, 'w') as f:
        f.write('set\th\teps\tgraphs\tcorrected_total\trel_error\thours\tTFLOP_h\n')
        for r in out:
            f.write('%s\t%.2f\t%.2g\t%d\t%.12f\t%.3e\t%.4f\t%.4f\n' % r)
    for r in out:
        print('%-28s h %.2f  eps %.2g  %3d graphs  %+.10f  rel %.2e  %7.3f h  %.4f TFLOP-h' % r)
    print('wrote', out_path)

if __name__ == '__main__':
    main()
