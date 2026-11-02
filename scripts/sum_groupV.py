#!/usr/bin/env python3
"""Group V total from a run of scripts/run_groupV.sh.

    python3 scripts/sum_groupV.py RESULTS LOG

RESULTS is the sincqed results file, LOG the matching sincqed.log.  For each g4 graph it takes the
value (mirror weight included) and adds the pruning loss left at eps from the graph's "eps extrapolation" (older logs: "eps shadow")
log line.  If a graph appears more than once, the run at the smallest eps is used.  Prints the
number of graphs, the raw sum, the summed shadow correction and the corrected total, against
Laporta's -2.176866027739 (arXiv:1704.06996, sets 1-6)."""
import re, sys
NUM = re.compile(r'^[+-]\d\.\d+e[+-]\d+$')
LAPORTA = -2.176866027739
res = {}
for l in open(sys.argv[1]):
    if l.startswith('#') or not l.strip(): continue
    t = l.split()
    if not t[0].startswith('g4_'): continue
    i = next(i for i in range(4, len(t)) if NUM.match(t[i]))      # the value; eps is t[i-2]
    eps = float(t[i - 2])
    if t[0] not in res or eps <= res[t[0]][0]: res[t[0]] = (eps, float(t[i]))
shadow = {}
for l in open(sys.argv[2]):
    m = re.match(r'# (g4_\d+) eps (?:shadow|extrapolation): .*pruning loss left at eps ~ ([+-][\d.e+-]+)', l)
    if m: shadow[m.group(1)] = float(m.group(2))     # the last line of a graph wins (its final run)
raw = sum(v for e, v in res.values())
sh = sum(shadow.get(g, 0.0) for g in res)
nosh = sorted(g for g in res if g not in shadow)
print('graphs %d (expect 269)%s' % (len(res), '; no shadow line for ' + ' '.join(nosh) if nosh else ''))
print('raw sum      %+.6f' % raw)
print('shadow       %+.6f   (sum |.| %.6f)' % (sh, sum(abs(shadow.get(g, 0.0)) for g in res)))
print('corrected    %+.6f   (Laporta %+.6f, difference %+.6f)' % (raw + sh, LAPORTA, raw + sh - LAPORTA))
