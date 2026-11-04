#!/usr/bin/env python3
"""The five-loop results in the paper: h extrapolation plus pruning correction, compared with Volkov (2019).

    python3 scripts/fit5.py [RESULTS_DIR ...]          (default data/g5; from the top of the tree)
    K13=0.49 python3 scripts/fit5.py                   (other eps-rule coefficient for D = 13 graphs)

Reads every top-0 result file g5_*_top0_*.txt (one line per run, as bin/sincqed writes them; the runs
come from scripts/run5_hscan.sh) and, for each graph:

  1. takes the value at each h at the matched threshold eps(h) = 1e-6 (h/2.4)^D, which keeps the pruned
     fraction about the same at every h.  A run within 5% of eps(h) is used directly; otherwise the value
     is interpolated in log eps between the two runs that bracket it; an h with neither is left out.
  2. fits S(h) = S_inf + A exp(-c/h) by least squares over all of those h (c on a grid of 0.01 in
     [0.5, 60], S_inf and A linear).  With exactly three points the fit is exact and ill-conditioned, so
     four or more h are wanted.
  3. adds the pruning correction delta = k (S(2.4, 1e-6) - S(2.4, 1e-5)), k = r/(1-r) with r = 0.24 the
     ratio of successive eps-decade steps measured on D = 12 graphs (Aitken on eps 1e-5, 1e-6, 1e-7 at
     h = 2.4 agrees to ~5e-4).  K12 and K13 override k for D = 12 and D = 13.
  4. compares with Volkov's value (data/g5/volkov2019.tsv), in absolute terms, percent and his sigma
     (the side of his asymmetric error that faces our value).
"""
import glob, math, os, re, sys
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
NUM = re.compile(r'^[+-]\d\.\d+e[+-]\d+$')
R = 0.24
K = {12: float(os.environ.get('K12', R / (1 - R))), 13: float(os.environ.get('K13', R / (1 - R)))}
CLASS = ['g5_3086', 'g5_1784', 'g5_372', 'g5_2490', 'g5_1078']      # class abcdeeabcd, Volkov sum -0.240488

def read(dirs):
    runs, Ds = {}, {}                        # (graph, h) -> {eps: (N, value)}
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, 'g5_*_top0_*.txt'))):
            for l in open(f, errors='replace'):
                t = l.split()
                if not t or not t[0].startswith('g5_'): continue
                i = next((k for k in range(4, len(t)) if NUM.match(t[k]) and len(t[k]) > 18), None)
                if i is None: continue
                h, N, eps, D = round(float(t[i - 4]), 3), int(t[i - 3]), float(t[i - 2]), int(t[i - 6])
                old = runs.setdefault((t[0], h), {}).get(eps)
                if old is None or N > old[0]: runs[(t[0], h)][eps] = (N, float(t[i]))
                Ds[t[0]] = D
    return runs, Ds

def matched(ev, tgt):
    es = sorted(ev)
    ex = [e for e in es if abs(e / tgt - 1) < 0.05]
    if ex: return ev[ex[0]][1]
    lo = [e for e in es if e < tgt]; hi = [e for e in es if e > tgt]
    if not (lo and hi): return None
    e1, e2 = lo[-1], hi[0]; x = math.log(e2 / tgt) / math.log(e2 / e1)
    return ev[e2][1] + (ev[e1][1] - ev[e2][1]) * x

def lsq(h, S):
    h, S = np.array(h), np.array(S); best = None
    for c in np.arange(0.5, 60, 0.01):
        X = np.c_[np.ones_like(h), np.exp(-c / h)]; co = np.linalg.lstsq(X, S, rcond=None)[0]
        r = float(np.sum((X @ co - S) ** 2))
        if best is None or r < best[0]: best = (r, c, co[0], co[1])
    return best[1], best[2], math.sqrt(best[0] / len(h))         # c, S_inf, rms residual

def volkov():
    v = {}
    for l in open(os.path.join(ROOT, 'data', 'g5', 'volkov2019.tsv')):
        if l.startswith('#'): continue
        t = l.split('\t'); v[t[0]] = (float(t[4]), float(t[5]), float(t[6]))
    return v

def main():
    dirs = sys.argv[1:] or [os.path.join(ROOT, 'data', 'g5')]
    runs, Ds = read(dirs); vk = volkov()
    graphs = sorted(Ds, key=lambda g: (Ds[g], int(g[3:])))
    print('%-8s %2s %-24s %6s %11s %10s %11s %11s %10s %8s %6s' % (
        'graph', 'D', 'h used', 'c', 'S_inf', 'eps corr', 'result', 'Volkov', 'diff', 'percent', 'sigma'))
    out = {}
    for g in graphs:
        D = Ds[g]; pts = []
        for (gg, h), ev in sorted(runs.items(), key=lambda x: -x[0][1]):
            if gg != g: continue
            v = matched(ev, 1e-6 * (h / 2.4) ** D)
            if v is not None: pts.append((h, v))
        if len(pts) < 3: print('%-8s %2d  only %d matched h' % (g, D, len(pts))); continue
        c, Si, rms = lsq([p[0] for p in pts], [p[1] for p in pts])
        ev = runs.get((g, 2.4), {})
        s6 = next((ev[e][1] for e in ev if abs(e / 1e-6 - 1) < 0.05), None)
        s5 = next((ev[e][1] for e in ev if abs(e / 1e-5 - 1) < 0.05), None)
        de = K[D] * (s6 - s5) if s5 is not None and s6 is not None else 0.0
        res = Si + de; out[g] = res
        hs = ','.join('%g' % p[0] for p in pts)
        row = '%-8s %2d %-24s %6.1f %+11.6f %+10.2e %+11.6f' % (g, D, hs, c, Si, de, res)
        if g in vk:
            V, up, dn = vk[g]; d = res - V; sg = up if d > 0 else dn
            row += ' %+11.6f %+10.2e %+7.3f%% %+6.1f' % (V, d, 100 * d / abs(V), d / sg)
        if s5 is None: row += '  (no eps 1e-5 run at h 2.4: no eps correction)'
        print(row)
    rel = [abs(out[g] - vk[g][0]) / abs(vk[g][0]) for g in out if g in vk]
    if rel: print('\nmean |relative difference| over %d graphs: %.2e' % (len(rel), sum(rel) / len(rel)))
    if all(g in out for g in CLASS):
        s = sum(out[g] for g in CLASS); print('class abcdeeabcd: %+.6f  (Volkov -0.240488, diff %+.2e)' % (s, s + 0.240488))
    print('eps rule k: D12 %.4f  D13 %.4f' % (K[12], K[13]))

if __name__ == '__main__':
    main()
