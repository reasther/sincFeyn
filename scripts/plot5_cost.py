#!/usr/bin/env python3
"""Fig. 5 of the paper: cost against relative accuracy for the eleven five-loop graphs, sinc against Volkov's
Monte Carlo (PRD 100, 096004).

    python3 scripts/plot5_cost.py [RUNSDIR [OUTPREFIX]]     (defaults data/g5 and g5_cost; from the top of the tree)

Writes OUTPREFIX.png/.pdf (per graph), OUTPREFIX_mean.png/.pdf (per-diagram means: the paper's figure) and
OUTPREFIX.json, and prints the tables behind them.  Projected sinc values are the recipe of scripts/fit5.py:
least squares S_inf + A exp(-c/h) over every matched-eps h, plus the eps correction from one h = 2.4, eps = 1e-5
run; every run that enters is charged.
Our flops: evals x flops/point.  flops/point = 163k for g5_3086 (count of the kernel's double path), scaled for
every other graph by its time per evaluation relative to g5_3086 in the same h = 2.3 batch; machine time is
counted flops / 0.187 TFLOP/s (measured, Mac Studio M2 Ultra) times its nominal FP64 peak SINCQED_MACPEAK [1.05].
Volkov: Calc_time (V100 GPU-seconds, data/g5/volkov2019.tsv) x 7.8 TFLOP/s (V100 FP64 peak)."""
import glob, math, os, re, sys, json
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
NUM = re.compile(r'^[+-]\d\.\d+e[+-]\d+$')
runsdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'data', 'g5')
out = sys.argv[2] if len(sys.argv) > 2 else 'g5_cost'
G = ['g5_3086', 'g5_1784', 'g5_1078', 'g5_372', 'g5_2490', 'g5_1663', 'g5_3024', 'g5_717', 'g5_2567', 'g5_1795', 'g5_1221']
D = {}
runs = {}                       # g -> {(h, eps): (N, value, evals, time)}; largest N per (h, eps)
for f in glob.glob(os.path.join(runsdir, 'g5_*top0*.txt')):
    for l in open(f):
        t = l.split()
        if not t or t[0] not in G: continue
        i = next(k for k in range(4, len(t)) if NUM.match(t[k]) and len(t[k]) > 18)
        D[t[0]] = int(t[i-6]); key = (round(float(t[i-4]), 3), float(t[i-2])); rec = (int(t[i-3]), float(t[i]), int(t[i+7]), float(t[i+8]))
        if key not in runs.setdefault(t[0], {}) or rec[0] > runs[t[0]][key][0]: runs[t[0]][key] = rec
# Volkov: value, mean of his asymmetric sigma, GPU-seconds
VK = {}
for l in open(os.path.join(ROOT, 'data', 'g5', 'volkov2019.tsv')):
    if l.startswith('#'): continue
    t = l.split('\t'); VK[t[0]] = (float(t[4]), (float(t[5]) + float(t[6]))/2, float(t[7]))
def volkov(g): return VK[g]
def tpe(g):                     # time per evaluation in the h = 2.3 batch
    (k, r), = [(k, r) for k, r in runs[g].items() if k[0] == 2.3]
    return r[3]/r[2]
F = {g: 163e3*tpe(g)/tpe('g5_3086') for g in G}
gfh = lambda g, r: r[2]*F[g]/1e9/3600
def fitls(pts):                 # least squares S = S_inf + A exp(-c/h) over all h, c on a grid
    hs = [h for h, _ in pts]; ss = [v for _, v in pts]; n = len(hs); best = None
    for i in range(50, 6000):
        c = i/100; x = [math.exp(-c/h) for h in hs]; mx = sum(x)/n; ms = sum(ss)/n
        sxx = sum((a - mx)**2 for a in x); A = sum((a - mx)*(b - ms) for a, b in zip(x, ss))/sxx; S0 = ms - A*mx
        r = sum((S0 + A*a - b)**2 for a, b in zip(x, ss))
        if best is None or r < best[0]: best = (r, c, S0)
    return best[2]
HS = (2.4, 2.3, 2.2, 2.1, 2.0)    # every matched-eps h point we have; each one is charged
res = []
for g in G:
    V, sig, gs = volkov(g); R = runs[g]; used = set(); pts = []; hc = {}
    for h in HS:
        tg = 1e-6*(h/2.4)**D[g]; es = sorted(e for (hh, e) in R if hh == h)
        if not es: continue
        ex = [e for e in es if abs(e/tg - 1) < 0.05]
        if not ex and not (min(es) < tg < max(es)): continue
        if ex: pts.append((h, R[(h, ex[0])][1])); used.add((h, ex[0])); hc[h] = gfh(g, R[(h, ex[0])]); continue
        e1 = max(e for e in es if e < tg); e2 = min(e for e in es if e > tg); x = math.log(e2/tg)/math.log(e2/e1)
        pts.append((h, R[(h, e2)][1] + (R[(h, e1)][1] - R[(h, e2)][1])*x)); used |= {(h, e1), (h, e2)}
        hc[h] = math.exp(math.log(gfh(g, R[(h, e2)])) + (math.log(gfh(g, R[(h, e1)])) - math.log(gfh(g, R[(h, e2)])))*x)   # log-interpolated cost
    # production recipe: one matched-eps run per h plus one h = 2.4 run at eps 1e-5; the eps correction is geometric
    # with the measured ratio r ~ 0.24 per decade, delta = r/(1-r) (S(1e-6) - S(1e-5)) at h = 2.4, applied to S_inf
    Sinf = fitls(pts); RR = 0.24
    delta = RR/(1 - RR)*(R[(2.4, 1e-6)][1] - R[(2.4, 1e-5)][1]); epsfix = 'cheap'; hc['eps1e-5'] = gfh(g, R[(2.4, 1e-5)])
    proj = Sinf + delta
    # single run: the lowest eps among the scan h values (ties: the smaller h)
    (hs, es), r = min(((k, r) for k, r in R.items() if k[0] in HS), key=lambda kr: (kr[0][1], kr[0][0]))
    res.append(dict(g=g, D=D[g], V=V, sig=sig, vk_gfh=gs/3600*7800, F=F[g],
                    single=dict(h=hs, eps=es, S=r[1], err=abs(r[1] - V)/abs(V), gfh=gfh(g, r)),
                    fixh={str(h): dict(S=v, err=abs(v - V)/abs(V), gfh=hc[h]) for h, v in pts}, proj=dict(S=proj, Sinf=Sinf, delta=delta, err=abs(proj - V)/abs(V), gfh=sum(hc.values()), eps=epsfix, nruns=len(used), hs=[h for h, _ in pts])))
json.dump(res, open(out + '.json', 'w'), indent=1)
print('%-8s %3s %7s %10s %9s %9s | %9s %9s | %9s %9s %s' % ('graph', 'D', 'F/pt', 'Volkov', 'sig_rel', 'MC GFh', 'single', 'GFh', 'proj', 'GFh', 'eps'))
for r in res:
    print('%-8s %3d %6.0fk %+10.6f %9.2e %9.0f | %9.2e %9.1f | %9.2e %9.1f %s' % (r['g'], r['D'], r['F']/1e3, r['V'], r['sig']/abs(r['V']), r['vk_gfh'],
          r['single']['err'], r['single']['gfh'], r['proj']['err'], r['proj']['gfh'], r['proj']['eps'] or 'h'))
# MC fit: cost = A * sigma_rel^-2
# nominal-peak convention (as for g3/g4): machine-hours x nominal FP64 peak for both.  Our counted flops ran at
# RATE TFLOP/s sustained (g5_3086 h2.3, current binary), so Mac-hours = flops/RATE; x MACPEAK (Mac Studio M2 Ultra CPU
# FP64 nominal: 16 P-cores x 16 flop/cycle x 3.5 GHz + 8 E-cores ~ 1.05 TFLOP/s; SINCQED_MACPEAK to change)
RATE = 0.187; MACPEAK = float(os.environ.get('SINCQED_MACPEAK', '1.05')); SC = MACPEAK/RATE
for r in res:
    r['vk_gfh'] /= 1e3; r['single']['gfh'] *= SC/1e3; r['proj']['gfh'] *= SC/1e3
    for v in r['fixh'].values(): v['gfh'] *= SC/1e3
A = math.exp(sum(math.log(r['vk_gfh']*(r['sig']/abs(r['V']))**2) for r in res)/len(res))
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': '#999', 'axes.labelcolor': '#333', 'xtick.color': '#555', 'ytick.color': '#555'})
fig, ax = plt.subplots(figsize=(7.5, 5.2), dpi=150)
B, O, Q = '#2a78d6', '#eb6834', '#1baf7a'
ax.grid(True, which='major', color='#e6e6e6', lw=0.6); ax.set_axisbelow(True)
cs = [1e-2, 300]; ax.plot(cs, [(A/c)**0.5 for c in cs], color='#888', lw=1.5, ls='--', label=r'Monte Carlo fit, $\sigma \propto$ cost$^{-1/2}$')
ax.scatter([r['vk_gfh'] for r in res], [r['sig']/abs(r['V']) for r in res], s=36, color='#888', marker='s', label='Volkov 2019 (V100)', zorder=3, edgecolor='white', lw=0.8)
ax.scatter([r['single']['gfh'] for r in res], [r['single']['err'] for r in res], s=36, color=O, marker='o', label='sinc, single run (lowest eps)', zorder=3, edgecolor='white', lw=0.8)
for m, lab, sel in (('D', 'sinc, projected: 4-5 h runs (2.4-2.1, +2.0 where run) + one 1e-5 run, all charged', 'cheap'),):
    rr = [r for r in res if r['proj']['eps'] == sel]
    ax.scatter([r['proj']['gfh'] for r in rr], [r['proj']['err'] for r in rr], s=150 if m == '*' else 40, color=B, marker=m, label=lab, zorder=4, edgecolor='white', lw=0.8)
for r in res:
    ax.plot([r['single']['gfh'], r['proj']['gfh']], [r['single']['err'], r['proj']['err']], color=B, lw=0.6, alpha=0.35, zorder=2)
    ax.annotate(r['g'][3:], (r['proj']['gfh'], r['proj']['err']), xytext=(4, 3), textcoords='offset points', fontsize=7, color='#555')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(1e-4, 0.5); ax.set_xlim(1e-2, 300)
ax.set_ylabel('relative error  |S - Volkov| / |Volkov|   (Volkov: his own sigma)'); ax.set_xlabel('cost (TFLOP-hours)')
ax.set_title('Five-loop g5 graphs: accuracy vs cost', fontsize=11, loc='left', color='#222')
ax.legend(fontsize=8, frameon=False, loc='upper right', borderaxespad=0.8)
fig.text(0.01, 0.005, 'Both at nominal FP64 peak: sinc Mac Studio hours x %.2f TFLOP/s, Volkov V100 hours x 7.8 TFLOP/s. Errors below ~1e-3 are within Volkov\'s own sigma.' % MACPEAK, fontsize=6.5, color='#777')
fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig(out + '.png'); fig.savefig(out + '.pdf')
print('MC fit A = %.3g TFLOP-h at sigma_rel = 1' % A)

# per-diagram averages: each fixed-h set (matched eps), the projected set, and Volkov
n = len(res); mean = lambda xs: sum(xs)/len(xs)
sets = [('h = %s' % h, mean([r['fixh'][h]['gfh'] for r in res]), mean([r['fixh'][h]['err'] for r in res])) for h in ('2.4', '2.3', '2.2', '2.1')]
pj = (mean([r['proj']['gfh'] for r in res]), mean([r['proj']['err'] for r in res]))
pj3 = [r for r in res if r['proj']['eps']]
vk = (mean([r['vk_gfh'] for r in res]), mean([r['sig']/abs(r['V']) for r in res]))
print('\nper-diagram mean over %d graphs:  set  mean TFLOP-h  mean |rel err|' % n)
for t, c, e in sets: print('  %-10s %8.4f  %.2e' % (t, c, e))
print('  projected  %8.4f  %.2e' % pj); print('  Volkov MC  %8.1f  %.2e (mean sigma_rel)' % vk)
fig, ax = plt.subplots(figsize=(7.5, 5.2), dpi=150)
ax.grid(True, which='major', color='#e6e6e6', lw=0.6); ax.set_axisbelow(True)
Am = vk[0]*vk[1]**2     # MC scaling through Volkov's mean point: sigma = (Am/cost)^(1/2)
cs = [1e-2, 300]; ax.plot(cs, [(Am/c)**0.5 for c in cs], color='#888', lw=1.5, ls='--', label=r'Monte Carlo through Volkov mean, $\sigma \propto$ cost$^{-1/2}$')
ax.scatter([vk[0]], [vk[1]], s=60, color='#888', marker='s', label='Volkov 2019, mean sigma (V100)', zorder=3, edgecolor='white', lw=0.8)
ax.plot([c for _, c, _ in sets], [e for _, _, e in sets], color=O, lw=1.5, marker='o', ms=7, mec='white', label='sinc, fixed h at matched eps', zorder=3)
for t, c, e in sets: ax.annotate(t, (c, e), xytext=(6, 2), textcoords='offset points', fontsize=8, color='#555')
ax.scatter([pj[0]], [pj[1]], s=90, color=B, marker='D', edgecolor='white', lw=0.8, zorder=4, label='sinc, projected: LSQ over h 2.4-2.1 (+2.0 for 5 graphs) + eps 1e-5 run, all charged')
ax.annotate('projected', pj, xytext=(6, 2), textcoords='offset points', fontsize=8, color='#555')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e-2, 300); ax.set_ylim(1e-4, 0.5)
ax.set_xlabel('mean cost per diagram (TFLOP-hours)'); ax.set_ylabel('mean relative error per diagram')
ax.set_title('Five-loop g5 graphs (11): per-diagram mean error vs cost', fontsize=11, loc='left', color='#222')
ax.legend(fontsize=8, frameon=False, loc='upper right', borderaxespad=0.8)
fig.text(0.01, 0.005, 'mean of |S - Volkov|/|Volkov| over the 11 graphs (Volkov: mean of his sigma/|value|). Both at nominal FP64 peak: Mac Studio %.2f, V100 7.8 TFLOP/s.' % MACPEAK, fontsize=6.5, color='#777')
fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig(out + '_mean.png'); fig.savefig(out + '_mean.pdf')
