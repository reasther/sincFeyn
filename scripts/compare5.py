"""Compare finished five-loop results with Volkov's 2019 values, after the runs (no VOLKOV_SUPP needed while running).

    python3 scripts/compare5.py SUPPDIR [RESULTS ...]

RESULTS are sincqed results files or folders of them (default runs5/); every data line in them counts, the
last one wins for a repeated (graph, h, eps).  SUPPDIR holds Volkov's class_*_aver.txt (supplement of
PRD 100, 096004).  Prints one line per graph and setting: value, prune estimate, evaluations, time, Volkov's
value with its errors, the difference (absolute, %, and in units of the prune estimate), and then for every
(h, eps) at which all five graphs of the class abcdeeabcd (g5_3086 g5_1784 g5_372 g5_2490 g5_1078) are done,
the class sum against Volkov's -0.240488.
"""
import sys, os, re, glob, subprocess

FLOAT = re.compile(r'^[-+][0-9]\.[0-9]+e[-+][0-9]+$')
CLASS = ['g5_3086', 'g5_1784', 'g5_372', 'g5_2490', 'g5_1078']

def files(args):
    for a in args or ['runs5']:
        if os.path.isdir(a): yield from sorted(glob.glob(os.path.join(a, '*.txt')))
        elif os.path.exists(a): yield a

def read(paths):
    res = {}
    for p in paths:
        for line in open(p, errors='replace'):
            if line.startswith('#') or not line.startswith('g5_'): continue
            t = line.split()
            i = next((k for k in range(1, len(t)) if FLOAT.match(t[k]) and len(t[k]) > 18), None)
            if i is None or i + 8 >= len(t): continue
            try: h, eps = float(t[i - 4]), float(t[i - 2])
            except ValueError: continue
            res[(t[0], h, eps)] = dict(v=float(t[i]), prune=float(t[i + 3]), evals=t[i + 7], time=float(t[i + 8]))
    return res

def volkov(name, supp):
    here = os.path.dirname(os.path.abspath(__file__))
    out = subprocess.run([sys.executable, os.path.join(here, 'g5info.py'), name, supp],
                         capture_output=True, text=True).stdout.strip().split('\t')
    if len(out) < 8 or out[5] == '-': return None
    return float(out[5]), out[6], out[7], out[2]

def main():
    if len(sys.argv) < 2: print(__doc__); sys.exit(1)
    supp = sys.argv[1]; res = read(files(sys.argv[2:]))
    if not res: print('no five-loop results found'); sys.exit(1)
    ref = {}
    print('%-8s %3s %4s %7s  %15s %9s %12s %9s  %15s  %10s %8s %8s' % ('graph', 'D', 'h', 'eps', 'value', 'prune',
          'evals', 'time_s', 'Volkov', 'diff', '%', '/prune'))
    for (g, h, eps) in sorted(res, key=lambda k: (int(k[0][3:]), -k[1], -k[2])):
        r = res[(g, h, eps)]
        if g not in ref: ref[g] = volkov(g, supp)
        vk = ref[g]
        if vk is None:
            print('%-8s %3s %4.1f %7.0e  %+15.9f %9.1e %12s %9.0f' % (g, '?', h, eps, r['v'], r['prune'], r['evals'], r['time'])); continue
        d = r['v'] - vk[0]
        print('%-8s %3s %4.1f %7.0e  %+15.9f %9.1e %12s %9.0f  %+15.6f  %+10.6f %+7.2f%% %+8.2f' % (g, vk[3], h, eps, r['v'],
              r['prune'], r['evals'], r['time'], vk[0], d, 100 * d / vk[0], d / r['prune'] if r['prune'] else 0))
    print()
    for (h, eps) in sorted({(k[1], k[2]) for k in res}, key=lambda x: (-x[0], -x[1])):
        if all((g, h, eps) in res for g in CLASS):
            s = sum(res[(g, h, eps)]['v'] for g in CLASS); p = sum(res[(g, h, eps)]['prune'] for g in CLASS)
            print('class abcdeeabcd  h=%.1f eps=%.0e  sum %+.6f  prune %.1e  Volkov -0.240488  diff %+.6f' % (h, eps, s, p, s + 0.240488))

if __name__ == '__main__':
    main()
