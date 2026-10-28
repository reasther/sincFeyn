"""Set the quadrature rule over the external-vertex split u in a generated kernel.

The rule integrates the split u in [0,1] of the electron line at the external
vertex (x1 = sigma*u, x2 = sigma*(1-u)).  At zero momentum transfer the point
function is linear in u, f = f0 + f1 u, times the simplex weight
(1-u)^p u^q (1, sa, sa*sa/2, sb*sa, ...).  So one point is exact: the weighted
mean u = (q+1)/(p+q+2), with weight B(p,q)/((1-u)^p u^q), B(p,q) =
p! q!/(p+q+1)!.  This is the default.  -n N sets an N-point Gauss-Legendre rule
instead (exact from N = 2 up to cubic weights; the generators emit 2 to 6).
Checked to roundoff on every kernel this script rewrites (tests/check_onept.sh).

  python3 tools/glpts.py kernels-ext/g4/*.c          # one point, in place
  python3 tools/glpts.py -n 2 kernels-ext/g4/*.c     # back to 2 Gauss-Legendre points

Rewrites the scalar table and loop bound, the SIMD lane tables (AVX 4 lanes,
SSE/NEON 2 lanes) and the header comment.  Idempotent; no other code changes.
A file is only rewritten when every GLX table in it has a matching evaluation
loop, so a table used by something this script does not understand is left
alone rather than silently desynchronised from its loop.
"""
import math, re, sys
from fractions import Fraction

def rule(n):                      # Gauss-Legendre on [0,1], Newton on P_n
    xs, ws = [], []
    for i in range(n):
        x = math.cos(math.pi*(i + 0.75)/(n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for k in range(2, n + 1): p0, p1 = p1, ((2*k - 1)*x*p1 - (k - 1)*p0)/k
            dp = n*(x*p1 - p0)/(x*x - 1) if n > 1 else 1.0
            dx = (p1 if n > 1 else x)/dp; x -= dx
            if abs(dx) < 1e-17: break
        xs.append((1 - x)/2); ws.append(1/((1 - x*x)*dp*dp))
    return [repr(v) for v in xs], [repr(v) for v in ws]

def mean_rule(src, t):
    """one point at the weighted mean of u for table GLX<t>, from the weight expression wg=(REAL)GLW<t>[g]..."""
    m = re.search(r'wg=\(REAL\)GLW%s\[gg?\]([^;,]*)[;,]' % re.escape(t), src)
    if not m: return None
    e = m.group(1)
    if re.sub(r'\bs[ab]\b|[*/]|\d+|\s', '', e): return None    # only products of sa, sb and integer divisors
    p = len(re.findall(r'\bsa\b', e)); q = len(re.findall(r'\bsb\b', e))
    u = Fraction(q + 1, p + q + 2)
    B = Fraction(math.factorial(p)*math.factorial(q), math.factorial(p + q + 1))
    return [repr(float(u))], [repr(float(B/((1 - u)**p * u**q)))]

def lanes(xs, ws, vw):
    ng = -(-len(xs)//vw)
    X, W = [], []
    for g in range(ng):
        cx, cw = xs[g*vw:(g+1)*vw], ws[g*vw:(g+1)*vw]
        cx += ['0.5']*(vw - len(cx)); cw += ['0']*(vw - len(cw))
        X.append('{' + ','.join(cx) + '}'); W.append('{' + ','.join(cw) + '}')
    return ng, ','.join(X), ','.join(W)

def transform(src, n):
    """Return (new source, note).  The source is None when nothing was done."""
    tabs = re.findall(r'^static const double GLX(\w*)\[\d+\]=', src, re.M)
    if not tabs: return None, 'no GL table'
    rules = {}
    for t in tabs:                                    # every table needs its loop
        if not re.search(r'for\(int (gg?)=0;\1<\d+;\1\+\+\)\{ const REAL u=\(REAL\)GLX%s\[\1\]' % re.escape(t), src):
            return None, 'GLX%s has no recognised loop' % t
        rules[t] = rule(n) if n else mean_rule(src, t)
        if rules[t] is None: return None, 'GLX%s: weight not a product of sa, sb' % t
    for t in tabs:
        xs, ws = rules[t]; n_t = len(xs)
        src = re.sub(r'^static const double GLX%s\[\d+\]=\{[^}]*\}, GLW%s\[\d+\]=\{[^}]*\};' % (re.escape(t), re.escape(t)),
                     lambda _: 'static const double GLX%s[%d]={%s}, GLW%s[%d]={%s};' % (t, n_t, ','.join(xs), t, n_t, ','.join(ws)), src, count=1, flags=re.M)
        src = re.sub(r'for\(int (gg?)=0;\1<\d+;\1\+\+\)\{ const REAL u=\(REAL\)GLX%s\[\1\]' % re.escape(t),
                     lambda mm: re.sub(r'<\d+;', '<%d;' % n_t, mm.group(0)), src)
        def vec(mm):                                  # vecg.py form: NGV with the lane count
            vw = 4 if mm.group(1) == '4' else 2
            ng, X, W = lanes(list(xs), list(ws), vw)
            return '%s\n#define NGV%s %d\nstatic const V4 GLXV%s[%d]={%s}, GLWV%s[%d]={%s};' % (mm.group(0).split('\n')[0], t, ng, t, ng, X, t, ng, W)
        src = re.sub(r'^#(?:if|else)(?: SQ_VW == (4))?\n#define NGV%s \d+\nstatic const V4 GLXV%s\[\d+\]=.*;$' % (re.escape(t), re.escape(t)),
                     vec, src, flags=re.M)
        def vec4(mm):                                 # vec4.py form: lane count as n/SQ_VW in the loop
            vw = 4 if mm.group(1) == '4' else 2
            ng, X, W = lanes(list(xs), list(ws), vw)
            return '%s\nstatic const V4 GLXV%s[%d]={%s}, GLWV%s[%d]={%s};' % (mm.group(0).split('\n')[0], t, ng, X, t, ng, W)
        src = re.sub(r'^#(?:if|else)(?: SQ_VW == (4))?\nstatic const V4 GLXV%s\[\d+\]=.*;$' % re.escape(t),
                     vec4, src, flags=re.M)
        src = re.sub(r'for\(int (gg?)=0;\1<[\d()+* /A-Z_-]*;\1\+\+\)\{ const V4 u=GLXV%s\[\1\]' % re.escape(t),
                     lambda mm: re.sub(r'<[^;]*;', '<(%d+SQ_VW-1)/SQ_VW;' % n_t, mm.group(0)), src)
    k = n or 1
    src = re.sub(r'GL points \d+ \*/', 'GL points %d */' % k, src, count=1)
    return src, ('%d GL points' % n) if n else 'one point at the weighted mean'


if __name__ == '__main__':
    args = sys.argv[1:]; n = 0          # 0: one point at the weighted mean
    if args[:1] == ['-n']: n = int(args[1]); args = args[2:]
    for fn in args:
        s = open(fn).read(); t, note = transform(s, n)
        if t is None: print('%s: %s, skipped' % (fn, note)); continue
        if t != s: open(fn, 'w').write(t)
        print('%s: %s' % (fn, note))
