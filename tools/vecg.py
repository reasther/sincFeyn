#!/usr/bin/env python3
"""vecg.py KERNEL.c > OUT.c, or vecg.py -i KERNEL.c ... (in place): SIMD Gauss-Legendre nodes for kernels whose point function calls
helper functions (four-loop layout).  The double half gets a vector copy of every helper and of
pt_ (suffix _v); f_ runs the nodes as SQ_VW lanes, padding the last group with zero weights."""
import re, sys
hdr = '''#ifndef SQ_V4
#define SQ_V4
/* the Gauss-Legendre nodes of the double half run as SIMD lanes: 4 with AVX, else 2 (SSE, NEON) */
#if defined(__AVX__)
#define SQ_VW 4
#else
#define SQ_VW 2
#endif
typedef double V4 __attribute__((vector_size(8*SQ_VW)));
typedef long long V4L __attribute__((vector_size(8*SQ_VW)));
#define VABS(x) ((V4)((V4L)(x) & 0x7fffffffffffffffLL))
#define VZERO ((V4){0})
#define VONE (VZERO+1)
#define VB(x) ((x)+VZERO)
static inline double VSUM(V4 v){ double t = 0; for(int i = 0; i < SQ_VW; i++) t += v[i]; return t; }
#endif
'''
def transform(src):
    dh0 = src.index('#define REAL double'); dh1 = src.index('#define REAL QREAL')
    dh = src[dh0:dh1]
    mf = re.search(r'static REAL FN\(f_(\w+)\)\(const REAL \*sf, REAL \*ab\)\{\n(.*?)\n\}\n', dh, re.S)
    if not mf or 'pt4_' in dh or '#define FNV' in dh or mf.group(2).count('FN(pt_') != 1:
        return src
    name = mf.group(1); f = mf.group(2)
    ng = int(re.search(r'for\(int g=0;g<(\d+);g\+\+\)', f).group(1))
    # vector copy of everything between the FN define and f_
    a = dh.index('#define FN(x) x##_d\n') + len('#define FN(x) x##_d\n'); b = dh.index(mf.group(0))
    body = dh[a:b]
    v = re.sub(r'\bREAL\b', 'V4', body).replace('RABS(', 'VABS(').replace('FN(', 'FNV(')
    # a vector cannot be set from a bare scalar constant: wrap constant right-hand sides
    num = r'[-+]?\(?[-+0-9.eE*/() ]*\)?'
    def wrap(m):
        if re.fullmatch(r'\s*', m.group(2)): return m.group(0)
        if re.search(r'\bint\s+$', v[max(0, m.start() - 8):m.start()]): return m.group(0)   # loop counters
        return '%s=VB(%s)' % (m.group(1), m.group(2))
    v = re.sub(r'(?<![\w\]\[<>=!])(\b[A-Za-z_]\w*(?:\[\d+\])*)=(' + num + r')(?=[,;])', wrap, v)
    vec = '#define FNV(x) x##_v\n' + v + '#undef FNV\n'
    # f_: lanes over GL nodes
    ngr4 = (ng + 3) // 4; ngr2 = (ng + 1) // 2
    mg = re.search(r'GLX_%s\[%d\]=\{([^}]*)\}, GLW_%s\[%d\]=\{([^}]*)\}' % (name, ng, name, ng), src)
    x = mg.group(1).split(','); w = mg.group(2).split(',')
    def groups(k, n):
        xs = x + ['0.5'] * (k * n - ng); ws = w + ['0'] * (k * n - ng)
        return ','.join('{%s}' % ','.join(xs[i*n:(i+1)*n]) for i in range(k)), ','.join('{%s}' % ','.join(ws[i*n:(i+1)*n]) for i in range(k))
    gx4, gw4 = groups(ngr4, 4); gx2, gw2 = groups(ngr2, 2)
    consts = ('#if SQ_VW == 4\n#define NGV_%s %d\nstatic const V4 GLXV_%s[%d]={%s}, GLWV_%s[%d]={%s};\n#else\n' % (name, ngr4, name, ngr4, gx4, name, ngr4, gw4) +
              '#define NGV_%s %d\nstatic const V4 GLXV_%s[%d]={%s}, GLWV_%s[%d]={%s};\n#endif\n' % (name, ngr2, name, ngr2, gx2, name, ngr2, gw2))
    f2 = f.replace('REAL s[', 'V4 s[')
    f2 = re.sub(r'(s\[\d+\])=([^;]+);', r'\1=VB(\2);', f2)
    f2 = f2.replace('REAL acc=0, aa=0, a;', 'V4 acc=VZERO, aa=VZERO, a;')
    f2 = re.sub(r'for\(int g=0;g<%d;g\+\+\)\{ const REAL u=\(REAL\)GLX_%s\[g\], sb=sig\*u, sa=sig\*\(1-u\);' % (ng, name),
                'for(int g=0;g<NGV_%s;g++){ const V4 u=GLXV_%s[g], sb=sig*u, sa=sig*(1-u);' % (name, name), f2)
    f2 = f2.replace('const REAL wg=(REAL)GLW_%s[g]' % name, 'const V4 wg=GLWV_%s[g]' % name)
    f2 = f2.replace('FN(pt_%s)(s,0,&a)' % name, 'pt_%s_v(s,VZERO,&a)' % name).replace('RABS(wg)', '(wg)')
    f2 = f2.replace('*ab=RABS(nrm)*aa; return nrm*acc;', '*ab=RABS(nrm)*VSUM(aa); return nrm*VSUM(acc);')
    if 'GLX_' in f2 or 'GLW_' in f2 or 'FN(pt_' in f2: return src   # another f_ layout (lbl4): leave it
    newf = 'static REAL FN(f_%s)(const REAL *sf, REAL *ab){\n%s\n}\n' % (name, f2)
    # the scalar double helpers stay (pt_ is in the graph descriptor)
    dh = dh.replace(mf.group(0), consts + vec + newf)
    return src[:dh0] + hdr + dh + src[dh1:]

if sys.argv[1] == '-i':
    for fn in sys.argv[2:]:
        src = open(fn).read(); new = transform(src)
        if new != src: open(fn, 'w').write(new); print('vectorised', fn)
else:
    sys.stdout.write(transform(open(sys.argv[1]).read()))
