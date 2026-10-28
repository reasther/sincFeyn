#!/usr/bin/env python3
"""vec4.py KERNEL.c > OUT.c : evaluate the 4 Gauss-Legendre nodes of the double half as one
4-wide vector (GCC vector extension).  Kernels with another node count are left unchanged."""
import re, sys
src = open(sys.argv[1]).read()
def unchanged():
    sys.stdout.write(src); sys.exit(0)
dh0 = src.index('#define REAL double'); dh1 = src.index('#define REAL QREAL')
dh = src[dh0:dh1]
mf = re.search(r'static REAL FN\(f_(\w+)\)\(const REAL \*sf, REAL \*ab\)\{\n(.*?)\n\}\n', dh, re.S)
if not mf or 'for(int g=0;g<4;g++)' not in mf.group(2):
    sys.stdout.write(src); sys.exit(0)
name = mf.group(1)
if mf.group(2).count('FN(pt_') != 1: unchanged()
mp = re.search(r'static REAL FN\(pt_%s\)\(const REAL \*x, const REAL Vx, REAL \*ab\)\{\n(.*?)\n\}\n' % re.escape(name), dh, re.S)
body = mp.group(1)
assert 'FN(' not in body
body = re.sub(r'\bREAL\b', 'V4', body).replace('RABS(', 'VABS(')
body = re.sub(r'(?<!int )(?<![\w\]])(\b[A-Za-z_]\w*)=0(?=[,;])', r'\1=VZERO', body)
body = re.sub(r'(?<!int )(?<![\w\]])(\b[A-Za-z_]\w*)=1(?=[,;])', r'\1=VONE', body)
body = re.sub(r'(?<!int )(?<![\w\]])(\b[A-Za-z_]\w*)=(-?\d+)(?=[,;])', r'\1=VB(\2)', body)   # FORM: E0=4;
pt4 = 'static V4 pt4_%s(const V4 *x, const double Vx, V4 *ab){\n%s\n}\n' % (name, body)
f = mf.group(2)
f = f.replace('REAL s[9];', 'V4 s[9];')
f = re.sub(r'(s\[\d+\])=([^;]+);', r'\1=VB(\2);', f)
f = f.replace('REAL acc=0, aa=0, a;', 'V4 acc=VZERO, aa=VZERO, a;')
f = re.sub(r'for\(int g=0;g<4;g\+\+\)\{ const REAL u=\(REAL\)GLX_%s\[g\]' % re.escape(name), 'GLXV_%s' % name, f)
f = f.replace('GLXV_%s' % name, 'for(int g=0;g<4/SQ_VW;g++){ const V4 u=GLXV_%s[g]' % name, 1)
f = f.replace('const REAL wg=(REAL)GLW_%s[g]' % name, 'const V4 wg=GLWV_%s[g]' % name)
f = f.replace('FN(pt_%s)(s,0,&a)' % name, 'pt4_%s(s,0,&a)' % name).replace('RABS(wg)', '(wg)')
f = f.replace('*ab=RABS(nrm)*aa; return nrm*acc;', '*ab=RABS(nrm)*VSUM(aa); return nrm*VSUM(acc);')
assert 'GLX_' not in f and 'GLW_' not in f and 'FN(pt_' not in f, f
mg = re.search(r'GLX_%s\[4\]=\{([^}]*)\}, GLW_%s\[4\]=\{([^}]*)\}' % (name, name), src)
x = mg.group(1).split(','); w = mg.group(2).split(',')
consts = ('#if SQ_VW == 4\nstatic const V4 GLXV_%s[1]={{%s}}, GLWV_%s[1]={{%s}};\n#else\n' % (name, ','.join(x), name, ','.join(w)) +
          'static const V4 GLXV_%s[2]={{%s},{%s}}, GLWV_%s[2]={{%s},{%s}};\n#endif\n' % (name, ','.join(x[:2]), ','.join(x[2:]), name, ','.join(w[:2]), ','.join(w[2:])))
newf = 'static REAL FN(f_%s)(const REAL *sf, REAL *ab){\n%s\n}\n' % (name, f)
dh = dh.replace(mf.group(0), consts + pt4 + newf)
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
out = src[:dh0] + hdr + dh + src[dh1:]
sys.stdout.write(out)
