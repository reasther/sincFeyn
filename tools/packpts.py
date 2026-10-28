#!/usr/bin/env python3
"""Add a packed-point wrapper fv_<name> to a vectorised kernel, in place:
  python3 tools/packpts.py kernels/g3/*.c kernels-ext/g4/*.c

The vector point function of the double half (pt4_ from vec4.py, pt_..._v from vecg.py) takes every
Sinc variable as a vector and computes lane by lane.  f_ fills the lanes with the Gauss-Legendre
nodes of one lattice point; with the one-point rule (glpts.py) that leaves one lane in use.  fv_
fills lane l with lattice point j0+l instead (u at the weighted mean in every lane), so all SQ_VW
lanes do work.  It is f_ for n points at once, lane for lane the same arithmetic:

  void fv_<name>(long n, const double *sf, long ld, double *out, double *ab)

point j has its Sinc variables at sf[j*ld ...]; out[j], ab[j] are f_'s value and |.| sum.  The graph
descriptor gets fv_<name> as its last field; the driver uses it in eval_batch when it is there.
Needs one GL point (run glpts.py first).  Idempotent; a file it does not recognise is left alone.
"""
import re, sys

def transform(src):
    """Return (new source, note).  The source is None when nothing was done."""
    md = re.search(r'^const graph_t G_(\w+)=\{"[^"]*","[^"]*",\d+,(\d+),[^}]*\};$', src, re.M)
    if not md: return None, 'no graph descriptor'
    name, D = md.group(1), int(md.group(2))
    if 'fv_' + name + '(' in src: return None, 'already packed'
    dh0 = src.find('#define REAL double'); dh1 = src.find('#define REAL QREAL')
    if dh0 < 0 or dh1 < 0: return None, 'no double half'
    dh = src[dh0:dh1]
    mf = re.search(r'static REAL FN\(f_%s\)\(const REAL \*sf, REAL \*ab\)\{\n(.*?)\n\}\n' % re.escape(name), dh, re.S)
    if not mf or 'GLXV_%s[g]' % name not in mf.group(1): return None, 'f_ is not vectorised'
    if not re.search(r'^static const double GLX_%s\[1\]=' % re.escape(name), src, re.M): return None, 'more than one GL point'
    b = mf.group(1)
    b = re.sub(r'\bsf\[', 'sv[', b)
    b = b.replace('(REAL)', '(double)')
    b = re.sub(r'\bREAL (sig|nrm)=', r'V4 \1=', b)
    b = re.sub(r'\bREAL jac=(\d+);', r'V4 jac=VB(\1);', b)
    b = re.sub(r'for\(int g=0;g<[^;]*;g\+\+\)', 'for(int g=0;g<1;g++)', b)
    b = b.replace('GLXV_%s[g]' % name, 'VB(GLX_%s[g])' % name).replace('GLWV_%s[g]' % name, 'VB(GLW_%s[g])' % name)
    tail = '*ab=RABS(nrm)*VSUM(aa); return nrm*VSUM(acc);'
    if b.count(tail) != 1: return None, 'unrecognised f_ ending'
    b = b.replace(tail, 'const V4 R=nrm*acc, A=VABS(nrm)*aa;')
    if re.search(r'\bREAL\b|RABS\(|VSUM\(|\bsf\b|\breturn\b', b): return None, 'unrecognised f_ body'
    body = '\n'.join('  ' + l for l in b.split('\n'))
    fv = ('/* packed points (tools/packpts.py): lane l holds lattice point j0+l */\n'
          'static void fv_%s(long n, const double *sf, long ld, double *out, double *ab){\n'
          '  for(long j0=0;j0<n;j0+=SQ_VW){\n'
          '    V4 sv[%d];\n'
          '    for(int i=0;i<%d;i++) for(int l=0;l<SQ_VW;l++) sv[i][l]=sf[(j0+l<n ? j0+l : n-1)*ld+i];\n'
          '%s\n'
          '    for(int l=0;l<SQ_VW && j0+l<n;l++){ out[j0+l]=R[l]; ab[j0+l]=A[l]; }\n'
          '  }\n'
          '}\n') % (name, D, D, body)
    dh = dh.replace(mf.group(0), mf.group(0) + fv)
    src = src[:dh0] + dh + src[dh1:]
    src = src.replace(md.group(0), md.group(0)[:-2] + ',fv_%s};' % name)
    return src, 'packed'

if __name__ == '__main__':
    for fn in sys.argv[1:]:
        s = open(fn).read(); t, note = transform(s)
        if t is not None and t != s: open(fn, 'w').write(t)
        print('%s: %s' % (fn, note))
