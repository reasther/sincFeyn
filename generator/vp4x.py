"""
Two-loop (fourth-order) vacuum-polarization insertion Pi4 on one photon of any vertex graph,
optionally together with one-loop insertions Pi2 on photons of the same graph.

Generalizes vp4.py (which does Pi4 in the one-loop vertex only).  The Pi4 loop (topologies 'V'
and 'S' of vp4.diagrams(), with their internal U subtractions) gives, per term k and Wick grade w,
(C0 + C1 k^2) e^{c k^2}.  On the photon j of the base graph (Schwinger parameter s_j):
    main family    : s_j -> s_j + c  minus unshifted (Pi(k^2) - Pi(0)), Gamma index +1-w
    contact family : D(k) k^2 e^{c k^2} = -e^{c k^2}: s_j is not integrated but set to c
                     (plus any Pi2 shifts on the same photon), Gamma index -w.
Cheng-Wu: the main family keeps the base wrapper (photon of arc 0 fixed to 1); the contact family
fixes the external-class sum sigma = 1 instead, so it works also when j is the arc-0 photon.
Pi2 insertions enter exactly as in vp.py (shift minus unshifted, weight omega, exponent a+b).

    python3 vp4x.py OUTDIR [check]
      writes w4_*.c and vp4xtable_4.h: Pi2 Pi4 in the one-loop vertex, Pi4 in the two-loop vertices.
      'check' instead writes w3_*.c: Pi4 in the one-loop vertex through this code, to compare with
      v3_p4* (vp4.py).
"""
import os, sys, pickle, math
from fractions import Fraction
import graphs as gr, build, vp4, vp

NW = build.NW

def yset(cls, src):
    o = []
    for t, c in enumerate(cls):
        for i in c: o.append("  y[%d]=%s[%d]/%d;" % (i, src, t, len(c)))
        if len(c) > 1: o.append("  jac*=%s[%d]/%d;" % (src, t, math.factorial(len(c)-1)))
    return o

def point_fn(g, cname, base, j, p2, key, d, tabs, contact, inner=False):
    """Point function FN(p_<cname>)(s, v, &ab): v = Pi4 class variables, then (a,b) per Pi2 on the
    base graph, then (a,b) of a Pi2 on the internal photon of the Pi4 loop if inner (the Pi4 loop
    is then evaluated with that photon shifted, minus unshifted)."""
    code, Ms = vp4.gen_terms_code(d, tabs)
    cls = d['classes']; nv4 = len(cls); np2 = len(p2); nt = len(d['terms'])
    h = ["static void FN(vp_%s)(const REAL *y, REAL *T, REAL *Ta, REAL *V, REAL *c){" % cname]
    h.append("  const REAL %s;" % ", ".join("y%d=y[%d]" % (i, i) for i in range(5)))
    h.append(code)
    for k in range(nt):
        h.append("  for(int m=0;m<2;m++) for(int w=0;w<%d;w++){ T[(%d*2+m)*%d+w]=m<%d?T%d[m][w]:0; Ta[(%d*2+m)*%d+w]=m<%d?T%da[m][w]:0; }" % (NW, k, NW, Ms[k], k, k, NW, Ms[k], k))
        h.append("  V[%d]=V%d; c[%d]=c%d;" % (k, k, k, k))
    h.append("}")
    npass = 2 if inner else 1
    o = ["static REAL FN(p_%s)(const REAL *x, const REAL *v, REAL *ab){" % cname]
    o.append("  REAL y[5], jac=1;")
    o += yset(cls, 'v')
    o.append("  REAL om=%d*jac, Vx=0, pc2_[%d];" % (d['weight'], max(np2, 1)))
    for t in range(np2):
        o.append("  { const REAL a=v[%d], b=v[%d], S=a+b; pc2_[t_]=a*b/S; om*=(REAL)(%d)*a*b/(S*S*S*S); Vx+=S; }".replace('t_', str(t)) % (nv4 + 2*t, nv4 + 2*t + 1, vp.OMEGA))
    o.append("  REAL T[%d][%d], Ta[%d][%d], V[%d][%d], c[%d][%d];" % (npass, nt*2*NW, npass, nt*2*NW, npass, nt, npass, nt))
    if inner:
        o.append("  { const REAL a=v[%d], b=v[%d], S=a+b; om*=(REAL)(%d)*a*b/(S*S*S*S); Vx+=S;" % (nv4 + 2*np2, nv4 + 2*np2 + 1, vp.OMEGA))
        o.append("    const REAL y4=y[4]; y[4]=y4+a*b/S; FN(vp_%s)(y,T[0],Ta[0],V[0],c[0]); y[4]=y4; FN(vp_%s)(y,T[1],Ta[1],V[1],c[1]); }" % (cname, cname))
    else:
        o.append("  FN(vp_%s)(y,T[0],Ta[0],V[0],c[0]);" % cname)
    o.append("  REAL s[%d], tot=0, tota=0, t, ta;" % g.n)
    pj = 2*g.L + j
    m = 1 if contact else 0
    for P in range(npass):
        psg = 1 if P == 0 else -1
        for k in range(nt):
            if contact and Ms[k] < 2: continue
            o.append("  for(int w=0;w<%d;w++){ const REAL Tw=T[%d][(%d*2+%d)*%d+w], Taw=Ta[%d][(%d*2+%d)*%d+w]; if(Tw==0 && Taw==0) continue;" % (NW, P, k, m, NW, P, k, m, NW))
            for mask in range(2**np2):
                sgn = (-1)**(np2 - bin(mask).count('1'))
                shifts = [t for t in range(np2) if (mask >> t) & 1]
                fams = [(-1, None)] if contact else [(1, True), (-1, False)]
                for fs, shifted in fams:
                    o.append("    for(int i=0;i<%d;i++) s[i]=x[i];" % g.n)
                    if contact: o.append("    s[%d]=c[%d][%d];" % (pj, P, k))
                    elif shifted: o.append("    s[%d]+=c[%d][%d];" % (pj, P, k))
                    for t in shifts: o.append("    s[%d]+=pc2_[%d];" % (2*g.L + p2[t], t))
                    dp = "-w" if contact else "1-w"
                    o.append("    t=FN(pt_%s)(s,V[%d][%d]+Vx,%s,&ta); tot%s=Tw*t; tota+=(Taw+RABS(Tw))*ta;" % (base, P, k, dp, '+' if sgn*fs*psg > 0 else '-'))
            o.append("  }")
    o.append("  *ab=RABS(om)*tota; return om*tot;\n}")
    return "\n".join(h) + "\n" + "\n".join(o), nv4

def wrapper_contact(g, cname, inner, nextra, j):
    """Like build.gen_wrapper, but Cheng-Wu sigma = 1 and photon j not integrated (set by inner)."""
    (bef, aft), others, fixed = build.classes(g)
    others = [c for c in others if c != [2*g.L + j]]
    if fixed != 2*g.L + j: others.append([fixed])
    D = len(others) + nextra; k1, k2 = len(bef), len(aft)
    deg = g.ne + k1 + k2 - 2; ngl = (deg + 2)//2
    xs, ws = build.gauss_legendre01(ngl)
    norm = Fraction(-1, 4)**g.L / 16
    glc = "static const double GLX_%s[%d]={%s}, GLW_%s[%d]={%s};" % (cname, ngl, ",".join(repr(float(x)) for x in xs), cname, ngl, ",".join(repr(float(w)) for w in ws))
    o = ["static REAL FN(f_%s)(const REAL *sf, REAL *ab){" % cname]
    o.append("  REAL s[%d]; s[%d]=0; REAL jac=1;" % (g.n, 2*g.L + j))
    o.append("  const REAL *v=sf+%d;" % (D - nextra))
    for d, c in enumerate(others):
        kk = len(c)
        for i in c: o.append("  s[%d]=sf[%d]/%d;" % (i, d, kk))
        if kk > 1: o.append("  jac*=%s;" % "*".join(["sf[%d]" % d]*(kk-1) + ["((REAL)1/%d)" % math.factorial(kk-1)]))
    o.append("  REAL acc=0, aa=0, a;")
    o.append("  for(int g=0;g<%d;g++){ const REAL u=(REAL)GLX_%s[g], sb=u, sa=1-u;" % (ngl, cname))
    for i in bef: o.append("    s[%d]=sb/%d;" % (i, k1))
    for i in aft: o.append("    s[%d]=sa/%d;" % (i, k2))
    wexpr = "(REAL)GLW_%s[g]" % cname
    if k1 > 1: wexpr += "*" + "*".join(["sb"]*(k1-1)) + "/%d" % math.factorial(k1-1)
    if k2 > 1: wexpr += "*" + "*".join(["sa"]*(k2-1)) + "/%d" % math.factorial(k2-1)
    o.append("    const REAL wg=%s;" % wexpr)
    o.append("    acc+=wg*FN(%s)(s,v,&a); aa+=RABS(wg)*a; }" % inner)
    o.append("  const REAL nrm=%s*jac;" % build.cnum(norm))
    o.append("  *ab=RABS(nrm)*aa; return nrm*acc;\n}")
    return "\n".join(o), D, ngl, glc

def emit(g, w0, prefix, j, p2, cache, out, title, Ltot, inner=False):
    """Files <prefix><key><m|c>.c for the Pi4 topologies key = V, S.  Returns the names."""
    base = prefix + '_b'
    pt, info = build.gen_graph(g, base, cache, dpvar=True)
    names = []
    for key, d in vp4.diagrams().items():
        tabs = {n.name: vp4.table(n, cache) for n in d['nodes']}
        # contact family only if some Pi node has a k^2 (C1) coefficient; the V topology has none
        hasc = any(len(v) > 1 and any(x != 0 for x in v[1:]) for n in d['nodes'] if n.kind == 'PI' for v in tabs[n.name].values())
        for contact in ((False, True) if hasc else (False,)):
            cn = '%s%s%s' % (prefix, key, 'c' if contact else 'm')
            pf, nv4 = point_fn(g, cn, base, j, p2, key, d, tabs, contact, inner)
            nextra = nv4 + 2*len(p2) + (2 if inner else 0)
            if contact: wr, D, ngl, glc = wrapper_contact(g, cn, 'p_%s' % cn, nextra, j)
            else: wr, D, ngl, glc = build.gen_wrapper(g, cn, inner='p_%s' % cn, nextra=nextra)
            tt = '%s Pi4%s%s(%d-%d) %s' % (title, key, '[+Pi2 inside]' if inner else '', g.arcs[j][0]+1, g.arcs[j][1]+1, 'contact' if contact else 'main')
            src = ["/* Generated by vp4x.py: %s  (base graph %s, weight %d, D=%d) */" % (tt, g.volkov(), w0, D),
                   '#include "common.h"', glc,
                   '#define REAL double', '#define FN(x) x##_d', pt, pf, wr, '#undef REAL', '#undef FN',
                   '#define REAL QREAL', '#define FN(x) x##_q', pt, pf, wr, '#undef REAL', '#undef FN',
                   'const graph_t G_%s={"%s","%s",%d,%d,%d,%d,%d,f_%s_d,f_%s_q,0,0};' %
                   (cn, cn, tt, Ltot, D, w0, g.n, 2*g.L, cn, cn)]
            open(os.path.join(out, cn + '.c'), 'w').write("\n".join(src) + "\n")
            names.append(cn); print(cn, w0, D, tt, flush=True)
    return names

def main(out, check=False):
    cf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache_dirac.pkl')
    cache = pickle.load(open(cf, 'rb'))
    (g1, w1), = gr.enumerate_graphs(1)
    names = []
    if check:
        names += emit(g1, 1, 'w3_', 0, [], cache, out, 'VP: 1-loop vertex', 3)
        L = 3
    else:
        L = 4
        names += emit(g1, 2, 'w4_00', 0, [0], cache, out, 'VP: 1-loop vertex +Pi2(1-2)', 4)   # 2 orderings of Pi2, Pi4 on the photon
        names += emit(g1, 1, 'w4_0i', 0, [], cache, out, 'VP: 1-loop vertex', 4, inner=True)   # Pi2 inside the Pi4
        i = 1
        for g, w in gr.enumerate_graphs(2):
            for ms, wt in vp.photon_sets(g, w, 1):
                jj = ms[0]
                names += emit(g, wt, 'w4_%02d' % i, jj, [], cache, out, 'VP: %s' % g.volkov(), 4); i += 1
    tb = ['/* Generated by vp4x.py */', '#include "common.h"'] + ['extern const graph_t G_%s;' % c for c in names]
    tb.append('static const graph_t *VP4XGRAPHS_%d[]={%s};' % (L, ",".join('&G_%s' % c for c in names)))
    tb.append('#define NVP4XGRAPHS_%d %d' % (L, len(names)))
    open(os.path.join(out, 'vp4xtable_%d.h' % L), 'w').write("\n".join(tb) + "\n")
    pickle.dump(cache, open(cf + '.tmp', 'wb')); os.replace(cf + '.tmp', cf)
    print('Pi4 insertions:', names)

if __name__ == '__main__':
    main(sys.argv[1], len(sys.argv) > 2 and sys.argv[2] == 'check')
