"""
One-loop vertex with a two-loop (fourth-order) vacuum-polarization insertion: 3 graphs
(vertex-type VP 'V', self-energy-type VP 'S' counted twice by symmetry).

The VP loop is renormalized inside by Volkov's U operators on its divergent subgraphs (the two
overlapping one-loop vertices of V, the one-loop self-energy of S; U satisfies a = -v, the Ward
identity, so the internal counterterms are consistent) and overall by Pi(k^2) - Pi(0).
Pi is read off the gauge-invariant k^mu k^nu coefficient of the loop, which for two loops is a
polynomial C_0 + C_1 k^2 times e^{c k^2}.  In the parent photon line D(k) = int ds e^{s k^2}:
    C_0 part :  s -> s + c   (minus unshifted: the Pi(0) subtraction)          'main' family
    C_1 part :  D(k) k^2 e^{c k^2} = -e^{c k^2}: photon parameter fixed to c    'contact' family
The two families are separate integrals (the contact one has one dimension less); both are
registered and must be summed.  Weight of the VP: -(-1/4)^2 C / U^2 (fermion loop sign and loop
factors; the one-loop analogue reproduces 119/36 - pi^2/3).

    python3 vp4.py OUTDIR       -> v3_p4*.c, vp4table.h, and pi4slope.c (standalone check of the
                                   renormalized Pi^(4) slope, expected -41/162 in these units,
                                   where the one-loop slope is -1/15)
"""
import os, sys, pickle, math, hashlib
from fractions import Fraction
import sympy as sp
import dirac, build, graphs as gr
from symb import symbolic

class N:            # minimal node object for build.emit_node / emit_pi
    pass

def mk(name, kind, A, B, m2, gidx, chain, nslots, slotline=None, selist=(), loop=None):
    n = N(); n.name, n.kind, n.A, n.B, n.m2, n.gidx = name, kind, A, B, m2, gidx
    n.chain = chain; n.slots = [None]*nslots; n.slotline = slotline or list(range(nslots))
    n.selist = list(selist); n.egaps = list(range(nslots - len(selist))); n.L = len(A[0]); n.seg = name; n.kids = ()
    n.loop = loop
    return n

def table(n, cache):
    key = repr((n.chain, len(n.slots), n.kind, n.loop))
    h = hashlib.md5(key.encode()).hexdigest()
    if h not in cache:
        if n.kind == 'PI': cache[h] = dirac.table_pi(n.loop, len(n.slots))
        else: cache[h] = dirac.table(n.chain, len(n.slots), n.kind)
    return cache[h]

def emit_pi(n, tab):
    """Arrays <nm>_W[m][pat][w] (k^2 power m), <nm>_U, <nm>_V (= sum of electron parameters),
    <nm>_c (k^2 coefficient of the exponent)."""
    s, U, F, coefU, CU = symbolic(n.A, n.B, n.m2)
    sub = {s[i]: sp.Symbol('y%d' % n.gidx[i]) for i in range(len(s))}
    cc = lambda e: build.cexpr(sp.expand(sp.sympify(e).subs(sub)))
    nm = n.name; npat = 2**len(n.selist); M = max(len(v) for v in tab.values())
    V0 = sum(s[i]*n.m2[i] for i in range(len(s)))
    o = ["  /* VP node %s */" % nm]
    o.append("  REAL %s_W[%d][%d][%d]={{{0}}}, %s_Wa[%d][%d][%d]={{{0}}}, %s_U, %s_V, %s_c;" % (nm, M, npat, build.NW, nm, M, npat, build.NW, nm, nm, nm))
    o.append("  { const REAL U=%s, iU=1/U;" % cc(U))
    o.append("    const REAL V0=%s, c=((%s)*V0*0+(%s))*iU;" % (cc(V0), '0', cc(sp.expand(U*V0 - F))))
    nl = len(n.egaps); used = set(); pairs = set()
    # line identifications (build._monomials): ck by equal (A, B) rows, w by equal A rows
    ra, rab = build._reps(n) if build._ident() else (None, None)
    rk = (lambda i: rab[i]) if ra else (lambda i: i)
    rw = (lambda i, j: tuple(sorted((ra[i], ra[j])))) if ra else (lambda i, j: tuple(sorted((i, j))))
    for (combo, pr) in tab:
        for j, ch in enumerate(combo):
            if ch == 'k': used.add(rk(n.slotline[j]))
        for a, b in pr: pairs.add(rw(n.slotline[a], n.slotline[b]))
    for i in sorted(used): o.append("    const REAL ck%d=(%s)*iU;" % (i, cc(coefU[i][0])))
    for (i, j) in sorted(pairs): o.append("    const REAL w%d_%d=(%s)*iU;" % (i, j, cc(CU[i][j])))
    o.append("    REAL q;")
    maxw = 0; merged = {}
    for (combo, pr), vals in tab.items():
        w = len(pr); maxw = max(maxw, w)
        fac = ['ck%d' % rk(n.slotline[j]) for j, ch in enumerate(combo) if ch == 'k']
        fac += ['w%d_%d' % rw(n.slotline[a], n.slotline[b]) for a, b in pr]
        pat = 0
        for t in range(len(n.selist)):
            if combo[nl + t] != '1': pat |= 1 << t
        pw = Fraction(-1, 2)**w
        if ra:
            v = merged.setdefault((tuple(sorted(fac)), pat, w), [Fraction(0)]*len(vals))
            for m, x in enumerate(vals): v[m] += Fraction(x)*pw
            continue
        o.append("    q=%s;" % ('*'.join(fac) if fac else '1'))
        for m, x in enumerate(vals):
            if x != 0:
                o.append("    %s_W[%d][%d][%d]+=%s*q; %s_Wa[%d][%d][%d]+=RABS(%s*q);" % (nm, m, pat, w, build.cnum(x*pw), nm, m, pat, w, build.cnum(x*pw)))
    for (fac, pat, w), vals in merged.items():
        if not any(vals): continue
        o.append("    q=%s;" % ('*'.join(fac) if fac else '1'))
        for m, x in enumerate(vals):
            if x != 0:
                o.append("    %s_W[%d][%d][%d]+=%s*q; %s_Wa[%d][%d][%d]+=RABS(%s*q);" % (nm, m, pat, w, build.cnum(x), nm, m, pat, w, build.cnum(x)))
    o.append("    const REAL iU2=iU*iU;")
    o.append("    for(int m=0;m<%d;m++) for(int p=0;p<%d;p++) for(int w=0;w<%d;w++){ %s_W[m][p][w]*=iU2; %s_Wa[m][p][w]*=iU2; }" % (M, npat, build.NW, nm, nm))
    o.append("    %s_U=U; %s_V=V0; %s_c=c; }" % (nm, nm, nm))
    n.maxw = maxw; n.M = M
    return "\n".join(o)

def emit_inner(n, tab):
    """V0 / SE node with parameters named y<i>: reuse build.emit_node on renamed symbols."""
    code = build.emit_node(n, tab)
    import re
    return re.sub(r'\bx(\d+)\b', r'y\1', code)

# ---------------- the two VP topologies ----------------
# params y: V: 0 l1(mu->x) 1 l2(x->nu) 2 l3(nu->y) 3 l4(y->mu) 4 photon ; loop momenta (l, q), external k
def diagrams():
    D = {}
    # --- vertex type
    A = [[1, 0], [1, 1], [1, 1], [1, 0], [0, 1]]; B = [[1, 0], [1, 0], [0, 0], [0, 0], [0, 0]]; m2 = [1, 1, 1, 1, 0]
    loop = [('s', 3), ('v', 'g'), ('s', 2), ('v', 'nu'), ('s', 1), ('v', 'g'), ('s', 0), ('v', 'mu')]
    full = mk('pV', 'PI', A, B, m2, [0, 1, 2, 3, 4], None, 4, loop=loop)
    # U on vertex at mu: lines l4 (slot0), l1 (slot1), photon ; V0 chain
    uA = [[-1], [-1], [1]]; uB = [[1, 0], [1, 0], [0, 0]]; um2 = [1, 1, 0]
    ch = [('v', 'g'), ('s', 1), ('v', 'mu'), ('s', 0), ('v', 'g')]
    uM = mk('uVm', 'V0', uA, uB, um2, [3, 0, 4], ch, 2)
    rM = mk('rVm', 'PI', [[1], [1]], [[1, 0], [0, 0]], [1, 1], [1, 2], None, 2, loop=[('s', 1), ('v', 'nu'), ('s', 0), ('v', 'mu')])
    uN = mk('uVn', 'V0', uA, uB, um2, [1, 2, 4], ch, 2)
    rN = mk('rVn', 'PI', [[1], [1]], [[1, 0], [0, 0]], [1, 1], [0, 3], None, 2, loop=[('s', 1), ('v', 'nu'), ('s', 0), ('v', 'mu')])
    D['V'] = dict(nodes=[full, uM, rM, uN, rN], terms=[(1, full, None), (-1, rM, uM), (-1, rN, uN)],
                  classes=[[0], [1], [2], [3], [4]], weight=1)
    # --- self-energy type: 0 l1(mu->a) 1 l2(a->b) 2 l3(b->nu) 3 l4(nu->mu) 4 photon
    A = [[1, 0], [1, -1], [1, 0], [1, 0], [0, 1]]; B = [[1, 0], [1, 0], [1, 0], [0, 0], [0, 0]]
    loop = [('s', 3), ('v', 'nu'), ('s', 2), ('v', 'g'), ('s', 1), ('v', 'g'), ('s', 0), ('v', 'mu')]
    full = mk('pS', 'PI', A, B, m2, [0, 1, 2, 3, 4], None, 4, loop=loop)
    se = mk('uSs', 'SE', [[-1], [1]], [[1, 0], [0, 0]], [1, 0], [1, 4], [('v', 'g'), ('s', 0), ('v', 'g')], 1)
    # reduced: slots 0 l1, 1 l3, 2 l4, 3 insertion (momentum of l1)
    red = mk('rS', 'PI', [[1], [1], [1]], [[1, 0], [1, 0], [0, 0]], [1, 1, 1], [0, 2, 3], None, 4,
             slotline=[0, 1, 2, 0], selist=['se'], loop=[('s', 2), ('v', 'nu'), ('s', 1), ('s', 3), ('s', 0), ('v', 'mu')])
    D['S'] = dict(nodes=[full, se, red], terms=[(1, full, None), (-1, red, se)],
                  classes=[[0, 2], [1], [3], [4]], weight=2)
    return D

def gen_terms_code(d, tabs):
    """C code computing, for every VP term, T<k>[m][w] (with -(-1/4)^2 and the term sign), V<k>, c<k>."""
    o = []
    for n in d['nodes']:
        o.append(emit_pi(n, tabs[n.name]) if n.kind == 'PI' else emit_inner(n, tabs[n.name]))
    K = Fraction(-1, 16)
    NW = build.NW
    out = []
    for k, (sg, pin, inner) in enumerate(d['terms']):
        M = pin.M
        o.append("  REAL T%d[%d][%d]={{0}}, T%da[%d][%d]={{0}}, V%d, c%d;" % (k, M, NW, k, M, NW, k, k))
        o.append("  { const REAL f=%s;" % build.cnum(K*sg))
        if inner is None:
            o.append("    for(int m=0;m<%d;m++) for(int w=0;w<%d;w++){ T%d[m][w]=f*%s_W[m][0][w]; T%da[m][w]=RABS(f)*%s_Wa[m][0][w]; }" % (M, NW, k, pin.name, k, pin.name))
            o.append("    V%d=%s_V; c%d=%s_c; }" % (k, pin.name, k, pin.name))
        elif inner.kind == 'V0':      # U (a) of the vertex times the reduced VP
            o.append("    for(int m=0;m<%d;m++) for(int i=0;i<%d;i++) for(int j=0;i+j<%d;j++){" % (M, NW, NW))
            o.append("      T%d[m][i+j]+=f*%s_W[m][0][i]*%s_W[0][0][j]; T%da[m][i+j]+=RABS(f)*%s_Wa[m][0][i]*%s_Wa[0][0][j]; }" % (k, pin.name, inner.name, k, pin.name, inner.name))
            o.append("    V%d=%s_V+%s_V; c%d=%s_c; }" % (k, pin.name, inner.name, k, pin.name))
        else:                         # SE insertion: pattern 0 -> u, 1 -> v
            o.append("    for(int m=0;m<%d;m++) for(int p=0;p<2;p++) for(int i=0;i<%d;i++) for(int j=0;i+j<%d;j++){" % (M, NW, NW))
            o.append("      T%d[m][i+j]+=f*%s_W[m][p][i]*%s_W[p][0][j]; T%da[m][i+j]+=RABS(f)*%s_Wa[m][p][i]*%s_Wa[p][0][j]; }" % (k, pin.name, inner.name, k, pin.name, inner.name))
            o.append("    V%d=%s_V+%s_V; c%d=%s_c; }" % (k, pin.name, inner.name, k, pin.name))
        out.append(M)
    return "\n".join(o), out

def main(out):
    cf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache_dirac.pkl')
    cache = pickle.load(open(cf, 'rb'))
    (g1, w1), = gr.enumerate_graphs(1)
    pt, info = build.gen_graph(g1, 'p4par', cache, dpvar=True)
    xs, ws = build.gauss_legendre01(3)
    glc = "static const double GLX_p4[3]={%s}, GLW_p4[3]={%s};" % (",".join(repr(float(x)) for x in xs), ",".join(repr(float(w)) for w in ws))
    norm_par = Fraction(-1, 4)/16
    names = []; slope = []
    for key, d in diagrams().items():
        tabs = {n.name: table(n, cache) for n in d['nodes']}
        code, Ms = gen_terms_code(d, tabs)
        nt = len(d['terms']); ny = 5
        cls = d['classes']; nv = len(cls)
        def yset(off):
            o = []
            for t, c in enumerate(cls):
                for i in c: o.append("  y[%d]=sf[%d]/%d;" % (i, off + t, len(c)))
                if len(c) > 1: o.append("  jac*=sf[%d]/%d;" % (off + t, math.factorial(len(c)-1)))
            return o
        body = ["  const REAL %s;" % ", ".join("y%d=y[%d]" % (i, i) for i in range(ny))] + [code]
        # main family: parent photon = 1 (+c), parent ext class sigma = sf[0] (GL split)
        cm = 'v3_p4%sm' % key
        f = ["static REAL FN(f_%s)(const REAL *sf, REAL *ab){" % cm, "  REAL y[5], s[3], jac=1; const REAL sig=sf[0];"] + yset(1) + body
        f.append("  REAL acc=0, aa=0, a, v;")
        f.append("  for(int g=0;g<3;g++){ const REAL u=(REAL)GLX_p4[g], wg=(REAL)GLW_p4[g]; s[0]=sig*u; s[1]=sig*(1-u);")
        for k in range(nt):
            f.append("    for(int w=0;w<%d;w++){ if(T%d[0][w]==0 && T%da[0][w]==0) continue;" % (build.NW, k, k))
            f.append("      s[2]=1+c%d; v=FN(pt_p4par)(s,V%d,1-w,&a); acc+=wg*T%d[0][w]*v; aa+=wg*(T%da[0][w]*a+RABS(T%d[0][w])*a);" % (k, k, k, k, k))
            f.append("      s[2]=1;     v=FN(pt_p4par)(s,V%d,1-w,&a); acc-=wg*T%d[0][w]*v; aa+=wg*(T%da[0][w]*a+RABS(T%d[0][w])*a); }" % (k, k, k, k))
        f.append("  }")
        f.append("  const REAL nrm=%s*%d*sig*jac;" % (build.cnum(norm_par), d['weight']))
        f.append("  *ab=RABS(nrm)*aa; return nrm*acc;\n}")
        # contact family: parent photon parameter = c, parent sigma = 1 (Cheng-Wu), VP classes free
        cc_ = 'v3_p4%sc' % key
        g = ["static REAL FN(f_%s)(const REAL *sf, REAL *ab){" % cc_, "  REAL y[5], s[3], jac=1;"] + yset(0) + body
        g.append("  REAL acc=0, aa=0, a, v;")
        g.append("  for(int gg=0;gg<3;gg++){ const REAL u=(REAL)GLX_p4[gg], wg=(REAL)GLW_p4[gg]; s[0]=u; s[1]=1-u;")
        for k in range(nt):
            if Ms[k] < 2: continue
            f_ = "    if(T%d[1][0]!=0 || T%da[1][0]!=0){ s[2]=c%d; v=FN(pt_p4par)(s,V%d,0,&a); acc-=wg*T%d[1][0]*v; aa+=wg*(T%da[1][0]*a+RABS(T%d[1][0])*a); }" % (k, k, k, k, k, k, k)
            g.append(f_)
        g.append("  }")
        g.append("  const REAL nrm=%s*%d*jac;" % (build.cnum(norm_par), d['weight']))
        g.append("  *ab=RABS(nrm)*aa; return nrm*acc;\n}")
        fm = "\n".join(f); fc = "\n".join(g)
        for nm_, src_, D in ((cm, fm, 1 + nv), (cc_, fc, nv)):
            title = 'VP: 1-loop vertex + Pi4(%s) %s' % (key, 'main' if nm_.endswith('m') else 'contact')
            txt = ["/* Generated by vp4.py: %s (weight %d included in the integrand) */" % (title, d['weight']),
                   '#include "common.h"', glc,
                   '#define REAL double', '#define FN(x) x##_d', pt, src_, '#undef REAL', '#undef FN',
                   '#define REAL QREAL', '#define FN(x) x##_q', pt, src_, '#undef REAL', '#undef FN',
                   'const graph_t G_%s={"%s","%s",3,%d,1,8,2,f_%s_d,f_%s_q,0,0};' % (nm_, nm_, title, D, nm_, nm_)]
            open(os.path.join(out, nm_ + '.c'), 'w').write("\n".join(txt) + "\n")
            names.append(nm_)
        # slope integrand (standalone test), photon of the VP fixed = 1
        sl = ["static double slope_%s(const double *sf){" % key, "  double y[5], jac=1;"]
        free = [c for c in cls if c != [4]]
        for t, c in enumerate(free):
            for i in c: sl.append("  y[%d]=sf[%d]/%d;" % (i, t, len(c)))
            if len(c) > 1: sl.append("  jac*=sf[%d]/%d;" % (t, math.factorial(len(c)-1)))
        sl.append("  y[4]=1;")
        sl += [b.replace('REAL', 'double').replace('RABS', 'fabs') for b in body]
        sl.append("  double r=0;")
        for k in range(nt):
            sl.append("  { double iV=1/V%d; r+=T%d[0][0]*c%d*iV*iV + T%d[0][1]*c%d*iV;" % (k, k, k, k, k))
            if Ms[k] >= 2: sl.append("    r+=T%d[1][0]*iV; }" % k)
            else: sl.append("  }")
        sl.append("  return %d*r*jac;\n}" % d['weight'])
        slope.append(("\n".join(sl), len(free), key))
    tb = ['/* Generated by vp4.py */', '#include "common.h"'] + ['extern const graph_t G_%s;' % c for c in names]
    tb.append('static const graph_t *VP4GRAPHS[]={%s};' % ",".join('&G_%s' % c for c in names))
    tb.append('#define NVP4GRAPHS %d' % len(names))
    open(os.path.join(out, 'vp4table.h'), 'w').write("\n".join(tb) + "\n")
    # standalone slope check
    t = ["/* pi4slope.c: slope of the renormalized Pi^(4)(k^2) at k^2 = 0 from the generated VP terms;",
         " * expected -41/162 = -0.2530864 (one-loop analogue -1/15).  Plain Sinc sums. */",
         "#include <stdio.h>", "#include <math.h>", "#include <stdlib.h>", "#define RABS(x) fabs(x)"]
    for s_, dim, key in slope: t.append(s_)
    t.append("""static double sum(double (*f)(const double*), int D, double h, int K){
    long n=1; for(int i=0;i<D;i++) n*=2*K+1; double tot=0;
    #pragma omp parallel for reduction(+:tot) schedule(dynamic,1024)
    for(long id=0;id<n;id++){ double sf[8], w=1; long r=id;
        for(int i=0;i<D;i++){ int k=(int)(r%(2*K+1))-K; r/=(2*K+1); sf[i]=exp(k*h); w*=sf[i]*h; }
        double v=f(sf)*w; if(isfinite(v)) tot+=v; }
    return tot;
}
int main(int argc,char**argv){ double h=argc>1?atof(argv[1]):0.5; int K=argc>2?atoi(argv[2]):40;""")
    for s_, dim, key in slope:
        t.append('  double s%s=sum(slope_%s,%d,h,K); printf("slope %s: %%.10f\\n", s%s);' % (key, key, dim, key, key))
    t.append('  printf("total: %%.10f   expected -41/162 = %%.10f\\n", %s, -41.0/162); return 0; }' % "+".join('s%s' % k for _, _, k in slope))
    open(os.path.join(out, 'pi4slope.c'), 'w').write("\n".join(t) + "\n")
    pickle.dump(cache, open(cf + '.tmp', 'wb')); os.replace(cf + '.tmp', cf)
    print('VP4:', names)

if __name__ == '__main__':
    main(sys.argv[1])
