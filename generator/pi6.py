"""
Proper vacuum-polarization loops with n internal photons (a single electron loop) inserted in the
photon of the one-loop vertex.  n = 2: the sixth-order VP Pi6, Kinoshita's group I(d) (15 graphs,
Laporta set 20, 0.049513202559526).  n = 1 reproduces vp4.py (used as a check).

The loop is a cycle of 2n+2 vertices in fermion-flow order, v0 = mu (k enters), v_P = nu (k
leaves), the other 2n vertices paired by the internal photons (chords).  Lines: fermion f_i from
v_i to v_{i+1} (parameter i), chord j (parameter 2n+2+j).

Internal renormalization: Volkov's U on every UV-divergent subgraph (vertex or self-energy), which
are the cyclic intervals of the cycle with one external photon (mu, nu or one chord end: vertex)
or none (self-energy), summed over forests of nested-or-disjoint intervals with sign (-1)^|F|.
Each term is the reduced loop G/(maximal elements) (a 'PI' node, trace with mu, nu) times, per
maximal element, U of it with its own children (V0 / SE nodes from build.Node, recursively).

Overall renormalization Pi(k^2) - Pi(0).  The k^mu k^nu coefficient of each term is
(C0 + C1 k^2 + C2 k^4 + ...) e^{c k^2}; in the parent photon D(k) = int dt e^{t k^2}, D k^2 = -1:
    C0 : t -> t + c minus unshifted            main family, Gamma index n_vp - 2L_vp - w
    C1 : t = c, not integrated                  contact family, index n_vp - 2L_vp - 1 - w
    C2 : -d/dt at t = c of the parent            index n_vp - 2L_vp - 2 - w: C_m all scale alike
         (k^2 e^{t k^2} = d/dt e^{t k^2})          (lambda^(-2L-w)), and d/dt P(d) carries one more
                                                  power of the scale than the integrand needs
The derivative is a 4-point central difference in t (step 1e-3 t, relative error ~1e-12).
Loop weight: -(-1/4)^L_vp (fermion loop and loop factors, as vp4.py).

    python3 pi6.py OUTDIR            -> w4_p6*.c and pi6table_4.h (Pi6 in the one-loop vertex)
    python3 pi6.py OUTDIR check      -> pi4 check: T arrays of this code and of vp4.py (C program)
"""
import os, sys, pickle, math, itertools, re
from fractions import Fraction
import graphs as gr, build, vp4, vp4x

NW = build.NW

class Loop:
    def __init__(self, n, P, chords):
        self.n, self.P = n, P
        self.NV = 2*n + 2
        self.chords = [tuple(sorted(c)) for c in chords]
        self.partner = {}
        for j, (a, b) in enumerate(self.chords):
            self.partner[a] = (j, b); self.partner[b] = (j, a)
        self.npar = self.NV + n
    def key(self): return (self.P, tuple(sorted(self.chords)))
    def mirror(self):
        m = lambda v: (-v) % self.NV
        return Loop(self.n, m(self.P), [(m(a), m(b)) for a, b in self.chords])
    def label(self, v):
        if v == 0: return 'mu'
        if v == self.P: return 'nu'
        return 'a%d' % self.partner[v][0]
    def name(self):
        return 'P%d_' % self.P + '_'.join('%d%d' % c for c in sorted(self.chords))
    # ---- divergent subgraphs: cyclic intervals (start, length) ----
    def verts(self, iv):
        s, l = iv; return [(s + i) % self.NV for i in range(l)]
    def classify(self, iv):
        vs = self.verts(iv); S = set(vs)
        for p in (vs[0], vs[-1]):
            if p in (0, self.P) or self.partner[p][1] not in S: return None
        ext = []
        for v in vs:
            if v in (0, self.P) or self.partner[v][1] not in S: ext.append(v)
        if len(ext) > 1: return None
        off = {v: i for i, v in enumerate(vs)}
        inner = [(off[a], off[b]) for a, b in self.chords if a in S and b in S]
        inner = [tuple(sorted(c)) for c in inner]
        for i in range(len(vs) - 1):
            if not any(a <= i < b for a, b in inner): return None
        return 'SE' if not ext else 'V'
    def subgraphs(self):
        out = []
        for s in range(self.NV):
            for l in range(2, self.NV):
                k = self.classify((s, l))
                if k: out.append((s, l))
        return out
    def lines_in(self, iv):
        s, l = iv; return [(s + i) % self.NV for i in range(l - 1)]
    def chords_in(self, iv):
        S = set(self.verts(iv)); return [j for j, (a, b) in enumerate(self.chords) if a in S and b in S]

def pairings(pts):
    if not pts: yield (); return
    a = pts[0]
    for i in range(1, len(pts)):
        for m in pairings(pts[1:i] + pts[i+1:]): yield ((a, pts[i]),) + m

def loops(n):
    """One representative per mirror pair, with weight."""
    NV = 2*n + 2; seen = set(); out = []
    for P in range(1, NV):
        for pr in pairings([v for v in range(1, NV) if v != P]):
            L = Loop(n, P, pr)
            if L.key() in seen: continue
            m = L.mirror(); seen.add(L.key()); seen.add(m.key())
            out.append((L, 1 if m.key() == L.key() else 2))
    return out

def nested(a, b, L):   # a contains b
    return set(L.verts(b)) <= set(L.verts(a))
def compatible(a, b, L):
    A, B = set(L.verts(a)), set(L.verts(b))
    return A <= B or B <= A or not (A & B)

def forests(L):
    subs = L.subgraphs(); out = []
    for r in range(len(subs) + 1):
        for comb in itertools.combinations(subs, r):
            if all(compatible(a, b, L) for a, b in itertools.combinations(comb, 2)): out.append(comb)
    return out

def children(iv, F, L):
    inner = [s for s in F if s != iv and nested(iv, s, L)]
    return [s for s in inner if not any(t != s and nested(t, s, L) for t in inner)]

def maximal(F, L):
    return [s for s in F if not any(t != s and nested(t, s, L) for t in F)]

# ---------------- inner (open-line) nodes via build.Node ----------------
def inner_node(L, iv, kids, tag):
    """U of the interval iv with children kids: a build.Node on the open line obtained by cutting
    the loop at the external vertex not in iv, parameters renamed to the loop's."""
    vs = set(L.verts(iv))
    cut = L.P if 0 in vs else 0
    other = 0 if cut == L.P else L.P
    pos = {(cut + 1 + i) % L.NV: i for i in range(L.NV - 1)}      # vertex -> position
    arcs = [(pos[a], pos[b]) for a, b in L.chords]
    g = gr.Graph(L.n, pos[other], arcs)
    seg = lambda x: (pos[L.verts(x)[0]], pos[L.verts(x)[-1]])
    s0 = seg(iv); assert s0[0] < s0[1] and g.segment(*s0) is not None, (L.key(), iv)
    ks = tuple(sorted(seg(k) for k in kids))
    kind = 'SE' if L.classify(iv) == 'SE' else 'V0'
    nd = build.Node(g, s0, ks, kind)
    # parameter renaming: g electron i (positions i -> i+1) is fermion f_{vertex(i)}; g photon 2n+k
    vert = {i: v for v, i in pos.items()}
    gch = [L.chords.index(tuple(sorted((vert[a], vert[b])))) for a, b in g.arcs]
    ren = lambda i: vert[i] if i < 2*L.n else L.NV + gch[i - 2*L.n]
    nd.gidx = [ren(i) for i in nd.gidx]
    nd._pos = pos
    nd.name = 'u%s_%s' % (tag, '_'.join('%d%d' % (x[0], x[1]) for x in [iv] + sorted(kids)))
    return nd

# ---------------- reduced loop (PI node) ----------------
def pi_node(L, F, tag):
    """The loop with the maximal elements of F shrunk: V -> a vertex, SE -> an insertion."""
    mx = maximal(F, L)
    owner = {}
    for iv in mx:
        for v in L.verts(iv): owner[v] = iv
    # flow walk from v0 over reduced nodes; items: ('n', label, members) or ('se', iv)
    items = []; v = L.verts(owner[0])[0] if 0 in owner else 0; start = v
    while True:
        iv = owner.get(v)
        if iv is None:
            items.append(('n', L.label(v), [v])); nxt = (v + 1) % L.NV
        else:
            vs = L.verts(iv)
            if L.classify(iv) == 'V':
                ext = [u for u in vs if u in (0, L.P) or L.partner[u][1] not in set(vs)][0]
                items.append(('n', L.label(ext), vs))
            else:
                items.append(('se', iv, vs))
            nxt = (vs[-1] + 1) % L.NV
        if nxt == start: break
        v = nxt
    # lines: between consecutive reduced items the fermion line leaving the last vertex of the item
    lines = []                      # global fermion index, in flow order
    for it in items: lines.append(it[-1][-1])       # line f_{last vertex}
    # chords remaining
    inside = set(j for iv in mx for j in L.chords_in(iv))
    rem = [j for j in range(L.n) if j not in inside]
    # reduced node index of each vertex
    nodeidx = {}
    for t, it in enumerate(items):
        for u in it[-1]: nodeidx[u] = t
    # momentum: line after item t (t = 0..len-1) ; k from node of v0 (t=0) to node of nu
    tnu = nodeidx[L.P]
    A = []; B = []; m2 = []; gidx = []
    for t, f in enumerate(lines):
        row = [1]
        for j in rem:
            a, b = sorted((nodeidx[L.chords[j][0]], nodeidx[L.chords[j][1]]))
            row.append(1 if a <= t < b else 0)
        A.append(row); B.append([1 if t < tnu else 0, 0]); m2.append(1); gidx.append(f)
    for jj, j in enumerate(rem):
        A.append([0] + [1 if kk == jj else 0 for kk in range(len(rem))]); B.append([0, 0]); m2.append(0); gidx.append(L.NV + j)
    # loop tokens in matrix order (reverse flow); slots: lines 0..nl-1 then SE insertions
    nl = len(lines); selist = []; slotline = list(range(nl)); flow = []
    for t, it in enumerate(items):
        if it[0] == 'n': flow.append(('v', it[1]))
        else:
            selist.append(it[1]); slotline.append(t - 1 if t > 0 else nl - 1)
            flow.append(('s', nl + len(selist) - 1))
        flow.append(('s', t))
    # an SE item sits between line t-1 (into it) and line t (out of it): same momentum
    loop = list(reversed(flow))
    name = 'pi%s_%s' % (tag, '_'.join('%d%d' % x for x in sorted(mx)) or 'x')
    nd = vp4.mk(name, 'PI', A, B, m2, gidx, None, nl + len(selist), slotline=slotline, selist=selist, loop=loop)
    return nd, mx

def terms(L, tag):
    """[(sign, pinode, maximal elements, {iv: inner node})]"""
    out = []
    for fi, F in enumerate(forests(L)):
        pn, mx = pi_node(L, F, tag + 'f%d' % fi)
        nodes = {}
        def mknode(iv):
            ks = children(iv, F, L)
            nodes[iv] = inner_node(L, iv, ks, tag + 'f%d' % fi)
            for k in ks: mknode(k)
        for iv in mx: mknode(iv)
        out.append(((-1)**len(F), pn, mx, nodes, F))
    return out

def series_classes(L, T):
    """Fermion lines with identical momentum rows in the full loop, split if some term separates
    them; chords singletons."""
    rows = {}
    for f in range(L.NV):
        r = tuple([1 if a <= f < b else 0 for a, b in L.chords] + [1 if f < L.P else 0])
        rows.setdefault(r, []).append(f)
    cls = []
    for c in rows.values():
        parts = {}
        for f in c:
            sig = []
            for sg, pn, mx, nodes, F in T:
                where = [nd.name for nd in [pn] + list(nodes.values()) if f in nd.gidx]
                sig.append(tuple(where))
            parts.setdefault(tuple(sig), []).append(f)
        cls += list(parts.values())
    return cls + [[L.NV + j] for j in range(L.n)]

def gen_terms_code(L, T, tabs, K):
    """C code: for each term k, T<k>[m][w] (with K and the term sign), V<k>, c<k>."""
    o = []; done = set()
    for sg, pn, mx, nodes, F in T:
        if pn.name not in done: o.append(vp4.emit_pi(pn, tabs[pn.name])); done.add(pn.name)
        for nd in nodes.values():
            if nd.name not in done: o.append(vp4.emit_inner(nd, tabs[nd.name])); done.add(nd.name)
    em = build.Emitter(); Ms = []
    def value(iv, nodes):
        nd = nodes[iv]
        kv = {c: value(c, nodes) for c in nd.kids_iv}
        outs = [0] if nd.kind == 'V0' else [0, 1]
        res = []
        for k in outs:
            dst = em.new()
            for pat in range(2**len(nd.selist)):
                cur, cura = '%s_W[%d][%d]' % (nd.name, k, pat), '%s_Wa[%d][%d]' % (nd.name, k, pat)
                for t, c in enumerate(nd.se_iv):
                    sel = kv[c][1] if (pat >> t) & 1 else kv[c][0]
                    nxt = em.new(); em.conv_into(nxt, cur, sel, cura, sel + 'a'); cur, cura = nxt, nxt + 'a'
                for c in nd.v_iv:
                    sel = kv[c][0]
                    nxt = em.new(); em.conv_into(nxt, cur, sel, cura, sel + 'a'); cur, cura = nxt, nxt + 'a'
                em.lines.append("  for(int w=0;w<%d;w++){ %s[w]+=%s[w]; %sa[w]+=%s[w]; }" % (NW, dst, cur, dst, cura))
            res.append(dst)
        return res
    for k, (sg, pn, mx, nodes, F) in enumerate(T):
        M = pn.M; Ms.append(M)
        em.lines.append("  REAL T%d[%d][%d]={{0}}, T%da[%d][%d]={{0}}, V%d, c%d;" % (k, M, NW, k, M, NW, k, k))
        vals = {iv: value(iv, nodes) for iv in mx}
        f = build.cnum(K*sg)
        for m in range(M):
            for pat in range(2**len(pn.selist)):
                cur, cura = '%s_W[%d][%d]' % (pn.name, m, pat), '%s_Wa[%d][%d]' % (pn.name, m, pat)
                for t, c in enumerate(pn.selist):
                    sel = vals[c][1] if (pat >> t) & 1 else vals[c][0]
                    nxt = em.new(); em.conv_into(nxt, cur, sel, cura, sel + 'a'); cur, cura = nxt, nxt + 'a'
                for c in mx:
                    if L.classify(c) != 'V': continue
                    sel = vals[c][0]
                    nxt = em.new(); em.conv_into(nxt, cur, sel, cura, sel + 'a'); cur, cura = nxt, nxt + 'a'
                em.lines.append("  for(int w=0;w<%d;w++){ T%d[%d][w]+=%s*%s[w]; T%da[%d][w]+=RABS(%s)*%s[w]; }" % (NW, k, m, f, cur, k, m, f, cura))
        Vs = "+".join(['%s_V' % pn.name] + ['%s_V' % nd.name for nd in nodes.values()])
        em.lines.append("  V%d=%s; c%d=%s_c;" % (k, Vs, k, pn.name))
    return "\n".join(o) + "\n" + "\n".join(em.lines), Ms

def prepare(L, tag, cache):
    T = terms(L, tag)
    tabs = {}
    for sg, pn, mx, nodes, F in T:
        tabs[pn.name] = vp4.table(pn, cache)
        for iv, nd in nodes.items():
            # children bookkeeping in loop coordinates for value()
            ks = children(iv, F, L)
            nd.se_iv = [c for c in ks if L.classify(c) == 'SE']
            nd.v_iv = [c for c in ks if L.classify(c) == 'V']
            nd.kids_iv = ks
            # build.Node orders selist by segment position; match ours to it
            order = sorted(nd.se_iv, key=lambda c: nd.selist.index(_seg_of(nd, L, c)))
            nd.se_iv = order
            tabs[nd.name] = nd.table(cache)
    return T, tabs

def _seg_of(nd, L, c):
    """Segment (in the node's open-line positions) of loop interval c."""
    vs = L.verts(c)
    return (nd._pos[vs[0]], nd._pos[vs[-1]])

def q_function(fname, L, T, tabs, K, ny):
    """Test function: sum_k sum_m sum_w T[m][w] kap^m beta^w exp(-V + c kap)."""
    code, Ms = gen_terms_code(L, T, tabs, K)
    o = ["static double %s(const double *y, double kap, double beta){" % fname]
    o.append("  const double %s;" % ", ".join("y%d=y[%d]" % (i, i) for i in range(ny)))
    o.append(code)
    o.append("  double r=0;")
    for k, M in enumerate(Ms):
        o.append("  for(int m=0;m<%d;m++) for(int w=0;w<%d;w++) r+=T%d[m][w]*pow(kap,m)*pow(beta,w)*exp(-V%d+c%d*kap);" % (M, NW, k, k, k))
    o.append("  return r;\n}")
    return "\n".join(o)

def check(out):
    """Pi4 through this code against vp4.py, pointwise."""
    cf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache_dirac.pkl')
    cache = pickle.load(open(cf, 'rb'))
    D = vp4.diagrams()
    src = ['#include <stdio.h>', '#include <math.h>', '#include <stdlib.h>', '#define REAL double', '#define RABS(x) fabs(x)', '#define FN(x) x##_d']
    pairs = []
    for key, L in (('V', Loop(1, 2, [(1, 3)])), ('S', Loop(1, 3, [(1, 2)]))):
        T, tabs = prepare(L, 'c' + key, cache)
        src.append(q_function('mine_' + key, L, T, tabs, Fraction(-1, 16), 5))
        d = D[key]; t2 = {n.name: vp4.table(n, cache) for n in d['nodes']}
        code, Ms = vp4.gen_terms_code(d, t2)
        o = ["static double ref_%s(const double *y, double kap, double beta){" % key,
             "  const double %s;" % ", ".join("y%d=y[%d]" % (i, i) for i in range(5)), code, "  double r=0;"]
        for k, M in enumerate(Ms):
            o.append("  for(int m=0;m<%d;m++) for(int w=0;w<%d;w++) r+=T%d[m][w]*pow(kap,m)*pow(beta,w)*exp(-V%d+c%d*kap);" % (M, NW, k, k, k))
        o.append("  return r;\n}")
        src.append("\n".join(o)); pairs.append(key)
    src.append("int main(){ srand(1); for(int t=0;t<5;t++){ double y[5]; for(int i=0;i<5;i++) y[i]=0.2+2.0*rand()/RAND_MAX;")
    for key in pairs:
        src.append('  printf("%s  mine %%+.15e  ref %%+.15e\\n", mine_%s(y,-0.7,0.3), ref_%s(y,-0.7,0.3));' % (key, key, key))
    src.append("} return 0; }")
    open(os.path.join(out, 'pi6check.c'), 'w').write("\n".join(src) + "\n")

def point_fn(L, T, tabs, K, weight, cls, cname, base, family):
    """FN(p_<cname>)(x, v, &ab): x = parent parameters (photon 2), v = loop class variables.
    family 0: C0, shifted minus unshifted; 1: C1 contact; 2: C2 contact, -d/dt."""
    code, Ms = gen_terms_code(L, T, tabs, K)
    nt = len(T); nv = len(cls); ny = L.npar
    o = ["static REAL FN(p_%s)(const REAL *x, const REAL *v, REAL *ab){" % cname]
    o.append("  REAL y[%d], jac=1;" % ny)
    for t, c in enumerate(cls):
        for i in c: o.append("  y[%d]=v[%d]/%d;" % (i, t, len(c)))
        if len(c) > 1: o.append("  jac*=%s;" % "*".join(["v[%d]" % t]*(len(c)-1) + ["((REAL)1/%d)" % math.factorial(len(c)-1)]))
    o.append("  const REAL %s;" % ", ".join("y%d=y[%d]" % (i, i) for i in range(ny)))
    o.append(code)
    o.append("  const REAL om=jac; REAL s[3], tot=0, tota=0, t, ta;")      # weight: graph_t field
    p0 = L.npar - 2*(L.n + 1)
    for k in range(nt):
        if family > 0 and Ms[k] <= family: continue
        m = family
        o.append("  for(int w=0;w<%d;w++){ const REAL Tw=T%d[%d][w], Taw=T%da[%d][w]; if(Tw==0 && Taw==0) continue;" % (NW, k, m, k, m))
        o.append("    s[0]=x[0]; s[1]=x[1];")
        if family == 0:
            o.append("    s[2]=x[2]+c%d; t=FN(pt_%s)(s,V%d,%d-w,&ta); tot+=Tw*t; tota+=(Taw+RABS(Tw))*ta;" % (k, base, k, p0))
            o.append("    s[2]=x[2];     t=FN(pt_%s)(s,V%d,%d-w,&ta); tot-=Tw*t; tota+=(Taw+RABS(Tw))*ta;" % (base, k, p0))
        elif family == 1:
            o.append("    s[2]=c%d; t=FN(pt_%s)(s,V%d,%d-w,&ta); tot-=Tw*t; tota+=(Taw+RABS(Tw))*ta;" % (k, base, k, p0 - 1))
        else:
            # -d/dt P(t) at t = c: (-P(t+2e)+8P(t+e)-8P(t-e)+P(t-2e))/(12e)
            o.append("    const REAL e=c%d*(REAL)1e-3; REAL d=0, da=0;" % k)
            for mult, wt in ((2, -1), (1, 8), (-1, -8), (-2, 1)):
                o.append("    s[2]=c%d+(%d)*e; t=FN(pt_%s)(s,V%d,%d-w,&ta); d+=(%d)*t; da+=%d*ta;" % (k, mult, base, k, p0 - 2, wt, abs(wt)))
            o.append("    d/=12*e; da/=12*e; tot-=Tw*d; tota+=(Taw+RABS(Tw))*da;")
        o.append("  }")
    o.append("  *ab=RABS(om)*tota; return om*tot;\n}")
    return "\n".join(o), nv, Ms

def emit(L, w0, cname, cache, out, g1, pt, base):
    T, tabs = prepare(L, cname, cache)
    cls = series_classes(L, T)
    K = -Fraction(-1, 4)**(L.n + 1)
    names = []
    _, _, Ms = point_fn(L, T, tabs, K, w0, cls, cname + 'm', base, 0)
    fams = [0] + [f for f in (1, 2) if any(M > f for M in Ms)]
    for fam in fams:
        cn = '%s%s' % (cname, 'mcd'[fam])
        pf, nv, _ = point_fn(L, T, tabs, K, w0, cls, cn, base, fam)
        if fam == 0: wr, D, ngl, glc = build.gen_wrapper(g1, cn, inner='p_%s' % cn, nextra=nv)
        else: wr, D, ngl, glc = vp4x.wrapper_contact(g1, cn, 'p_%s' % cn, nv, 0)
        tt = 'VP: 1-loop vertex + Pi%d(%s) %s' % (2*L.n + 2, L.name(), ('main', 'contact', 'contact k^4')[fam])
        src = ["/* Generated by pi6.py: %s  (weight %d, D=%d) */" % (tt, w0, D),
               '#include "common.h"', glc,
               '#define REAL double', '#define FN(x) x##_d', pt, pf, wr, '#undef REAL', '#undef FN',
               '#define REAL QREAL', '#define FN(x) x##_q', pt, pf, wr, '#undef REAL', '#undef FN',
               'const graph_t G_%s={"%s","%s",%d,%d,%d,%d,%d,f_%s_d,f_%s_q,0,0};' %
               (cn, cn, tt, L.n + 2, D, w0, g1.n, 2*g1.L, cn, cn)]
        open(os.path.join(out, cn + '.c'), 'w').write("\n".join(src) + "\n")
        names.append(cn); print(cn, w0, D, tt, 'terms', len(T), 'classes', cls, flush=True)
    return names

def main(out, n=2):
    cf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache_dirac.pkl')
    cache = pickle.load(open(cf, 'rb'))
    (g1, w1), = gr.enumerate_graphs(1)
    pre = 'w4_p6' if n == 2 else 'w3_p4'
    base = pre + '_b'
    pt, info = build.gen_graph(g1, base, cache, dpvar=True)
    names = []
    for i, (L, w) in enumerate(loops(n)):
        names += emit(L, w, '%s%d' % (pre, i), cache, out, g1, pt, base)
    Ltot = n + 2
    tb = ['/* Generated by pi6.py */', '#include "common.h"'] + ['extern const graph_t G_%s;' % c for c in names]
    tb.append('static const graph_t *PI6GRAPHS_%d[]={%s};' % (Ltot, ",".join('&G_%s' % c for c in names)))
    tb.append('#define NPI6GRAPHS_%d %d' % (Ltot, len(names)))
    open(os.path.join(out, 'pi6table_%d.h' % Ltot), 'w').write("\n".join(tb) + "\n")
    pickle.dump(cache, open(cf + '.tmp', 'wb')); os.replace(cf + '.tmp', cf)
    print('Pi%d insertions:' % (2*n + 2), names)

if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[2] == 'check': check(sys.argv[1])
    elif len(sys.argv) > 2 and sys.argv[2] == 'pi4': main(sys.argv[1], 1)
    else: main(sys.argv[1])
