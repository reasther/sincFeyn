"""
Build the reduced-graph nodes of every subtraction term of a graph and emit C.

Every term of Volkov's formula is a product over the forest elements G_j of an operator applied
to the reduced graph G_j / (maximal forest elements inside G_j).  In Schwinger parameters each
factor is a one- or two-loop parametric integrand in its own lines; the product lives in the
parameter space of G and shares the overall scale, so the overall scale integral gives
Gamma(n-2L-w) / (sum_j V_j)^(n-2L-w) with w the total number of Wick pairs.

Node kinds:  F2  (operator A: F2(0) via the t-derivative of the projector numerator)
             V0  (vertex at q=0 on shell: U -> a, L -> a+b+c)
             SE  (self-energy on shell: U -> u + v pslash, inserted into the parent chain)
"""
import sympy as sp, re, math, pickle, os, hashlib
from fractions import Fraction
import graphs as gr, dirac
from symb import symbolic

NW = 4    # graded arrays carry Wick counts 0..NW-1

class Node:
    def __init__(self, g, seg, kids, kind):
        self.g, self.seg, self.kids, self.kind = g, seg, tuple(kids), kind
        u, v = seg
        inkid = lambda pos: next((c for c in kids if c[0] <= pos <= c[1]), None)
        gapinkid = lambda i: any(c[0] <= i < c[1] for c in kids)
        karcs = set(k for c in kids for k in g.segment(*c)['arcs'])
        self.arcs = [k for k in g.segment(u, v)['arcs'] if k not in karcs]
        self.egaps = [i for i in range(u, v) if not gapinkid(i)]
        self.selist = [c for c in kids if g.segment(*c)['kind'] == 'SE']
        self.vlist = [c for c in kids if g.segment(*c)['kind'] == 'V']
        # slots: electron lines then SE insertions
        self.slots = [('e', i) for i in self.egaps] + [('x', c) for c in self.selist]
        self.gidx = [i for i in self.egaps] + [2*g.L + k for k in self.arcs]   # global params of node lines
        # chain (matrix order: from v down to u)
        def phname(pos):
            if pos == g.ext: return 'mu'
            k, q = g.partner[pos]
            return 'a%d' % k if k in self.arcs else 'mu'
        chain = []; pos = v
        while pos >= u:
            c = inkid(pos)
            if c is None: chain.append(('v', phname(pos))); lo = pos
            else:
                sc = g.segment(*c)
                if sc['kind'] == 'V': chain.append(('v', phname(sc['extpos'])))
                else: chain.append(('s', len(self.egaps) + self.selist.index(c)))
                lo = c[0]
            if lo > u: chain.append(('s', self.egaps.index(lo-1)))
            pos = lo - 1
        self.chain = chain
        nmu = sum(1 for t in chain if t == ('v', 'mu'))
        assert nmu == (0 if kind == 'SE' else 1), (seg, kids, kind, chain)
        # momentum routing: loop momentum per node arc; electron gap i covered by arc (a,b) if a<=i<b
        A = []; B = []; m2 = []
        for i in self.egaps:
            A.append([-1 if g.arcs[k][0] <= i < g.arcs[k][1] else 0 for k in self.arcs])
            B.append([1, 0] if (kind != 'F2' or i < g.ext) else [0, 1]); m2.append(1)
        for k in self.arcs:
            A.append([1 if kk == k else 0 for kk in self.arcs]); B.append([0, 0]); m2.append(0)
        self.A, self.B, self.m2 = A, B, m2
        self.L = len(self.arcs)
        # slot -> local line index providing its momentum coefficients
        self.slotline = list(range(len(self.egaps)))
        for c in self.selist:
            self.slotline.append(self.egaps.index(c[0]-1))
        self.name = 'n_%s_%d_%d_%s' % (kind, u, v, '_'.join('%d%d' % c for c in kids) or 'x')

    def table(self, cache):
        key = repr((self.chain, len(self.slots), self.kind)) + (repr(self.loop) if getattr(self, 'loop', None) else '')
        h = hashlib.md5(key.encode()).hexdigest()
        if h not in cache:
            cache[h] = dirac.table(self.chain, len(self.slots), self.kind, loop=getattr(self, 'loop', None))
        return cache[h]

def cexpr(e):
    s = sp.ccode(e)
    s = re.sub(r'pow\((x\d+), (\d+)\)', lambda m: '(' + '*'.join([m.group(1)]*int(m.group(2))) + ')', s)
    assert 'pow' not in s, s
    return s

def cnum(f):
    f = Fraction(f)
    if f.denominator == 1: return '%d' % f.numerator
    return '((REAL)%d/%d)' % (f.numerator, f.denominator)

CHUNK = int(os.environ.get('GEN_CHUNK', '20000'))   # max characters per factored expression before it becomes its own function (0 = off)
TRIE = True      # emit node polynomials in factored form (False: flat sum, one line per structure)

def _factored(terms):
    """terms: list of (Fraction coefficient, list of factor names).  Returns (expr, abs-bound expr)."""
    import collections
    const = sum((c for c, f in terms if not f), Fraction(0))
    rest = [(c, f) for c, f in terms if f]
    parts, aparts = [], []
    if const != 0: parts.append(cnum(const)); aparts.append(cnum(abs(const)))
    while rest:
        cnt = collections.Counter(x for _, f in rest for x in set(f))
        best = max(cnt.items(), key=lambda kv: (kv[1], kv[0]))[0]   # deterministic tie-break
        inn = []; out = []
        for c, f in rest:
            if best in f:
                g = list(f); g.remove(best); inn.append((c, g))
            else: out.append((c, f))
        e, ea = _factored(inn)
        parts.append("%s*(%s)" % (best, e)); aparts.append("A%s*(%s)" % (best, ea))
        rest = out
    if not parts: return "0", "0"
    return "+".join(parts), "+".join(aparts)

class Chunker:
    """Splits factored expressions longer than `limit` characters into separate C functions
    FN(ch_<prefix><j>)(F, r) (r[0] = value, r[1] = abs bound) reading the factors from an array F,
    so that no single function is too large for the compiler (gcc's memory grows with the largest
    function)."""
    def __init__(self, prefix, limit=20000):
        self.prefix, self.limit, self.funcs, self.n = prefix, limit, [], 0
    def chunk(self, e, ea, pend):
        j = self.n; self.n += 1; nm = 'ch_%s%d' % (self.prefix, j)
        import re
        used = sorted(set(re.findall(r'\bA?(?:cpp|cp|cP|w)\d+(?:_\d+)?\b', e + ' ' + ea)))
        o = ["static void FN(%s)(const REAL *F, REAL *r){" % nm]
        if used: o.append("  const REAL %s;" % ", ".join("%s=F[%d]" % (u, self.index[u]) for u in used))
        for p in pend: o.append("  REAL c%s[2]; FN(ch_%s%d)(F, c%s);" % (p, self.prefix, p, p))
        o.append("  r[0]=%s; r[1]=%s;" % (e, ea)); o.append("}")
        self.funcs.append("\n".join(o))
        return 'c%d[0]' % j, 'c%d[1]' % j, [j]

def _fac2(terms, ch):
    """As _factored, but subexpressions longer than ch.limit become chunk functions.
    Returns (expr, abs expr, chunk ids to be computed before use)."""
    import collections
    const = sum((c for c, f in terms if not f), Fraction(0))
    rest = [(c, f) for c, f in terms if f]
    parts, aparts, pend = [], [], []
    if const != 0: parts.append(cnum(const)); aparts.append(cnum(abs(const)))
    while rest:
        cnt = collections.Counter(x for _, f in rest for x in set(f))
        best = max(cnt.items(), key=lambda kv: (kv[1], kv[0]))[0]
        inn = []; out = []
        for c, f in rest:
            if best in f:
                g = list(f); g.remove(best); inn.append((c, g))
            else: out.append((c, f))
        e, ea, p = _fac2(inn, ch)
        parts.append("%s*(%s)" % (best, e)); aparts.append("A%s*(%s)" % (best, ea)); pend += p
        rest = out
    if not parts: return "0", "0", []
    e, ea = "+".join(parts), "+".join(aparts)
    if len(e) > ch.limit: return ch.chunk(e, ea, pend)
    return e, ea, pend

# Line identifications.  The coefficients of a node are functions of the momentum routing only:
# a line's shift coefficients cp, cpp (cP = cp + cpp) are B_i + A_i.X, its contraction coefficients
# w_ij = A_i adj(M) A_j / U, with (A_i, B_i) its loop and external rows.  Lines with the same
# (A, B) rows (one series class) therefore have identical cp, cpp, cP and w, and lines with the same
# A row (e.g. the two lines at the external vertex) identical w; with IDENT=2 also cp_i = cp_r +
# (B_i - B_r) for the first line r of the same A row.  Writing the Dirac polynomial in one symbol
# per distinct coefficient lets FORM (and the factored tree) merge the monomials before optimising.
# GEN_IDENT=0 keeps one symbol per slot (the kernels before 2026-10-25), 1 only the identities.
IDENT = None
def _ident():
    if IDENT is not None: return IDENT
    return int(os.environ.get('GEN_IDENT', '2'))

def _reps(nd):
    if getattr(nd, '_rep', None) is None:
        ra, rab = {}, {}; nd._rep = ([], [])
        for i in range(len(nd.A)):
            nd._rep[0].append(ra.setdefault(tuple(nd.A[i]), i))
            nd._rep[1].append(rab.setdefault((tuple(nd.A[i]), tuple(nd.B[i])), i))
    return nd._rep

def _monomials(nd, combo, pr):
    """The factors of one table structure, as a list of (Fraction multiplier, [factor names]):
    one entry, unless IDENT=2 expands a shift coefficient into its representative plus a constant."""
    lev = _ident()
    ra, rab = _reps(nd) if lev else (None, None)
    lin = []                                     # each factor: list of (coefficient, name or None)
    for j, ch in enumerate(combo):
        if ch not in ('p', 'pp', 'P'): continue
        i = nd.slotline[j]
        if not lev: lin.append([(1, 'c%s%d' % (ch, i))]); continue
        if lev == 1: lin.append([(1, 'c%s%d' % (ch, rab[i]))]); continue
        r = ra[i]; Bi, Br = nd.B[i], nd.B[r]
        d = {'p': Bi[0] - Br[0], 'pp': Bi[1] - Br[1], 'P': Bi[0] + Bi[1] - Br[0] - Br[1]}[ch]
        lin.append([(1, 'c%s%d' % (ch, r))] + ([(d, None)] if d else []))
    fac = []
    for a, b in pr:
        i, j = nd.slotline[a], nd.slotline[b]
        if lev: i, j = ra[i], ra[j]
        i, j = sorted((i, j)); fac.append('w%d_%d' % (i, j))
    out = [(Fraction(1), fac)]
    for opts in lin:
        out = [(m*c, f + ([nm] if nm else [])) for m, f in out for c, nm in opts]
    return out

def _emit_node_py(nd, tab, chk=None):
    """C code block computing node arrays: <nm>_U, <nm>_V, <nm>_W[kind][pat][w], <nm>_Wa[...].
    chk: a Chunker; very long expressions then become separate functions (ch.funcs)."""
    s, U, F, coefU, CU = symbolic(nd.A, nd.B, nd.m2)
    sub = {s[i]: sp.Symbol('x%d' % nd.gidx[i]) for i in range(len(s))}
    cc = lambda e: cexpr(sp.expand(sp.sympify(e).subs(sub)))
    nm = nd.name; nk = 1 if nd.kind == 'F2' else 2; npat = 2**len(nd.selist)
    o = ["  /* node %s: seg %s kids %s  L=%d */" % (nd.kind, nd.seg, nd.kids, nd.L)]
    o.append("  REAL %s_W[%d][%d][%d]={{{0}}}, %s_Wa[%d][%d][%d]={{{0}}}, %s_U, %s_V;" % (nm, nk, npat, NW, nm, nk, npat, NW, nm, nm))
    o.append("  { const REAL U=%s, iU=1/U;" % cc(U))
    o.append("    const REAL V=(%s)*iU;" % cc(F))
    nl = len(nd.egaps)
    used = set(); pairs = set()
    for (combo, pr) in tab:
        for _, fac in _monomials(nd, combo, pr):
            for f in fac:
                m = re.match(r'c(pp|p|P)(\d+)$', f)
                if m: used.add((int(m.group(2)), m.group(1)))
                else: pairs.add(tuple(map(int, f[1:].split('_'))))
    for (i, ch) in sorted(used):
        if ch == 'p': e = coefU[i][0]
        elif ch == 'pp': e = coefU[i][1]
        else: e = coefU[i][0] + coefU[i][1]
        o.append("    const REAL c%s%d=(%s)*iU;" % (ch, i, cc(e)))
    for (i, j) in sorted(pairs):
        o.append("    const REAL w%d_%d=(%s)*iU;" % (i, j, cc(CU[i][j])))
    maxw = 0
    groups = {}                      # (k, pat, w) -> list of (coefficient, [factors])
    for (combo, pr), val in tab.items():
        w = len(pr); maxw = max(maxw, w)
        pat = 0
        for t, c in enumerate(nd.selist):
            if combo[nl + t] != '1': pat |= 1 << t
        pw = Fraction(-1, 2)**w
        vals = [val*pw] if nd.kind == 'F2' else [x*pw for x in val]
        for m, fac in _monomials(nd, combo, pr):
            for k, x in enumerate(vals):
                if x != 0: groups.setdefault((k, pat, w), []).append((x*m, fac))
    if _ident():                     # identified monomials summed first (exact; also the |.| bound)
        for key, terms in groups.items():
            d = {}
            for x, fac in terms: fs = tuple(sorted(fac)); d[fs] = d.get(fs, 0) + x
            groups[key] = [(x, list(fs)) for fs, x in d.items() if x != 0]
    if not TRIE:
        o.append("    REAL q;")
        for (k, pat, w), terms in groups.items():
            for x, fac in terms:
                o.append("    q=%s;" % ('*'.join(fac) if fac else '1'))
                o.append("    %s_W[%d][%d][%d]+=%s*q; %s_Wa[%d][%d][%d]+=RABS(%s*q);" % (nm, k, pat, w, cnum(x), nm, k, pat, w, cnum(x)))
    else:
        # factored (Horner-like) form: pull out the most common factor recursively.  The bound
        # sum |c * prod f| is computed on the same tree with |f| and |c| (exactly equal).
        facs = sorted({f for t in groups.values() for _, fa in t for f in fa})
        o.append("    const REAL %s;" % ", ".join("A%s=RABS(%s)" % (f, f) for f in facs) if groups else "")
        if chk is not None and groups:
            names = facs + ['A' + f for f in facs]
            chk.index = {u: i for i, u in enumerate(names)}
            o.append("    const REAL F_[%d]={%s};" % (len(names), ",".join(names)))
        for (k, pat, w), terms in groups.items():
            if chk is None:
                e, ea = _factored(terms)
            else:
                e, ea, pend = _fac2(terms, chk)
                for p in pend: o.append("    REAL c%d[2]; FN(ch_%s%d)(F_, c%d);" % (p, chk.prefix, p, p))
            o.append("    %s_W[%d][%d][%d]+=%s; %s_Wa[%d][%d][%d]+=%s;" % (nm, k, pat, w, e, nm, k, pat, w, ea))
    o.append("    const REAL iU2=iU*iU;")
    o.append("    for(int k=0;k<%d;k++) for(int p=0;p<%d;p++) for(int w=0;w<%d;w++){ %s_W[k][p][w]*=iU2; %s_Wa[k][p][w]*=iU2; }" % (nk, npat, NW, nm, nm))
    o.append("    %s_U=U; %s_V=V; }" % (nm, nm))
    nd.maxw = maxw
    return "\n".join(o)

class Emitter:
    """Generates the graded-array algebra for one graph's terms."""
    def __init__(self):
        self.lines = []; self.cnt = 0; self.memo = {}
    def new(self):
        self.cnt += 1; v = 't%d' % self.cnt
        self.lines.append("  REAL %s[%d]={0}, %sa[%d]={0};" % (v, NW, v, NW)); return v
    def conv_into(self, dst, a, b, aa, ba):
        # dst += a (*) b  (graded, truncated at NW)
        for i in range(NW):
            for j in range(NW - i):
                self.lines.append("  %s[%d]+=%s[%d]*%s[%d]; %sa[%d]+=%s[%d]*%s[%d];" % (dst, i+j, a, i, b, j, dst, i+j, aa, i, ba, j))

def node_function(nd, tab, prefix=''):
    """Node as its own C function FN(nd_<prefix><name>)(x, W, Wa, V): the compiler then never sees one
    huge function (gcc's peak memory scales with the largest function, not with the file).
    Returns (function code, declaration + call code for the caller)."""
    import re
    nm = nd.name; fn = 'nd_%s%s' % (prefix, nm)
    ch = Chunker('%s%s_' % (prefix, nm), CHUNK) if CHUNK else None
    body = emit_node(nd, tab, ch)
    nk = 1 if nd.kind == 'F2' else 2; npat = 2**len(nd.selist); sz = nk*npat*NW
    f = (ch.funcs if ch else []) + ["static void FN(%s)(const REAL *x, REAL *W, REAL *Wa, REAL *V){" % fn]
    used = sorted({int(t) for t in re.findall(r'\bx(\d+)\b', body)})
    if used: f.append("  const REAL %s;" % ", ".join("x%d=x[%d]" % (i, i) for i in used))
    f.append(body)
    f.append("  { const REAL *w=&%s_W[0][0][0], *wa=&%s_Wa[0][0][0]; for(int i=0;i<%d;i++){ W[i]=w[i]; Wa[i]=wa[i]; } *V=%s_V; }" % (nm, nm, sz, nm))
    f.append("}")
    call = ["  REAL %s_W[%d][%d][%d], %s_Wa[%d][%d][%d], %s_V;" % (nm, nk, npat, NW, nm, nk, npat, NW, nm),
            "  FN(%s)(x, &%s_W[0][0][0], &%s_Wa[0][0][0], &%s_V);" % (fn, nm, nm, nm)]
    return "\n".join(f), "\n".join(call)

def gen_graph(g, cname, cache, log=print, dpvar=False):
    """dpvar: the point function takes an extra int dp added to the Gamma-function index
    (used when the graph is a factor of a larger parameter space, e.g. two-loop VP insertions)."""
    """Returns (C code for pt_<cname>, info dict)."""
    T = gr.terms(g); n = g.n; L = g.L
    if os.environ.get('GEN_BARE'): T = T[:1]
    nodes = {}
    def node_for(seg, F, op):
        kids = gr.children(seg, F)
        kind = 'F2' if op == 'A' else ('SE' if g.segment(*seg)['kind'] == 'SE' else 'V0')
        key = (seg, kids, kind)
        if key not in nodes: nodes[key] = Node(g, seg, kids, kind)
        return nodes[key]
    termspecs = []
    for sign, ops in T:
        F = tuple(ops.keys())
        termspecs.append((sign, ops, F))
        for seg, op in ops.items(): node_for(seg, F, op)
    code = []; tabs = {}; funcs = []
    split = os.environ.get('GEN_SPLIT', '1') != '0'       # one C function per node (lower compile memory)
    for key, nd in nodes.items():
        tabs[key] = nd.table(cache)
        if split:
            fn, call = node_function(nd, tabs[key], cname + '_'); funcs.append(fn); code.append(call)
        else:
            code.append(emit_node(nd, tabs[key]))
    em = Emitter(); maxtotal = 0
    # recursive value of (seg, subtree ops) -> (names of graded arrays: list per output kind)
    def value(seg, ops, F):
        op = ops[seg]; nd = node_for(seg, F, op)
        mk = (seg, tuple(sorted((s, o) for s, o in ops.items() if gr.contains(seg, s))))
        if mk in em.memo: return em.memo[mk]
        kids_vals = {c: value(c, ops, F) for c in nd.kids}
        if nd.kind == 'F2': outs = [(0, 'A')]
        elif nd.kind == 'V0': outs = [(0 if op == 'U' else 1, op)]
        else: outs = [(0, 'u'), (1, 'v')]
        res = []; mw = []
        for (k, lab) in outs:
            dst = em.new(); mwk = 0
            for pat in range(2**len(nd.selist)):
                cur = '%s_W[%d][%d]' % (nd.name, k, pat); cura = '%s_Wa[%d][%d]' % (nd.name, k, pat)
                cw = nd.maxw
                for t, c in enumerate(nd.selist):
                    uv = kids_vals[c]; sel = uv[1] if (pat >> t) & 1 else uv[0]
                    nxt = em.new(); em.conv_into(nxt, cur, sel[0], cura, sel[0]+'a'); cur, cura = nxt, nxt+'a'; cw += sel[1]
                for c in nd.vlist:
                    sel = kids_vals[c][0]
                    nxt = em.new(); em.conv_into(nxt, cur, sel[0], cura, sel[0]+'a'); cur, cura = nxt, nxt+'a'; cw += sel[1]
                em.lines.append("  for(int w=0;w<%d;w++){ %s[w]+=%s[w]; %sa[w]+=%s[w]; }" % (NW, dst, cur, dst, cura))
                mwk = max(mwk, cw)
            res.append((dst, mwk))
        em.memo[mk] = res
        return res
    em.lines.append("  REAL tot=0, tota=0;")
    p0 = n - 2*L
    for sign, ops, F in termspecs:
        G = (0, 2*L)
        (arr, mw), = value(G, ops, F)
        nonlocal_max = mw
        maxtotal = max(maxtotal, mw)
        Vs = "+".join(["%s_V" % node_for(s, F, o).name for s, o in ops.items()] + ["Vx"])
        em.lines.append("  { const REAL iV=1/(%s); REAL r=0, ra=0, ivp=1;" % Vs)
        if dpvar:
            em.lines.append("    for(int w=0;w<%d;w++){ const int q=%d+dp-w; if(q<1) continue; REAL f=1; ivp=1;" % (NW, p0))
            em.lines.append("      for(int j=1;j<q;j++) f*=j; for(int j=0;j<q;j++) ivp*=iV; r+=f*%s[w]*ivp; ra+=f*%sa[w]*ivp; }" % (arr, arr))
        else:
          # sum_w Gamma(p0-w) arr[w] iV^(p0-w)
          for w in range(p0-1, -1, -1):
            em.lines.append("    ivp*=iV;")
            if w < NW:
                em.lines.append("    r+=%d*%s[%d]*ivp; ra+=%d*%sa[%d]*ivp;" % (math.factorial(p0-w-1), arr, w, math.factorial(p0-w-1), arr, w))
        em.lines.append("    tot%s=r; tota+=ra; }" % ('+' if sign > 0 else '-'))
    if maxtotal >= p0:
        log("  WARNING: graded Wick count may reach %d >= n-2L=%d; check that those terms vanish" % (maxtotal, p0))
    body = funcs + ["static REAL FN(pt_%s)(const REAL *x, const REAL Vx, %sREAL *ab){" % (cname, "const int dp, " if dpvar else "")]
    body += ["  const REAL x%d=x[%d];" % (i, i) for i in range(n)]
    body += code + em.lines
    body.append("  *ab=tota; return tot;\n}")
    info = dict(nterms=len(T), nnodes=len(nodes), maxw=maxtotal,
                nstruct=sum(len(t) for t in tabs.values()))
    return "\n".join(body), info

def cw_arc(g):
    """The arc whose photon carries Cheng-Wu (parameter = 1): arc 0, which contains position 0, except in
    a self-mirror graph whose arc 0 is not its own mirror.  There the first arc that is its own mirror is
    used if there is one, so that the mirror permutes the remaining Sinc axes and the driver can fold the
    sum (sum.c, symmetry fold).  GEN_CWARC0=1 keeps arc 0 (the kernels before 2026-10-21)."""
    if os.environ.get('GEN_CWARC0', '0') == '1': return 0
    V = 2*g.L
    mir = lambda a: (V - a[1], V - a[0])
    arcs = [tuple(a) for a in g.arcs]
    if V - g.ext != g.ext or sorted(map(mir, arcs)) != sorted(arcs) or mir(arcs[0]) == arcs[0]: return 0
    inv = [k for k, a in enumerate(arcs) if mir(a) == a]
    return inv[0] if inv else 0

def classes(g):
    """Series classes of lines (identical momentum at q=0).  Returns (extclass, others) where
    extclass = (before-ext lines, after-ext lines) and others = list of lists of line indices.

    A row other than the external vertex's own may have lines on both sides of the external
    vertex (momenta p-l and p'-l, e.g. the first and last segment under an outer photon).  At q=0
    U, V and the contraction coefficients see only the sum of their parameters, and the point
    function is at most linear in the split between the two sides (checked on every such pair at
    two and three loops and on samples at four).  The split integral is then exact at the weighted
    mean, which puts every line of the row at sigma/K: the row is one series class.
    GEN_SPLITSIDES=1 keeps the two sides apart (the kernels before 2026-10-21)."""
    split = os.environ.get('GEN_SPLITSIDES', '0') == '1'
    exrow = tuple(1 if a <= g.ext-1 < b else 0 for a, b in g.arcs)
    rows = {}
    for i in range(g.ne):          # identical momentum at q=0: same loop combination
        row = tuple(1 if a <= i < b else 0 for a, b in g.arcs)
        key = (row, i < g.ext) if (split or row == exrow) else (row, None)
        rows.setdefault(key, []).append(i)
    before = rows.pop((exrow, True)); after = rows.pop((exrow, False))
    others = list(rows.values()) + [[2*g.L + k] for k in range(g.L)]
    fixed = 2*g.L + cw_arc(g)   # Cheng-Wu = 1 on this photon
    others = [c for c in others if c != [fixed]]
    return (before, after), others, fixed

def gen_wrapper(g, cname, inner=None, nextra=0, extranorm=1):
    """inner: name of a point function FN(inner)(s, v, &a) taking nextra extra Sinc variables v[]
    (appended after the graph's own); default: the graph's own pt with Vx = 0."""
    (bef, aft), others, fixed = classes(g)
    D = 1 + len(others) + nextra; k1, k2 = len(bef), len(aft)
    deg = g.ne + k1 + k2 - 2
    ngl = (deg + 2)//2
    xs, ws = gauss_legendre01(ngl)
    norm = Fraction(-1, 4)**g.L / 16 * extranorm
    glc = "static const double GLX_%s[%d]={%s}, GLW_%s[%d]={%s};" % (cname, ngl, ",".join(repr(float(x)) for x in xs), cname, ngl, ",".join(repr(float(w)) for w in ws))
    o = ["static REAL FN(f_%s)(const REAL *sf, REAL *ab){" % cname]
    o.append("  REAL s[%d]; s[%d]=1; const REAL sig=sf[0];" % (g.n, fixed))
    o.append("  REAL jac=1;")
    if nextra: o.append("  const REAL *v=sf+%d;" % (D - nextra))
    for d, c in enumerate(others):
        kk = len(c)
        for i in c: o.append("  s[%d]=sf[%d]/%d;" % (i, d+1, kk))
        if kk > 1: o.append("  jac*=%s;" % "*".join(["sf[%d]" % (d+1)]*(kk-1) + ["((REAL)1/%d)" % math.factorial(kk-1)]))
    o.append("  REAL acc=0, aa=0, a;")
    o.append("  for(int g=0;g<%d;g++){ const REAL u=(REAL)GLX_%s[g], sb=sig*u, sa=sig*(1-u);" % (ngl, cname))
    for i in bef: o.append("    s[%d]=sb/%d;" % (i, k1))
    for i in aft: o.append("    s[%d]=sa/%d;" % (i, k2))
    wexpr = "(REAL)GLW_%s[g]" % cname
    if k1 > 1: wexpr += "*" + "*".join(["sb"]*(k1-1)) + "/%d" % math.factorial(k1-1)
    if k2 > 1: wexpr += "*" + "*".join(["sa"]*(k2-1)) + "/%d" % math.factorial(k2-1)
    o.append("    const REAL wg=%s;" % wexpr)
    if inner: o.append("    acc+=wg*FN(%s)(s,v,&a); aa+=RABS(wg)*a; }" % inner)
    else: o.append("    acc+=wg*FN(pt_%s)(s,0,&a); aa+=RABS(wg)*a; }" % cname)
    o.append("  const REAL nrm=%s*sig*jac;" % cnum(norm))
    o.append("  *ab=RABS(nrm)*aa; return nrm*acc;\n}")
    return "\n".join(o), D, ngl, glc

def gauss_legendre01(n):
    import numpy as np
    x, w = np.polynomial.legendre.leggauss(n)
    return list((x+1)/2), list(w/2)


# ---------------- FORM output optimisation of the node polynomials ----------------
# With USE_FORM (set by gen.py when the FORM executable is found) each node's Dirac polynomial, and
# its |.| companion, is written as FORM 'Format O3' code (Horner scheme + common subexpressions),
# about 3x fewer operations than the factored tree + compact.py; the Schwinger polynomials still go
# through compact.py's factoring.  Every optimised expression is checked against the exact
# polynomial at a random rational point before it is used.
USE_FORM = False
FORM_LEVEL = os.environ.get('GEN_FORM_LEVEL', 'O3')

def emit_node(nd, tab, chk=None):
    if USE_FORM: return _emit_node_form(nd, tab)
    return _emit_node_py(nd, tab, chk)

def _fsym(f): return f.replace('_', 'x')

FORM_JOINT = os.environ.get('GEN_FORM_JOINT', '0') == '1'   # lblcommon.use_form turns it on for the lepton loops

def _crhs(rhs):
    """FORM Format C right-hand side -> our C (REAL rationals, no pow, x->_ in symbol names)."""
    rhs = re.sub(r'\s+', ' ', rhs).strip()
    rhs = re.sub(r'(\d+)\./(\d+)\.', r'((REAL)\1/\2)', rhs)
    rhs = re.sub(r'(?<![\w.])(\d+)\.(?![\d/])', r'\1', rhs)
    rhs = re.sub(r'pow\((\w+),(\d+)\)', lambda m: '(' + '*'.join([m.group(1)]*int(m.group(2))) + ')', rhs)
    rhs = re.sub(r'\b([a-zA-Z]+\d+)x(\d+)\b', r'\1_\2', rhs)
    assert 'pow' not in rhs and not re.search(r'\d\.', rhs), rhs
    return rhs

def _untag(stmts, ntag):
    """Straight-line code from FORM for E = sum_k tg<k> E_k  ->  code for the E_k separately.

    Every temporary is carried as its components (tag-free part, coefficient of each tag) in the
    ring with tg_i tg_j = 0.  Reducing modulo that ideal commutes with + and *, and E is linear in
    the tags, so the components of E are exactly the E_k, whatever order FORM's Horner scheme put the
    tags in.  Returns [(lhs, rhs)] in FORM's raw syntax with fresh temporaries Z<n>_ and outputs E<k>."""
    tok = re.compile(r'\s*(pow\(\w+,\d+\)|\d+\.?/\d+\.?|\d+\.?|\w+|[-+*])')
    comps = {}; out = []; cnt = [0]
    def atom(t):
        m = re.match(r'tg(\d+)$', t)
        if m: return {int(m.group(1)) + 1: '1'}
        m = re.match(r'pow\((Z\d+_),(\d+)\)$', t)
        if m:                                   # a power of a temporary: expand into a product
            return None, [m.group(1)]*int(m.group(2))
        if re.match(r'Z\d+_$', t): return comps[t]
        return {0: t}
    def mul(a, b):
        r = {}
        for ca, ea in a.items():
            for cb, eb in b.items():
                if ca and cb: continue
                c = ca or cb
                e = eb if ea == '1' else (ea if eb == '1' else '%s*%s' % (ea, eb))
                r.setdefault(c, []).append(e)
        return {c: (v[0] if len(v) == 1 else '(%s)' % '+'.join(v)) for c, v in r.items()}
    for lhs, rhs in stmts:
        toks = tok.findall(rhs); assert ''.join(toks).replace(' ', '') == rhs.replace(' ', ''), rhs
        terms = []; sign = 1; cur = None
        for t in toks + ['+']:
            if t in '+-':
                if cur is not None: terms.append((sign, cur))
                sign = 1 if t == '+' else -1; cur = None
            elif t == '*': continue
            else:
                a = atom(t)
                fs = [comps[z] for z in a[1]] if isinstance(a, tuple) else [a]
                for f in fs: cur = f if cur is None else mul(cur, f)
        acc = {}
        for s, d in terms:
            for c, e in d.items(): acc.setdefault(c, []).append(('-' if s < 0 else '+') + e)
        new = {}
        for c, es in acc.items():
            e = ''.join(es).lstrip('+')
            if re.match(r'(\w+|\d+\.?/\d+\.?|\d+\.?)$', e): new[c] = e; continue   # alias, no statement
            cnt[0] += 1; z = 'Z%d_' % cnt[0]; out.append((z, e)); new[c] = z
        comps[lhs] = new
    E = comps['E']
    assert set(E) <= set(range(1, ntag + 1)), 'tag-free part in a tagged expression'
    for k in range(ntag): out.append(('E%d' % k, E.get(k + 1, '0')))
    return out

def form_optimise(exprs, level=None, joint=None):
    """exprs: list of (name, [(Fraction, [factor names])]).  Returns [(lhs, rhs)] C statements
    (REAL rationals, no pow) computing every name, after checking them exactly.
    joint (GEN_FORM_JOINT=1, and the default of the lepton-loop generators, lblcommon.use_form):
    several expressions are optimised as one, E = sum_k tg_k E_k, so that FORM shares
    sub-expressions between them (it optimises each Local on its own), and split again exactly by
    _untag; the outputs are then named E0, E1, ... in order."""
    import subprocess, tempfile, random, formdirac
    level = level or FORM_LEVEL
    joint = (FORM_JOINT if joint is None else joint) and len(exprs) > 1
    syms = sorted({_fsym(f) for _, ts in exprs for _, fs in ts for f in fs})
    L = ['#-', 'Off statistics;', 'Format 250;']
    tags = ['tg%d' % k for k in range(len(exprs))] if joint else []
    if syms or tags: L.append('Symbols %s;' % ','.join(syms + tags))
    bodies = [''.join('%+d/%d*%s' % (c.numerator, c.denominator, '*'.join(_fsym(f) for f in fs) or '1') for c, fs in ts) for _, ts in exprs]
    if joint:
        L.append('Local E = %s;' % ('+'.join('tg%d*(%s)' % (k, b) for k, b in enumerate(bodies) if b) or '0'))
    else:
        for (nm, _), body in zip(exprs, bodies): L.append('Local %s = %s;' % (nm, body or '0'))
    L += ['.sort', 'Format C;', 'Format %s;' % level, 'Print;', '.end']
    with tempfile.TemporaryDirectory() as td:
        fn = os.path.join(td, 'o.frm'); open(fn, 'w').write('\n'.join(L) + '\n')
        r = subprocess.run([formdirac.FORM, '-q', fn], cwd=td, capture_output=True, text=True)
        if r.returncode: raise RuntimeError('FORM failed:\n' + r.stdout[-1500:])
    raw = [(lhs, re.sub(r'\s+', ' ', rhs).strip()) for lhs, rhs in re.findall(r'\s*(\w+)\s*=\s*(.*?);', r.stdout, flags=re.S)]
    if joint:
        raw = _untag(raw, len(exprs))
        ren = {'E%d' % k: nm for k, (nm, _) in enumerate(exprs)}
        raw = [(ren.get(l, l), r_) for l, r_ in raw]
    out = [(lhs, _crhs(rhs)) for lhs, rhs in raw]
    # exact check at a random rational point
    rnd = random.Random(12345)
    val = {f: Fraction(rnd.randint(-97, 97), rnd.randint(1, 61)) for _, ts in exprs for _, fs in ts for f in fs}
    env = dict(val)
    for lhs, rhs in out:
        e = re.sub(r'\(\(REAL\)(\d+)/(\d+)\)', r'_F(\1,\2)', rhs)
        env[lhs] = eval(e, {'_F': Fraction, '__builtins__': {}}, env)
    for nm, ts in exprs:
        exact = Fraction(0)
        for c, fs in ts:
            t = c
            for f in fs: t *= val[f]
            exact += t
        assert env.get(nm, Fraction(0)) == exact, ('FORM output check failed', nm)
    return out

def _emit_node_form(nd, tab):
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tools'))
    import compact
    base = _emit_node_py(nd, tab, None).split('\n')
    groups = {}
    for (combo, pr), val in tab.items():
        w = len(pr)
        pat = 0
        for t, c in enumerate(nd.selist):
            if combo[len(nd.egaps) + t] != '1': pat |= 1 << t
        pw = Fraction(-1, 2)**w
        vals = [val*pw] if nd.kind == 'F2' else [Fraction(x)*pw for x in val]
        for m, fac in _monomials(nd, combo, pr):
            for k, x in enumerate(vals):
                if x != 0: groups.setdefault((k, pat, w), []).append((x*m, fac))
    nm = nd.name; keys = sorted(groups)
    keep = [l for l in base if not re.match(r'\s*const REAL A\w+=RABS', l) and not re.match(r'\s*%s_W\[.*\+=.*; %s_Wa' % (re.escape(nd.name), re.escape(nd.name)), l)]
    ins = max(i for i, l in enumerate(keep) if re.match(r'\s*const REAL (V|c\w+|w\d+_\d+)=', l)) + 1
    def merged(ts):   # identical monomials summed first, as in the factored tree's |.| bound
        d = {}
        for c, fs in ts: d[tuple(sorted(fs))] = d.get(tuple(sorted(fs)), 0) + c
        return [(c, list(fs)) for fs, c in d.items() if c != 0]
    blk = []
    if keys:
        mg = [merged(groups[k]) for k in keys]
        val = form_optimise([('E%d' % n, mg[n]) for n in range(len(keys))])
        ab = form_optimise([('E%d' % n, [(abs(c), ['a' + f for f in fs]) for c, fs in mg[n]]) for n in range(len(keys))])
        def decl(st):
            ts = sorted({l for l, _ in st if l.startswith('Z')}, key=lambda z: int(z[1:-1]))
            return 'REAL %s; ' % ', '.join(ts + ['E%d' % n for n in range(len(keys))])
        blk.append('    { ' + decl(val))
        blk += ['      %s=%s;' % lr for lr in val]
        blk += ['      %s_W[%d][%d][%d]+=E%d;' % ((nm,) + k + (n,)) for n, k in enumerate(keys)]
        blk.append('    }')
        facs = sorted({f for k in keys for _, fs in groups[k] for f in fs})
        blk.append('    if(SQ_AB){ ' + decl(ab))
        if facs: blk.append('      const REAL %s;' % ', '.join('a%s=RABS(%s)' % (f, f) for f in facs))
        blk += ['      %s=%s;' % lr for lr in ab]
        blk += ['      %s_Wa[%d][%d][%d]+=E%d;' % ((nm,) + k + (n,)) for n, k in enumerate(keys)]
        blk.append('    }')
    keep[ins:ins] = blk
    return '\n'.join(compact.compact_x(keep, '  ', nd.name + '_'))   # tag: unique names when nodes are inline
