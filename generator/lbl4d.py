"""
Four-loop graphs with a light-by-light loop joined to the open line by four photons and the
external photon on the open line (Kinoshita group IV(d), 18 diagrams, Laporta sets 24+25).

Open line: vertices at positions 0..4, the external photon at position m (1, 2 or 3; m = 0 or 4
would make the graph one-particle reducible); the other four positions carry photons a1..a4 in
order.  Electron line e_i joins positions i and i+1 (momentum p before m, p' after, minus the
photon momenta emitted so far).  Loop: vertices L0..L3 with a1 at L0 and a permutation of a2..a4 at
L1..L3 (6 orderings, both orientations); line l_i runs L_i -> L_{i+1}.
Parameters: 0-3 e0..e3, 4-7 l0..l3, 8-11 a1..a4.  Loop momenta k1, k2, k3, l; k4 = -k1-k2-k3.

The loop's own divergence (a constant 4-photon tensor) cancels in the sum over orderings, so, as
in lbl.py, each ordering is evaluated minus the same loop at zero photon momenta, and the 6
orderings are summed pointwise.  For m=1 (and its mirror m=3) the loop with a2..a4 and the open
lines e2, e3 form a three-loop light-by-light vertex H around photon a1 (not containing mu),
log divergent: that term is A(G) - A(G/H) U(H), G/H being the one-loop vertex around mu.
(Found by the scale test; m=2 has no such subgraph.)

The two mirror images m=1 and m=3 are equal (weight 2), m=2 is its own mirror.

    python3 lbl4d.py OUTDIR [JOBS]    -> l4_d1.c, l4_d2.c   (tables cached in cache_lbl.pkl)
"""
import itertools, math, os, sys, pickle, hashlib, multiprocessing as mp, time
from fractions import Fraction
import build, dirac, lbl4c

NE, NL, NPH = 4, 4, 4

class DNode:
    def __init__(self, m, perm, zero=False):
        self.m, self.perm, self.zero = m, perm, zero            # perm[j] = photon (1..4) at loop vertex j
        att = [p for p in range(5) if p != m]                   # positions of a1..a4
        # photon momenta as loop-momentum rows (k1,k2,k3,l) and external rows (p,p')
        kA = {1: [1, 0, 0, 0], 2: [0, 1, 0, 0], 3: [0, 0, 1, 0], 4: [-1, -1, -1, 0]}
        A = []; B = []
        for i in range(NE):
            a = [0, 0, 0, 0]
            for t, pos in enumerate(att):
                if pos <= i: a = [x - y for x, y in zip(a, kA[t+1])]
            A.append(a); B.append([1, 0] if i < m else [0, 1])
        la = [0, 0, 0, 1]
        for j in range(NL):
            A.append(list(la)); B.append([0, 0])
            if j < NL-1 and not zero: la = [x + y for x, y in zip(la, kA[perm[j+1]])]
        for t in range(1, 5): A.append(kA[t]); B.append([0, 0])
        self.A, self.B = A, B
        self.m2 = [1]*(NE + NL) + [0]*NPH
        n = NE + NL + NPH
        # m=2 is its own mirror (a1 <-> a4, a2 <-> a3, e0 <-> e3), which maps loop line j of one ordering
        # onto a line of another ordering that sits at a rotated position.  Each ordering is its own
        # integral, so its loop lines may take the parameters in any order: the orderings whose loop runs
        # a1 -> perm[1] with perm[1] > perm[3] (the reversals of the other three) put line j on parameter
        # l_{j+1 mod 4}.  The mirror and the loop reversal are then both axis permutations of the summed
        # integrand (l1 <-> l3; and l0 <-> l2, l1 <-> l3 with the photons and open lines), and the driver
        # folds the sum over both.  GEN_LBLCW0=1: the same labels in every ordering (kernels before 2026-10-26).
        rot = m == 2 and perm[1] > perm[3] and os.environ.get('GEN_LBLCW0', '0') != '1'
        self.gidx = list(range(NE)) + [NE + ((j + 1) % NL if rot else j) for j in range(NL)] + list(range(NE + NL, n)); self.slotline = list(range(NE + NL)); self.selist = []; self.egaps = list(range(NE + NL))
        self.kind = 'F2'; self.L = 4; self.seg = ('lbl4d', m) + tuple(perm); self.kids = ()
        self.slots = [('e', i) for i in range(NE + NL)]
        name = lambda pos: 'mu' if pos == m else 'a%d' % (att.index(pos) + 1)
        ch = []
        for pos in range(4, -1, -1):
            ch.append(('v', name(pos)))
            if pos > 0: ch.append(('s', pos - 1))
        self.chain = ch
        lp = []
        for j in range(NL-1, -1, -1):
            lp.append(('s', NE + j)); lp.append(('v', 'a%d' % perm[j]))
        self.loop = lp
        self.name = 'n_d%d_%s' % (m, ''.join(str(x) for x in perm)) + ('_z' if zero else '')
    def key(self):
        return hashlib.md5((repr((self.chain, len(self.slots), self.kind)) + repr(self.loop)).encode()).hexdigest()
    def table(self, cache):
        return cache[self.key()]

PERMS = [(1,) + p for p in itertools.permutations([2, 3, 4])]

def _tab(nd):
    t = time.time(); tb = dirac.table(nd.chain, NE + NL, 'F2', loop=nd.loop)
    return nd.key(), tb, time.time() - t

def h_node(m, perm, zero):
    """U(H) for m=1: H = the loop + a2..a4 + open lines e2, e3, a three-loop light-by-light vertex
    around photon a1 (which is H's external photon).  V0 node at H's own zero-momentum kinematics.
    Photons are relabelled a2..a4 -> a1..a3 and a1 -> mu so the Dirac table is lbl4c's HU table."""
    assert m == 1
    d = DNode(m, perm, zero)
    rows = [2, 3] + [NE + j for j in range(NL)] + [NE + NL + t for t in (1, 2, 3)]
    A = [[d.A[r][c] for c in (1, 2, 3)] for r in rows]
    B = [[1, 0], [1, 0]] + [[0, 0]]*(len(rows) - 2)
    ren = lambda a: 'mu' if a == 1 else 'a%d' % (a - 1)
    chain = [('v', 'a3'), ('s', 1), ('v', 'a2'), ('s', 0), ('v', 'a1')]
    toks = []
    for j in range(NL-1, -1, -1):
        toks.append(('s', 2 + j)); toks.append(('v', ren(perm[j])))
    nd = lbl4c.node('n_d%d_%s_HU%s' % (m, ''.join(map(str, perm)), '_z' if zero else ''), 'V0', chain, toks, A, B, rows, 6)
    nd.zero = zero
    return nd

def gh_node():
    """A(G/H) for m=1: one-loop vertex around mu, photon a1 from position 0 to the contracted H."""
    return lbl4c.node('n_d1_GH', 'F2', [('v', 'g'), ('s', 1), ('v', 'mu'), ('s', 0), ('v', 'g')], None,
                      [[-1], [-1], [1]], [[1, 0], [0, 1], [0, 0]], [0, 1, NE + NL], 2)

def terms(m):
    """[(sign, zero, [(node, k)])]: A(G) per ordering; for m=1 also -A(G/H) U(H)."""
    T = []
    for zero in (False, True):
        for p in PERMS:
            T.append((1, zero, [(DNode(m, p, zero), 0)]))
            if m == 1: T.append((-1, zero, [(gh_node(), 0), (h_node(m, p, zero), 0)]))
    return T

def tab_of(nd, cache):
    return cache[nd.key()] if isinstance(nd, DNode) else cache[lbl4c.tkey(nd)]

def pt_code(cache, cname, m):
    n = NE + NL + NPH; p0 = n - 2*4
    import lblcommon
    nodes = {}
    for sg, zero, prod in terms(m):
        for nd, k in prod: nodes[nd.name] = nd
    funcs = []; code = []
    for nd in nodes.values():
        fn, call = lblcommon.node_function(nd, tab_of(nd, cache), n); funcs.append(fn); code.append(call)
    em = build.Emitter(); em.lines.append("  REAL tot=0, tota=0;")
    for sg, zero, prod in terms(m):
        s = sg*(-1 if zero else 1)
        nd0, k0 = prod[0]
        cur, cura = '%s_W[%d][0]' % (nd0.name, k0), '%s_Wa[%d][0]' % (nd0.name, k0)
        for nd, k in prod[1:]:
            nxt = em.new(); em.conv_into(nxt, cur, '%s_W[%d][0]' % (nd.name, k), cura, '%s_Wa[%d][0]' % (nd.name, k)); cur, cura = nxt, nxt + 'a'
        Vs = "+".join(["%s_V" % nd.name for nd, k in prod] + ["Vx"])
        em.lines.append("  { const REAL iV=1/(%s); REAL r=0, ra=0, ivp=1;" % Vs)
        for w in range(p0-1, -1, -1):
            em.lines.append("    ivp*=iV;")
            if w < build.NW:
                f = math.factorial(p0-w-1)
                em.lines.append("    r+=%d*%s[%d]*ivp; ra+=%d*%s[%d]*ivp;" % (f, cur, w, f, cura, w))
        em.lines.append("    tot%s=r; tota+=ra; }" % ('+' if s > 0 else '-'))
    pt = ["static REAL FN(pt_%s)(const REAL *x, const REAL Vx, REAL *ab){" % cname]
    pt += code + em.lines + ["  *ab=tota; return tot;\n}"]
    return "\n".join(funcs) + "\n" + "\n".join(pt)

def write(out, cache, m, w0):
    cname = 'l4_d%d' % m
    pt = pt_code(cache, cname, m)
    # ext class {e_{m-1}, e_m} (GL), photon a1 fixed = 1, all other parameters free.  m=2 is its own
    # mirror (a1 <-> a4, a2 <-> a3), which no photon survives: there Cheng-Wu fixes the ext-class sum
    # sigma = 1 instead and a1 is free, so the driver can fold the sum over the mirror as well as over
    # the loop reversal (GEN_LBLCW0=1: a1 = 1, the kernels before 2026-10-26).
    ea, eb = m - 1, m
    cwsig = m == 2 and os.environ.get('GEN_LBLCW0', '0') != '1'
    free = [i for i in range(NE + NL + NPH) if i not in (ea, eb) + (() if cwsig else (NE + NL,))]
    deg = (NE + NL) + 2 - 2; ngl = (deg + 2)//2
    xs, ws = build.gauss_legendre01(ngl)
    glc = "static const double GLX_%s[%d]={%s}, GLW_%s[%d]={%s};" % (cname, ngl, ",".join(repr(float(x)) for x in xs), cname, ngl, ",".join(repr(float(w)) for w in ws))
    norm = -Fraction(-1, 4)**4/16          # fermion loop sign
    D = len(free) + (0 if cwsig else 1)
    wr = ["static REAL FN(f_%s)(const REAL *sf, REAL *ab){" % cname,
          "  REAL s[%d]; const REAL sig=1;" % (NE + NL + NPH) if cwsig else
          "  REAL s[%d]; const REAL sig=sf[0]; s[%d]=1;" % (NE + NL + NPH, NE + NL)]
    wr += ["  s[%d]=sf[%d];" % (i, d + (0 if cwsig else 1)) for d, i in enumerate(free)]
    wr += ["  REAL acc=0, aa=0, a;",
           "  for(int g=0;g<%d;g++){ const REAL u=(REAL)GLX_%s[g]; s[%d]=sig*u; s[%d]=sig*(1-u);" % (ngl, cname, ea, eb),
           "    const REAL wg=(REAL)GLW_%s[g]; acc+=wg*FN(pt_%s)(s,0,&a); aa+=wg*a; }" % (cname, cname),
           "  const REAL nrm=%s*sig;" % build.cnum(norm),
           "  *ab=RABS(nrm)*aa; return nrm*acc;\n}"]
    wr = "\n".join(wr)
    title = 'LbL loop with 4 photons to the line, mu at %d (6 orderings)' % (m + 1)
    src = ["/* Generated by lbl4d.py: %s */" % title, '#include "common.h"', glc,
           '#define REAL double', '#define FN(x) x##_d', pt, wr, '#undef REAL', '#undef FN',
           '#define REAL QREAL', '#define FN(x) x##_q', pt, wr, '#undef REAL', '#undef FN',
           'const graph_t G_%s={"%s","%s",4,%d,%d,%d,%d,f_%s_d,f_%s_q,pt_%s_d,pt_%s_q};' %
           (cname, cname, title, D, w0, NE + NL + NPH, NE + NL, cname, cname, cname, cname)]
    open(os.path.join(out, cname + '.c'), 'w').write("\n".join(src) + "\n")
    return cname

_W = None
def _write(a): return write(_W[0], _W[1], *a)

def main(out, jobs):
    import lblcommon; lblcommon.use_form()
    cf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache_lbl.pkl')
    cache = pickle.load(open(cf, 'rb')) if os.path.exists(cf) else {}
    nodes = [DNode(m, p) for m in (1, 2) for p in PERMS]
    need = {nd.key(): nd for nd in nodes if nd.key() not in cache}
    print('tables needed:', len(need), flush=True)
    if need:
        with mp.Pool(jobs) as pool:
            for k, tb, dt in pool.imap_unordered(_tab, list(need.values())):
                cache[k] = tb; print('  table %s: %d structures, %.0fs' % (k[:8], len(tb), dt), flush=True)
                c2 = pickle.load(open(cf, 'rb')) if os.path.exists(cf) else {}
                c2.update(cache); pickle.dump(c2, open(cf + '.tmp', 'wb')); os.replace(cf + '.tmp', cf)
    global _W; _W = (out, cache)                          # the forked workers see it without pickling
    with mp.get_context('fork').Pool(jobs) as pool: names = pool.map(_write, [(1, 2), (2, 1)])
    import lblvp; lblvp.write_table(out)
    print(names)
    return names

if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 4)
