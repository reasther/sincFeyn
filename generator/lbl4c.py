"""
Four-loop light-by-light graphs with an extra photon r on the open line (Kinoshita group IV(c),
48 diagrams, Laporta set 21).

Open line: vertices 0..4; r joins positions (i, j), the three light-by-light photons a1..a3 sit at
the other positions in order.  (i, j) = (0,1) or (3,4) is an external-leg correction (excluded).
Mirror pairs (weight 2): (0,2)~(2,4), (0,3)~(1,4), (1,2)~(2,3); self-mirror (0,4), (1,3).
Loop: vertices L0 (external photon mu), L1..L3 carrying a1..a3 in one of 6 orderings.
Parameters: 0-3 e0..e3 (e_q joins positions q, q+1), 4-7 l0..l3, 8-10 a1..a3, 11 r.
Loop momenta (k1, k2, l, r); a3 = p - p' - k1 - k2.

Subtractions (Volkov, G' in I[G]; in every term the loop is evaluated minus itself at zero photon
momenta, per ordering, and the 6 orderings are summed pointwise, as in lbl.py):
  (1,2):          A(G) - A(G/S) U(S)               S = r self-energy on e1
  (0,2), (1,3):   A(G) - A(G/V) U(V)               V = r vertex around the photon inside
  (0,3):          A(G)
  (0,4):          A(G) - A(G/H) U(H) - L(G/H) A(H) + U(G/H) A(H)
                  H = the three-loop light-by-light vertex (positions 1..3 and the loop)

    python3 lbl4c.py OUTDIR [JOBS]    -> l4_c02.c l4_c03.c l4_c04.c l4_c12.c l4_c13.c
"""
import itertools, math, os, sys, pickle, hashlib, multiprocessing as mp, time
from fractions import Fraction
import build, dirac

PERMS = list(itertools.permutations([1, 2, 3]))      # perm[a-1] = loop vertex (1..3) of photon a
LP0, AP0, RP = 4, 8, 11
NPAR = 12
PLACEMENTS = [((0, 2), 2), ((0, 3), 2), ((0, 4), 1), ((1, 2), 2), ((1, 3), 1)]

class N:
    pass

def node(name, kind, chain, loop, A, B, gidx, nlines, selist=(), slotline=None):
    n = N(); n.name, n.kind, n.chain, n.loop, n.A, n.B, n.gidx = name, kind, chain, loop, A, B, gidx
    n.m2 = [1 if g < AP0 else 0 for g in gidx]
    n.selist = list(selist); n.egaps = list(range(nlines))
    n.slots = [None]*(nlines + len(selist))
    n.slotline = slotline if slotline is not None else list(range(nlines))
    n.L = len(A[0]); n.seg = name; n.kids = ()
    return n

def tkey(n):
    return hashlib.md5((repr((n.chain, len(n.slots), n.kind)) + (repr(n.loop) if n.loop else '')).encode()).hexdigest()

def table(n, cache):
    k = tkey(n)
    if k not in cache: cache[k] = dirac.table(n.chain, len(n.slots), n.kind, loop=n.loop)
    return cache[k]

KA = {1: [1, 0, 0, 0], 2: [0, 1, 0, 0], 3: [-1, -1, 0, 0], 'r': [0, 0, 0, 1]}
KB = {1: [0, 0], 2: [0, 0], 3: [1, -1], 'r': [0, 0]}
add = lambda a, b: [x + y for x, y in zip(a, b)]
sub = lambda a, b: [x - y for x, y in zip(a, b)]

class Placement:
    def __init__(self, i, j):
        self.i, self.j = i, j
        att = [p for p in range(5) if p not in (i, j)]
        self.ph = {i: 'r', j: 'r'}
        for t, p in enumerate(att): self.ph[p] = t + 1
    def open_row(self, q, v0=False):
        """Momentum of e_q: p minus what was emitted at positions <= q (r out at i, back at j)."""
        a, b = [0, 0, 0, 0], [1, 0]
        for p in range(q + 1):
            x = self.ph[p]
            if x == 'r': a = sub(a, KA['r']) if p == self.i else add(a, KA['r'])
            else: a = sub(a, KA[x]); b = sub(b, KB[x])
        return a, ([1, 0] if v0 else b)
    def seq(self, drop_pos=(), drop_lines=(), se_line=None, contract=None):
        """Chain (matrix order) over the open line; drop_pos: vertices removed; drop_lines: lines
        removed; se_line: line replaced by the SE slot 'S'; contract: (lo, hi) positions merged into
        the external vertex 'mu'.  Items: ('v', name) or ('s', ('e', q)) / ('s', 'S')."""
        out = []
        for p in range(4, -1, -1):
            if contract and contract[0] <= p <= contract[1]:
                if p == contract[0]: out.append(('v', 'mu'))
            elif p not in drop_pos:
                x = self.ph[p]; out.append(('v', 'g' if x == 'r' else 'a%d' % x))
            if p > 0:
                q = p - 1
                if contract and contract[0] <= q < contract[1]: continue
                if q == se_line: out.append(('s', 'S'))
                elif q not in drop_lines: out.append(('s', ('e', q)))
        return out

def loop_rows(perm, zero, cols):
    inv = {perm[a]: a + 1 for a in range(3)}
    A = []; B = []; la, lb = [0, 0, 1, 0], [0, 0]
    for q in range(4):
        A.append([la[c] for c in cols]); B.append(list(lb))
        if q < 3 and not zero: a = inv[q+1]; la = add(la, KA[a]); lb = add(lb, KB[a])
    return A, B, inv

def lbl_node(P, name, kind, perm, zero, items, lines, photons, v0=False):
    """Node with the loop: items = chain over the open line (from P.seq); lines = kept open lines
    (ascending); photons = photon labels (1,2,3,'r') of the node."""
    cols = [0, 1, 2] + ([3] if 'r' in photons else [])
    se = any(t == ('s', 'S') for t in items)
    nl = len(lines) + 4
    slot = {q: t for t, q in enumerate(lines)}
    chain = [('s', slot[x[1][1]]) if x[0] == 's' and x[1] != 'S' else (('s', nl) if x[0] == 's' else x) for x in items]
    A = []; B = []
    for q in lines:
        a, b = P.open_row(q, v0); A.append([a[c] for c in cols]); B.append(b)
    lA, lB, inv = loop_rows(perm, zero, cols)
    if v0: lB = [[0, 0]]*4
    A += lA; B += lB
    for x in photons:
        A.append([KA[x][c] for c in cols]); B.append([0, 0] if v0 else list(KB[x]))
    toks = []
    for q in range(3, -1, -1):
        toks.append(('s', len(lines) + q)); toks.append(('v', 'a%d' % inv[q] if q else 'mu'))
    gidx = list(lines) + [LP0 + q for q in range(4)] + [RP if x == 'r' else AP0 + x - 1 for x in photons]
    slotline = None
    if se:
        slotline = list(range(nl)) + [slot[P.i - 1]]
    return node(name, kind, chain, toks, A, B, gidx, nl, selist=(['S'] if se else ()), slotline=slotline)

def v1_node(name, kind, l1, l2):
    """One-loop vertex (lines l1 before, l2 after, photon r) as F2 or V0 node."""
    B = [[1, 0], [0, 1], [0, 0]] if kind == 'F2' else [[1, 0], [1, 0], [0, 0]]
    return node(name, kind, [('v', 'g'), ('s', 1), ('v', 'mu'), ('s', 0), ('v', 'g')], None,
                [[-1], [-1], [1]], B, [l1, l2, RP], 2)

def terms_for(i, j, perm, zero):
    """[(sign, [(node, k)], selist-kid or None)] for one ordering; nodes with the loop carry zero."""
    P = Placement(i, j); tag = '%d%d_%s%s' % (i, j, ''.join(map(str, perm)), '_z' if zero else '')
    T = []
    full = lbl_node(P, 'n_c%s_G' % tag, 'F2', perm, zero, P.seq(), [0, 1, 2, 3], [1, 2, 3, 'r'])
    T.append((1, [(full, 0)], None))
    if j == i + 1:
        keep = [q for q in range(4) if q != i]
        red = lbl_node(P, 'n_c%s_GS' % tag, 'F2', perm, zero, P.seq(drop_pos=(i, j), se_line=i), keep, [1, 2, 3])
        S = node('n_c%d%d_S' % (i, j), 'SE', [('v', 'g'), ('s', 0), ('v', 'g')], None, [[-1], [1]], [[1, 0], [0, 0]], [i, RP], 1)
        T.append((-1, [(red, 0)], S))
    elif j == i + 2:
        keep = [q for q in range(4) if q not in (i, i + 1)]
        red = lbl_node(P, 'n_c%s_GV' % tag, 'F2', perm, zero, P.seq(drop_pos=(i, j), drop_lines=(i, i + 1)), keep, [1, 2, 3])
        V = v1_node('n_c%d%d_V' % (i, j), 'V0', i, i + 1)
        T.append((-1, [(red, 0), (V, 0)], None))
    elif (i, j) == (0, 4):
        HA = lbl_node(P, 'n_c%s_HA' % tag, 'F2', perm, zero, P.seq(drop_pos=(0, 4), drop_lines=(0, 3)), [1, 2], [1, 2, 3])
        HU = lbl_node(P, 'n_c%s_HU' % tag, 'V0', perm, zero, P.seq(drop_pos=(0, 4), drop_lines=(0, 3)), [1, 2], [1, 2, 3], v0=True)
        GF = v1_node('n_c04_GHF', 'F2', 0, 3); GV = v1_node('n_c04_GHV', 'V0', 0, 3)
        T += [(-1, [(GF, 0), (HU, 0)], None), (-1, [(GV, 1), (HA, 0)], None), (1, [(GV, 0), (HA, 0)], None)]
    return T

def all_nodes(i, j):
    out = {}
    for perm in PERMS:
        for zero in (False, True):
            for sg, prod, se in terms_for(i, j, perm, zero):
                for nd, k in prod: out[nd.name] = nd
                if se is not None: out[se.name] = se
    return out

def _tab(nd):
    t = time.time(); tb = dirac.table(nd.chain, len(nd.slots), nd.kind, loop=nd.loop)
    return tkey(nd), tb, time.time() - t

def pt_code(cache, cname, i, j):
    nodes = all_nodes(i, j)
    import lblcommon
    funcs = []; code = []
    for nm, nd in nodes.items():
        fn, call = lblcommon.node_function(nd, table(nd, cache), NPAR); funcs.append(fn); code.append(call)
    em = build.Emitter(); p0 = NPAR - 2*4
    em.lines.append("  REAL tot=0, tota=0;")
    for perm in PERMS:
        for zero in (False, True):
            for sg, prod, se in terms_for(i, j, perm, zero):
                s = sg*(-1 if zero else 1)
                # graded product of the node outputs (with the SE patterns if any)
                dst = em.new()
                npat = 2 if se is not None else 1
                for pat in range(npat):
                    nd0, k0 = prod[0]
                    cur, cura = '%s_W[%d][%d]' % (nd0.name, k0, pat), '%s_Wa[%d][%d]' % (nd0.name, k0, pat)
                    for nd, k in prod[1:] + ([(se, pat)] if se is not None else []):
                        nxt = em.new(); em.conv_into(nxt, cur, '%s_W[%d][0]' % (nd.name, k), cura, '%s_Wa[%d][0]' % (nd.name, k)); cur, cura = nxt, nxt + 'a'
                    em.lines.append("  for(int w=0;w<%d;w++){ %s[w]+=%s[w]; %sa[w]+=%s[w]; }" % (build.NW, dst, cur, dst, cura))
                Vs = "+".join(["%s_V" % nd.name for nd, k in prod] + (["%s_V" % se.name] if se is not None else []) + ["Vx"])
                em.lines.append("  { const REAL iV=1/(%s); REAL r=0, ra=0, ivp=1;" % Vs)
                for w in range(p0-1, -1, -1):
                    em.lines.append("    ivp*=iV;")
                    if w < build.NW:
                        f = math.factorial(p0-w-1)
                        em.lines.append("    r+=%d*%s[%d]*ivp; ra+=%d*%sa[%d]*ivp;" % (f, dst, w, f, dst, w))
                em.lines.append("    tot%s=r; tota+=ra; }" % ('+' if s > 0 else '-'))
    pt = ["static REAL FN(pt_%s)(const REAL *x, const REAL Vx, REAL *ab){" % cname]
    pt += ["  const REAL x%d=x[%d];" % (q, q) for q in range(NPAR)]
    pt += code + em.lines + ["  *ab=tota; return tot;\n}"]
    return "\n".join(funcs) + "\n" + "\n".join(pt)

def write(out, cache, i, j, w0):
    cname = 'l4_c%d%d' % (i, j)
    pt = pt_code(cache, cname, i, j)
    # ext class {l0, l3} (GL); Cheng-Wu: photon a1 = 1, or a2 in the self-mirror placements (0,4) and
    # (1,3), so that the open-line mirror (a1 <-> a3) keeps it fixed and the driver can fold the sum
    # over it as well as over the loop reversal (GEN_LBLCW0=1: a1 always, the kernels before
    # 2026-10-26).  Open lines with the same momentum at q=0 (e0, e3 for (0,4); e0, e2 for (1,2))
    # enter the point function only through the sum of their parameters (checked numerically), so
    # they are one Sinc variable sigma with each line at sigma/2 (GEN_SPLITSIDES=1: separate).
    fixed = AP0 + 1 if (i, j) in ((0, 4), (1, 3)) and os.environ.get('GEN_LBLCW0', '0') != '1' else AP0
    P = Placement(i, j); rows = {}
    for q in range(4):
        key = q if os.environ.get('GEN_SPLITSIDES', '0') == '1' else repr(P.open_row(q, v0=True))
        rows.setdefault(key, []).append(q)
    classes = list(rows.values()) + [[q] for q in range(4, NPAR) if q not in (LP0, LP0 + 3, fixed)]
    ngl = (8 + 2)//2
    xs, ws = build.gauss_legendre01(ngl)
    glc = "static const double GLX_%s[%d]={%s}, GLW_%s[%d]={%s};" % (cname, ngl, ",".join(repr(float(x)) for x in xs), cname, ngl, ",".join(repr(float(w)) for w in ws))
    norm = -Fraction(-1, 4)**4/16
    D = 1 + len(classes)
    wr = ["static REAL FN(f_%s)(const REAL *sf, REAL *ab){" % cname,
          "  REAL s[%d]; const REAL sig=sf[0]; s[%d]=1;" % (NPAR, fixed)]
    for d, c in enumerate(classes):
        wr += ["  s[%d]=sf[%d]%s;" % (q, d + 1, '/%d' % len(c) if len(c) > 1 else '') for q in c]
    jac = ''.join('*sf[%d]' % (d + 1) for d, c in enumerate(classes) if len(c) > 1)
    assert all(len(c) <= 2 for c in classes)
    wr += ["  REAL acc=0, aa=0, a;",
           "  for(int g=0;g<%d;g++){ const REAL u=(REAL)GLX_%s[g]; s[%d]=sig*u; s[%d]=sig*(1-u);" % (ngl, cname, LP0, LP0 + 3),
           "    const REAL wg=(REAL)GLW_%s[g]; acc+=wg*FN(pt_%s)(s,0,&a); aa+=wg*a; }" % (cname, cname),
           "  const REAL nrm=%s*sig%s;" % (build.cnum(norm), jac),
           "  *ab=RABS(nrm)*aa; return nrm*acc;\n}"]
    wr = "\n".join(wr)
    title = 'LbL + photon on the line %d-%d (6 orderings)' % (i + 1, j + 1)
    src = ["/* Generated by lbl4c.py: %s */" % title, '#include "common.h"', glc,
           '#define REAL double', '#define FN(x) x##_d', pt, wr, '#undef REAL', '#undef FN',
           '#define REAL QREAL', '#define FN(x) x##_q', pt, wr, '#undef REAL', '#undef FN',
           'const graph_t G_%s={"%s","%s",4,%d,%d,%d,%d,f_%s_d,f_%s_q,pt_%s_d,pt_%s_q};' %
           (cname, cname, title, D, w0, NPAR, fixed, cname, cname, cname, cname)]
    open(os.path.join(out, cname + '.c'), 'w').write("\n".join(src) + "\n")
    return cname

_W = None
def _write(a): (i, j), w = a; return write(_W[0], _W[1], i, j, w)

def main(out, jobs, only=None):
    import lblcommon; lblcommon.use_form()
    cf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache_lbl.pkl')
    cache = pickle.load(open(cf, 'rb')) if os.path.exists(cf) else {}
    pl = [(ij, w) for ij, w in PLACEMENTS if only is None or ij in only]
    need = {}
    for (i, j), w in pl:
        for nd in all_nodes(i, j).values():
            if tkey(nd) not in cache: need[tkey(nd)] = nd
    print('tables needed:', len(need), flush=True)
    if need:
        with mp.get_context('fork').Pool(jobs) as pool:
            for k, tb, dt in pool.imap_unordered(_tab, sorted(need.values(), key=lambda n: -len(n.slots))):
                cache[k] = tb; print('  table %s: %d structures, %.0fs' % (k[:8], len(tb), dt), flush=True)
                c2 = pickle.load(open(cf, 'rb')) if os.path.exists(cf) else {}
                c2.update(cache); pickle.dump(c2, open(cf + '.tmp', 'wb')); os.replace(cf + '.tmp', cf)
    global _W; _W = (out, cache)                          # the forked workers see it without pickling
    with mp.get_context('fork').Pool(jobs) as pool: names = pool.map(_write, pl)
    import lblvp; lblvp.write_table(out)
    print(names)
    return names

if __name__ == '__main__':
    only = None
    if len(sys.argv) > 3: only = [tuple(int(c) for c in s) for s in sys.argv[3].split(',')]
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 4, only)
