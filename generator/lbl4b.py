"""
Four-loop light-by-light graphs with an extra photon r whose both ends sit on the lepton loop
(Kinoshita group IV(b), 60 diagrams, Laporta set 23, 0.822284485811034).

Open line: vertices 0..2 carrying a1, a2, a3 in order; lines e0, e1.
Loop: six vertices, position 0 = external photon mu, positions 1..5 = a1..a3 and both ends of r
in every arrangement (5!/2 = 60 sequences, both orientations); loop line l_q runs q -> q+1 (mod 6).
Parameters: 0-1 e0 e1, 2-7 l0..l5, 8-10 a1..a3, 11 r.  Loop momenta (k1, k2, l, r); a3 = p-p'-k1-k2.

Every sequence is evaluated minus the same graph at zero external photon momenta on the loop
(r kept), as in lbl.py; this removes the overall divergence of the two-loop 4-photon subgraph per
sequence, and it cancels in the sum over all 60 sequences, which is taken pointwise (the 4-photon
tensor is gauge invariant only in that sum, so all 60 are one integrand).  Subdivergences, with r
at loop positions p < q:
  q = p+1:        SE on l_p:              A(G) - A(G/S) U(S)
  q = p+2:        vertex around a_k:      A(G) - A(G/V) U(V)
  (p,q) = (1,5):  vertex around mu:       A(G) - A(G/Vm) U(Vm) - L(G/Vm) A(Vm) + U(G/Vm) A(Vm)
  gap 3:          A(G)
The zero-momentum subtraction applies to every node that carries the loop (sign flipped).

The point function is split by r placement into 10 compilation units l4_bp<p><q>.c (no graph
entry) summed by l4_b.c, so no file is too large for the compiler.

    python3 lbl4b.py OUTDIR [JOBS]    -> l4_b.c, l4_bp12.c ... l4_bp45.c (tables in cache_lbl.pkl)
"""
import itertools, math, os, sys, pickle, multiprocessing as mp, time
from fractions import Fraction
import build, dirac
from lbl4c import node, tkey, table, v1_node, KA, KB, add, sub

NE, NL = 2, 6
LP0, AP0, RP = 2, 8, 11
NPAR = 12
PAIRS = [(p, q) for p in range(1, 6) for q in range(p + 1, 6)]
PERMS = list(itertools.permutations([1, 2, 3]))

def sequence(p, q, perm):
    """Loop sequence: seq[0]='mu', 'r' at p and q, a-photons perm[0..2] on the other positions."""
    s = ['mu'] + [None]*5; s[p] = s[q] = 'r'
    rest = [x for x in range(1, 6) if x not in (p, q)]
    for x, a in zip(rest, perm): s[x] = a
    return s

def loop_node(name, kind, seq, zero, v0=False, drop_v=(), drop_l=(), se_line=None):
    """Node containing the loop (and the open line).  drop_v: loop positions removed; drop_l: loop
    lines removed; se_line: loop line replaced by the SE slot (its r ends also dropped)."""
    has_r = 'r' in [seq[x] for x in range(6) if x not in drop_v]
    cols = [0, 1, 2] + ([3] if has_r else [])
    lines = [0, 1] + [LP0 + q for q in range(6) if q not in drop_l and q != se_line]
    slot = {g: t for t, g in enumerate(lines)}
    nl = len(lines)
    A = []; B = []
    # open line: e0 = p - k1, e1 = p - k1 - k2
    oa, ob = [0, 0, 0, 0], [1, 0]
    for e in range(NE):
        oa = sub(oa, KA[e + 1]); A.append([oa[c] for c in cols]); B.append([1, 0] if v0 else list(ob))
    # loop lines
    # r flows only through its own arc: for r around mu (ends at 1 and 5) it starts on l0
    la, lb = [0, 0, 1, (1 if seq[1] == seq[5] == 'r' else 0)], [0, 0]; seen_r = False; rows = []
    for q in range(6):
        rows.append((list(la), list(lb)))
        x = seq[(q + 1) % 6]
        if q < 5:
            if x == 'r':
                la = add(la, KA['r']) if seen_r else sub(la, KA['r']); seen_r = True
            elif not zero:
                la = add(la, KA[x]); lb = add(lb, KB[x])
    for q in range(6):
        if q in drop_l or q == se_line: continue
        a, b = rows[q]; A.append([a[c] for c in cols]); B.append([0, 0] if v0 else b)
    photons = [1, 2, 3] + (['r'] if has_r else [])
    for x in photons:
        A.append([KA[x][c] for c in cols]); B.append([0, 0] if v0 else list(KB[x]))
    gidx = lines + [RP if x == 'r' else AP0 + x - 1 for x in photons]
    chain = [('v', 'a3'), ('s', slot[1]), ('v', 'a2'), ('s', slot[0]), ('v', 'a1')]
    vname = lambda x: 'mu' if x == 'mu' else ('g' if x == 'r' else 'a%d' % x)
    toks = []
    for q in range(5, -1, -1):
        if q == se_line: toks.append(('s', nl))
        elif q not in drop_l: toks.append(('s', slot[LP0 + q]))
        if q not in drop_v: toks.append(('v', vname(seq[q])))
    selist = ['S'] if se_line is not None else ()
    slotline = list(range(nl)) + [slot[LP0 + se_line - 1]] if se_line is not None else None
    return node(name, kind, chain, toks, A, B, gidx, nl, selist=selist, slotline=slotline)

def terms_for(p, q, perm, zero):
    """[(sign, [(node, k)], SE node or None)] for one sequence."""
    seq = sequence(p, q, perm)
    tag = '%d%d_%s%s' % (p, q, ''.join(map(str, perm)), '_z' if zero else '')
    T = [(1, [(loop_node('n_b%s_G' % tag, 'F2', seq, zero), 0)], None)]
    if q == p + 1:
        red = loop_node('n_b%s_GS' % tag, 'F2', seq, zero, drop_v=(p, q), se_line=p)
        S = node('n_b%d%d_S' % (p, q), 'SE', [('v', 'g'), ('s', 0), ('v', 'g')], None, [[-1], [1]], [[1, 0], [0, 0]], [LP0 + p, RP], 1)
        T.append((-1, [(red, 0)], S))
    elif q == p + 2:
        red = loop_node('n_b%s_GV' % tag, 'F2', seq, zero, drop_v=(p, q), drop_l=(p, p + 1))
        V = v1_node('n_b%d%d_V' % (p, q), 'V0', LP0 + p, LP0 + p + 1)
        T.append((-1, [(red, 0), (V, 0)], None))
    elif (p, q) == (1, 5):
        RF = loop_node('n_b%s_RF' % tag, 'F2', seq, zero, drop_v=(1, 5), drop_l=(5, 0))
        RV = loop_node('n_b%s_RV' % tag, 'V0', seq, zero, v0=True, drop_v=(1, 5), drop_l=(5, 0))
        VA = v1_node('n_b15_VA', 'F2', LP0 + 5, LP0); VU = v1_node('n_b15_VU', 'V0', LP0 + 5, LP0)
        T += [(-1, [(RF, 0), (VU, 0)], None), (-1, [(RV, 1), (VA, 0)], None), (1, [(RV, 0), (VA, 0)], None)]
    return T

def all_nodes(p, q):
    out = {}
    for perm in PERMS:
        for zero in (False, True):
            for sg, prod, se in terms_for(p, q, perm, zero):
                for nd, k in prod: out[nd.name] = nd
                if se is not None: out[se.name] = se
    return out

def _tab(nd):
    t = time.time(); tb = dirac.table(nd.chain, len(nd.slots), nd.kind, loop=nd.loop)
    return tkey(nd), tb, time.time() - t

def part_code(cache, p, q):
    """Point-function part for r at (p, q): all 6 a-orderings and their zero subtractions."""
    import lblcommon
    funcs = []; code = []
    for nm, nd in all_nodes(p, q).items():
        fn, call = lblcommon.node_function(nd, table(nd, cache), NPAR); funcs.append(fn); code.append(call)
    em = build.Emitter(); p0 = NPAR - 2*4
    em.lines.append("  REAL tot=0, tota=0;")
    for perm in PERMS:
        for zero in (False, True):
            for sg, prod, se in terms_for(p, q, perm, zero):
                s = sg*(-1 if zero else 1)
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
    pt = ["REAL FN(bp%d%d)(const REAL *x, const REAL Vx, REAL *ab){" % (p, q)]
    pt += code + em.lines + ["  *ab=tota; return tot;\n}"]
    return "\n".join(funcs) + "\n" + "\n".join(pt)

def write_part(out, cache, p, q):
    body = part_code(cache, p, q)
    src = ["/* Generated by lbl4b.py: IV(b) point-function part, r at loop positions %d,%d (no graph entry) */" % (p, q),
           '#include "common.h"',
           '#define REAL double', '#define FN(x) l4_b_##x##_d', body, '#undef REAL', '#undef FN',
           '#define REAL QREAL', '#define FN(x) l4_b_##x##_q', body, '#undef REAL', '#undef FN']
    open(os.path.join(out, 'l4_bp%d%d.c' % (p, q)), 'w').write("\n".join(src) + "\n")

def write_main(out):
    cname = 'l4_b'
    decl = ["REAL FN(bp%d%d)(const REAL *x, const REAL Vx, REAL *ab);" % pq for pq in PAIRS]
    pt = ["static REAL FN(pt)(const REAL *x, const REAL Vx, REAL *ab){", "  REAL tot=0, tota=0, a;"]
    pt += ["  tot+=FN(bp%d%d)(x,Vx,&a); tota+=a;" % pq for pq in PAIRS]
    pt += ["  *ab=tota; return tot;\n}"]
    # ext class {l5, l0} (the loop lines at mu, GL), photon a2 fixed = 1, all other parameters free.
    # a2 is the photon the open-line mirror (a1 <-> a3, with the loop reversed) keeps in place, so the
    # driver can fold the sum over the mirror as well as over the loop reversal (GEN_LBLCW0=1: a1,
    # the kernels before 2026-10-26).
    ea, eb = LP0 + 5, LP0
    fixed = AP0 if os.environ.get('GEN_LBLCW0', '0') == '1' else AP0 + 1
    free = [i for i in range(NPAR) if i not in (ea, eb, fixed)]
    ngl = (NE + NL + 2)//2
    xs, ws = build.gauss_legendre01(ngl)
    glc = "static const double GLX_%s[%d]={%s}, GLW_%s[%d]={%s};" % (cname, ngl, ",".join(repr(float(x)) for x in xs), cname, ngl, ",".join(repr(float(w)) for w in ws))
    norm = -Fraction(-1, 4)**4/16          # fermion loop sign
    D = 1 + len(free)
    wr = ["static REAL FN(f)(const REAL *sf, REAL *ab){",
          "  REAL s[%d]; const REAL sig=sf[0]; s[%d]=1;" % (NPAR, fixed)]
    wr += ["  s[%d]=sf[%d];" % (i, d + 1) for d, i in enumerate(free)]
    wr += ["  REAL acc=0, aa=0, a;",
           "  for(int g=0;g<%d;g++){ const REAL u=(REAL)GLX_%s[g]; s[%d]=sig*u; s[%d]=sig*(1-u);" % (ngl, cname, ea, eb),
           "    const REAL wg=(REAL)GLW_%s[g]; acc+=wg*FN(pt)(s,0,&a); aa+=wg*a; }" % cname,
           "  const REAL nrm=%s*sig;" % build.cnum(norm),
           "  *ab=RABS(nrm)*aa; return nrm*acc;\n}"]
    body = "\n".join(decl + pt + wr)
    title = 'LbL with a photon on the loop (60 sequences)'
    src = ["/* Generated by lbl4b.py: %s */" % title, '#include "common.h"', glc,
           '#define REAL double', '#define FN(x) l4_b_##x##_d', body, '#undef REAL', '#undef FN',
           '#define REAL QREAL', '#define FN(x) l4_b_##x##_q', body, '#undef REAL', '#undef FN',
           'const graph_t G_%s={"%s","%s",4,%d,1,%d,%d,l4_b_f_d,l4_b_f_q,l4_b_pt_d,l4_b_pt_q};' %
           (cname, cname, title, D, NPAR, fixed)]
    open(os.path.join(out, cname + '.c'), 'w').write("\n".join(src) + "\n")

def needed(cache):
    need = {}
    for p, q in PAIRS:
        for nd in all_nodes(p, q).values():
            if tkey(nd) not in cache: need[tkey(nd)] = nd
    return need

_W = None
def _write_part(pq): return write_part(_W[0], _W[1], *pq)

def main(out, jobs):
    import lblcommon; lblcommon.use_form()
    cf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache_lbl.pkl')
    cache = pickle.load(open(cf, 'rb')) if os.path.exists(cf) else {}
    need = needed(cache)
    print('tables needed:', len(need), flush=True)
    if need:
        with mp.Pool(jobs) as pool:
            for k, tb, dt in pool.imap_unordered(_tab, sorted(need.values(), key=lambda n: -len(n.slots))):
                cache[k] = tb; print('  table %s: %d structures, %.0fs' % (k[:8], len(tb), dt), flush=True)
                c2 = pickle.load(open(cf, 'rb')) if os.path.exists(cf) else {}
                c2.update(cache); pickle.dump(c2, open(cf + '.tmp', 'wb')); os.replace(cf + '.tmp', cf)
    global _W; _W = (out, cache)                          # the forked workers see it without pickling
    with mp.get_context('fork').Pool(jobs) as pool: pool.map(_write_part, PAIRS)
    write_main(out)
    import lblvp; lblvp.write_table(out)
    print('l4_b + %d parts' % len(PAIRS))

if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 4)
