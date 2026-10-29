"""
Three-loop light-by-light graphs (electron loop, 6 photon orderings), summed pointwise.

Open electron line: p -> o1 -e0-> o2 -e1-> o3 -> p'; photon a (a = 1,2,3) joins o_a to a loop vertex.
Electron loop: vertices L0 (external photon, index mu), L1, L2, L3; line l_i runs L_i -> L_{i+1}.
The 6 orderings (photon a at loop vertex perm[a]) share one labelling of the Schwinger parameters
    0:e0 1:e1 2:l0 3:l1 4:l2 5:l3 6:a1 7:a2 8:a3,
so the log UV divergence of the individual loops (a momentum-independent 4-photon tensor that
vanishes in the sum over orderings) cancels point by point.  No subgraph subtractions are needed.
Photon momenta k1, k2 are loop momenta, k3 = -k1 - k2 + p - p'; loop momentum l on l0.
"""
import itertools, math, os, sys, pickle, hashlib, multiprocessing as mp, time
from fractions import Fraction
import build, dirac

class LNode:
    def __init__(self, perm, zero=False):
        self.perm = perm; self.zero = zero                     # perm[a] = loop vertex (1..3) of photon a (0-based a)
        kA = {0: [1, 0, 0], 1: [0, 1, 0], 2: [-1, -1, 0]}
        kB = {0: [0, 0], 1: [0, 0], 2: [1, -1]}
        A = [[-1, 0, 0], [-1, -1, 0]]; B = [[1, 0], [1, 0]]
        inv = {perm[a]: a for a in range(3)}
        la, lb = [0, 0, 1], [0, 0]
        for i in range(4):
            A.append(list(la)); B.append(list(lb))
            if i < 3 and not zero:
                a = inv[i+1]; la = [x+y for x, y in zip(la, kA[a])]; lb = [x+y for x, y in zip(lb, kB[a])]
        for a in range(3): A.append(kA[a]); B.append(kB[a])
        self.A, self.B = A, B
        self.m2 = [1]*6 + [0]*3
        self.gidx = list(range(9)); self.slotline = list(range(6)); self.selist = []; self.egaps = list(range(6))
        self.kind = 'F2'; self.L = 3; self.seg = ('lbl',) + tuple(perm); self.kids = ()
        self.slots = [('e', i) for i in range(6)]
        self.chain = [('v', 'a3'), ('s', 1), ('v', 'a2'), ('s', 0), ('v', 'a1')]
        vn = lambda j: 'a%d' % (inv[j]+1)
        self.loop = [('s', 5), ('v', vn(3)), ('s', 4), ('v', vn(2)), ('s', 3), ('v', vn(1)), ('s', 2), ('v', 'mu')]
        self.name = 'n_lbl_%d%d%d' % tuple(p for p in perm) + ('_z' if zero else '')
    def table(self, cache):
        key = repr((self.chain, len(self.slots), self.kind)) + repr(self.loop)
        h = hashlib.md5(key.encode()).hexdigest()
        if h not in cache: cache[h] = dirac.table(self.chain, 6, 'F2', loop=self.loop)
        return cache[h]
    def key(self):
        return hashlib.md5((repr((self.chain, len(self.slots), self.kind)) + repr(self.loop)).encode()).hexdigest()

SUBTRACT = os.environ.get('LBL_NOSUB') is None
PERMS = [p for p in itertools.permutations([1, 2, 3])]

def _tab(nd):
    t = time.time(); tb = dirac.table(nd.chain, 6, 'F2', loop=nd.loop)
    return nd.key(), tb, time.time() - t

def generate(out, cache, cachefile, jobs):
    nodes = [LNode(p) for p in PERMS]
    znodes = [LNode(p, zero=True) for p in PERMS] if SUBTRACT else []
    need = [nd for nd in nodes if nd.key() not in cache]
    if need:
        with mp.get_context('fork').Pool(jobs) as pool:
            for k, tb, dt in pool.imap_unordered(_tab, need):
                cache[k] = tb; print('  lbl table %s: %d structures, %.0fs' % (k[:8], len(tb), dt), flush=True)
                pickle.dump(cache, open(cachefile + '.tmp', 'wb')); os.replace(cachefile + '.tmp', cachefile)
    cname = 'g3_lbl'
    pt = pt_code(cache, cname)
    return write_wrapper(out, cname, pt)

def pt_code(cache, cname):
    """Point function pt_<cname>(x, Vx, &ab) of the six orderings (and their zero-momentum
    subtractions), summed; the tables must be in cache."""
    nodes = [LNode(p) for p in PERMS]
    znodes = [LNode(p, zero=True) for p in PERMS] if SUBTRACT else []
    code = []; lines = ["  REAL tot=0, tota=0;"]
    p0 = 9 - 6
    for nd in nodes + znodes:
        code.append(build.emit_node(nd, nd.table(cache)))
        assert nd.maxw < p0
        lines.append("  { const REAL iV=1/(%s_V+Vx); REAL r=0, ra=0, ivp=1;" % nd.name)
        for w in range(p0-1, -1, -1):
            lines.append("    ivp*=iV;")
            if w < build.NW:
                f = math.factorial(p0-w-1)
                lines.append("    r+=%d*%s_W[0][0][%d]*ivp; ra+=%d*%s_Wa[0][0][%d]*ivp;" % (f, nd.name, w, f, nd.name, w))
        lines.append("    tot%s=r; tota+=ra; }" % ('-' if nd.zero else '+'))
    pt = ["static REAL FN(pt_%s)(const REAL *x, const REAL Vx, REAL *ab){" % cname]
    pt += ["  const REAL x%d=x[%d];" % (i, i) for i in range(9)]
    pt += code + lines + ["  *ab=tota; return tot;\n}"]
    return "\n".join(pt)

def write_wrapper(out, cname, pt):
    # wrapper: ext class {l0, l3} (GL, 4 points), photon a2 fixed = 1, free: e0 e1 l1 l2 a1 a3.  a2 is
    # the photon the open-line mirror (a1 <-> a3) keeps in place, so the driver folds the sum over the
    # mirror as well as the loop reversal (group of order 4).  GEN_LBLCW0=1: a1 fixed, the kernel
    # before 2026-10-27.
    cw0 = os.environ.get('GEN_LBLCW0', '0') == '1'
    xs, ws = build.gauss_legendre01(4)
    glc = "static const double GLX_%s[4]={%s}, GLW_%s[4]={%s};" % (cname, ",".join(repr(float(x)) for x in xs), cname, ",".join(repr(float(w)) for w in ws))
    norm = -Fraction(-1, 4)**3/16          # fermion loop sign
    wr = ["static REAL FN(f_%s)(const REAL *sf, REAL *ab){" % cname,
          ("  REAL s[9]; const REAL sig=sf[0]; s[0]=sf[1]; s[1]=sf[2]; s[3]=sf[3]; s[4]=sf[4]; s[6]=1; s[7]=sf[5]; s[8]=sf[6];" if cw0 else
           "  REAL s[9]; const REAL sig=sf[0]; s[0]=sf[1]; s[1]=sf[2]; s[3]=sf[3]; s[4]=sf[4]; s[6]=sf[5]; s[7]=1; s[8]=sf[6];"),
          "  REAL acc=0, aa=0, a;",
          "  for(int g=0;g<4;g++){ const REAL u=(REAL)GLX_%s[g]; s[2]=sig*u; s[5]=sig*(1-u);" % cname,
          "    const REAL wg=(REAL)GLW_%s[g]; acc+=wg*FN(pt_%s)(s,0,&a); aa+=wg*a; }" % (cname, cname),
          "  const REAL nrm=%s*sig;" % build.cnum(norm),
          "  *ab=RABS(nrm)*aa; return nrm*acc;\n}"]
    wr = "\n".join(wr)
    src = ["/* Generated by lbl.py: three-loop light-by-light graphs (6 photon orderings summed pointwise). */",
           '#include "common.h"', glc,
           '#define REAL double', '#define FN(x) x##_d', pt, wr, '#undef REAL', '#undef FN',
           '#define REAL QREAL', '#define FN(x) x##_q', pt, wr, '#undef REAL', '#undef FN',
           'const graph_t G_%s={"%s","LbL (6 orderings)",3,7,1,9,%d,f_%s_d,f_%s_q,pt_%s_d,pt_%s_q};' % (cname, cname, 6 if cw0 else 7, cname, cname, cname, cname)]
    open(os.path.join(out, cname + '.c'), 'w').write("\n".join(src) + "\n")
    return cname

if __name__ == '__main__':
    out = sys.argv[1]
    cf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache_lbl.pkl')
    cache = pickle.load(open(cf, 'rb')) if os.path.exists(cf) else {}
    print(generate(out, cache, cf, int(sys.argv[2]) if len(sys.argv) > 2 else 4))
