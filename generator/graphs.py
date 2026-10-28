"""
Vertex graphs without lepton loops, their UV-divergent subgraphs, forests, and Volkov's
subtraction terms (PRD 96, 096018; PRD 98, 076018):

    A~_G = sum_{F in Forests[G], G' in I[G] & F} (-1)^{|F|-1} M^G_{F,G'}

F: forest (pairwise nested-or-disjoint set) of UV-divergent 1PI subgraphs with G in F.
I[G]: vertexlike subgraphs containing the external vertex.
M^G_{F,G'} (Volkov eq. 11): A on G'; L on G'' in F with G' < G'' < G; (L - U) on G if G' != G;
U on the rest.

A graph of order L is an electron chain with positions 0..2L; the external photon sits at
position `ext`, the other 2L positions are paired by photon arcs.  Electron line i joins
positions i and i+1; photon line k is arc k.  Parameter index: electrons 0..2L-1, photons 2L+k.
"""
import itertools

def pairings(pts):
    if not pts: yield (); return
    a = pts[0]
    for i in range(1, len(pts)):
        rest = pts[1:i] + pts[i+1:]
        for m in pairings(rest): yield ((a, pts[i]),) + m

class Graph:
    def __init__(self, L, ext, arcs):
        self.L, self.ext = L, ext
        self.arcs = tuple(sorted(tuple(sorted(a)) for a in arcs))
        self.npos = 2*L + 1
        self.partner = {}
        for k, (a, b) in enumerate(self.arcs):
            self.partner[a] = (k, b); self.partner[b] = (k, a)
        self.ne = 2*L; self.nph = L; self.n = 3*L
    def covered(self, i, arcs=None):
        arcs = self.arcs if arcs is None else arcs
        return any(a <= i < b for a, b in arcs)
    def is_1pi(self):
        return all(self.covered(i) for i in range(self.ne))
    def mirror(self):
        m = lambda p: 2*self.L - p
        return Graph(self.L, m(self.ext), [(m(a), m(b)) for a, b in self.arcs])
    def key(self): return (self.ext, self.arcs)
    def volkov(self):
        """Volkov's notation: '<ext>; a-b, ...' with positions counted from 1."""
        return "%d; " % (self.ext+1) + ", ".join("%d-%d" % (a+1, b+1) for a, b in self.arcs)
    # ---------------- subgraphs ----------------
    def segment(self, u, v):
        """Subgraph on positions u..v (all arcs inside).  Returns dict or None if not a
        UV-divergent 1PI subgraph."""
        inside = [k for k, (a, b) in enumerate(self.arcs) if u <= a and b <= v]
        arcs = [self.arcs[k] for k in inside]
        for p in (u, v):
            if p == self.ext or p not in self.partner: return None
            k, q = self.partner[p]
            if not (u <= q <= v): return None
        extph = [p for p in range(u+1, v) if p == self.ext or not (u <= self.partner[p][1] <= v)]
        if len(extph) > 1: return None
        if not all(self.covered(i, arcs) for i in range(u, v)): return None
        return dict(u=u, v=v, arcs=tuple(inside), kind='SE' if not extph else 'V',
                    extpos=extph[0] if extph else None,
                    hasext=(u < self.ext < v))
    def subgraphs(self):
        out = []
        for u in range(self.npos):
            for v in range(u+1, self.npos):
                s = self.segment(u, v)
                if s is not None: out.append((u, v))
        return out          # includes (0, 2L) = G itself

def enumerate_graphs(L):
    """All 1PI vertex graphs without lepton loops; returns list of (rep, weight) with one
    representative per mirror pair (weight 2) or self-mirror graph (weight 1)."""
    npos = 2*L + 1; seen = set(); out = []
    for ext in range(npos):
        pts = [p for p in range(npos) if p != ext]
        for pr in pairings(pts):
            g = Graph(L, ext, pr)
            if not g.is_1pi(): continue
            if g.key() in seen: continue
            m = g.mirror(); seen.add(g.key()); seen.add(m.key())
            out.append((g, 1 if m.key() == g.key() else 2))
    return out

def all_graphs(L):
    npos = 2*L + 1; out = []
    for ext in range(npos):
        pts = [p for p in range(npos) if p != ext]
        for pr in pairings(pts):
            g = Graph(L, ext, pr)
            if g.is_1pi(): out.append(g)
    return out

# ---------------- forests and terms ----------------
def contains(a, b):     # segment a contains segment b (as subgraphs: nested segments)
    return a[0] <= b[0] and b[1] <= a[1]
def disjoint(a, b):
    return a[1] < b[0] or b[1] < a[0]
def compatible(a, b):
    return contains(a, b) or contains(b, a) or disjoint(a, b)

def forests(g):
    subs = g.subgraphs(); G = (0, 2*g.L)
    others = [s for s in subs if s != G]
    out = []
    for r in range(len(others)+1):
        for comb in itertools.combinations(others, r):
            if all(compatible(a, b) for a, b in itertools.combinations(comb, 2)):
                out.append((G,) + comb)
    return out

def terms(g):
    """Expanded subtraction terms: list of (sign, {segment: op}) with op in 'A','U','L'."""
    G = (0, 2*g.L); res = []
    for F in forests(g):
        n = len(F)
        for Gp in F:
            s = g.segment(*Gp)
            if not (s['kind'] == 'V' and s['hasext']): continue
            sign = (-1)**(n-1)
            fixed = {}; lu = []
            for Gs in F:
                if Gs == Gp: fixed[Gs] = 'A'
                elif contains(Gs, Gp) and Gs == G: lu.append(Gs)     # (L - U) only on G itself
                elif contains(Gs, Gp): fixed[Gs] = 'L'              # Volkov (11): L on G'' in I[G], G' < G'' < G
                else: fixed[Gs] = 'U'
            for choice in itertools.product(['L', 'U'], repeat=len(lu)):
                ops = dict(fixed); sg = sign
                for Gs, c in zip(lu, choice):
                    ops[Gs] = c
                    if c == 'U': sg = -sg
                res.append((sg, ops))
    return res

def children(seg, F):
    """Maximal elements of F strictly inside seg."""
    inner = [s for s in F if s != seg and contains(seg, s)]
    return tuple(sorted(s for s in inner if not any(t != s and contains(t, s) for t in inner)))

if __name__ == '__main__':
    for L in (1, 2, 3):
        a = all_graphs(L); r = enumerate_graphs(L)
        print("L=%d: %d graphs, %d up to mirror" % (L, len(a), len(r)))
    for g, w in enumerate_graphs(2):
        print(g.volkov(), 'w', w, 'subs', g.subgraphs(), 'terms', len(terms(g)))
        for sg, ops in terms(g): print('   ', sg, ops)
    import collections
    c = collections.Counter(len(terms(g)) for g, w in enumerate_graphs(3))
    print('3-loop term counts:', sorted(c.items()))
