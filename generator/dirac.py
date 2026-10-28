"""
Dirac-structure tables for vertex / self-energy (reduced) graphs in Schwinger-parametric form.

A chain is a list of tokens in matrix order (left = outgoing electron):
    ('v', name)   vertex gamma; an internal photon name appears twice, 'mu' = external index
    ('s', j)      numerator slot j: an electron line (p-slash + m) or an SE insertion (u + v p-slash)
Each slot is expanded into choices; line momentum = Q_j + a_j p + b_j p'.
    mode 'F2' : choices 'p','pp','1','Q' ; value = dN/dt at t=0 of the F2 projector numerator
                N(t) = (4-t) Tr[(p'+1) G^mu (p+1) gamma_mu] - (2+t) Tr[(p'+1) G^mu (p+1)] (p+p')_mu,
                F2(0) = N'(0)/16.  Integer.
    mode 'V0' : q = 0, p = p', p^2 = 1 ; choices 'P','1','Q' ; Gamma^mu = a g^mu + b p^mu + c p^mu pslash
                + d [pslash, g^mu] ; value = (a, a+b+c)   (U and L operators)
    mode 'SE' : p^2 = 1 ; choices 'P','1','Q' ; Sigma = u + v pslash ; value = (u, v)
Q choices are Wick-paired in all perfect matchings; the pairing factor (-1/2)^w is NOT included
(the caller multiplies it in, together with the contraction coefficients C_ij).
"""
import numpy as np, itertools
from fractions import Fraction

I2 = np.eye(2); Z2 = np.zeros((2, 2))
sx = np.array([[0, 1], [1, 0]]); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]])
g0 = np.block([[I2, Z2], [Z2, -I2]]).astype(complex)
GU = np.array([g0] + [np.block([[Z2, s], [-s, Z2]]).astype(complex) for s in (sx, sy, sz)])
eta = np.array([1., -1, -1, -1])
GL = np.array([eta[m]*GU[m] for m in range(4)])
One = np.eye(4, dtype=complex)
def slash(v): return sum(eta[m]*v[m]*GU[m] for m in range(4))

def matchings(l):
    if not l: yield (); return
    a = l[0]
    for k in range(1, len(l)):
        for m in matchings(l[1:k] + l[k+1:]): yield ((a, l[k]),) + m

def combos(nslots, choices, maxQ=None):
    for c in itertools.product(choices, repeat=nslots):
        Qs = [i for i in range(nslots) if c[i] == 'Q']
        if len(Qs) % 2: continue
        if maxQ is not None and len(Qs) > maxQ: continue
        for m in matchings(Qs): yield c, m

def _batch(chain, combo, pr, mats, loop=None):
    """mats: dict choice -> (T,4,4) array (T kinematic points). Returns M (T,B,4,4), idx, col."""
    wick = {}
    for n, (a, b) in enumerate(pr): wick[a] = 'w%d' % n; wick[b] = 'w%d' % n
    allt = chain + (loop or [])
    vnames = sorted(set(t[1] for t in allt if t[0] == 'v' and t[1] != 'mu'))
    hasmu = any(t == ('v', 'mu') for t in allt)
    names = vnames + sorted(set(wick.values())) + (['mu'] if hasmu else [])
    idx = np.array(list(itertools.product(range(4), repeat=len(names))), dtype=np.int64).reshape(-1, len(names))
    col = {n: i for i, n in enumerate(names)}
    T = next(iter(mats.values())).shape[0]; B = idx.shape[0]
    seen = set()
    def mul(M, X):
        return X.copy() if M is None else M @ X
    def product(toks):
        M = None
        for tok in toks:
            if tok[0] == 'v':
                nm = tok[1]
                if nm == 'mu': M = mul(M, np.broadcast_to(GU[idx[:, col['mu']]][None], (T, B, 4, 4))); continue
            else:
                j = tok[1]; ch = combo[j]
                if ch != 'Q': M = mul(M, np.broadcast_to(mats[ch][:, None], (T, B, 4, 4))); continue
                nm = wick[j]
            g = GL if nm in seen else GU; seen.add(nm)
            M = mul(M, np.broadcast_to(g[idx[:, col[nm]]][None], (T, B, 4, 4)))
        if M is None: M = np.broadcast_to(One, (T, B, 4, 4)).copy()
        return M
    M = product(chain)
    if loop is not None:
        Ml = product(loop)
        M = M * np.einsum('tbii->tb', Ml)[..., None, None]
    return M, idx, col

def _rat(x, den=4096, tol=1e-7):
    f = Fraction(x).limit_denominator(den)
    assert abs(float(f) - x) < tol*max(1, abs(x)), (x, f)
    return f

MO_OFF = False  # True: plain left-to-right contraction (cross-check of the middle-out order)
SLOW = False   # True: brute-force index enumeration for F2 (slower; kept as a cross-check)

USE_FORM = False   # set by the generators: tables come from FORM (formdirac.py), same format

def table(chain, nslots, mode, maxQ=None, loop=None):
    """loop: optional token list of a closed fermion loop (trace taken; sign not included)."""
    if USE_FORM and maxQ is None:
        import formdirac
        return formdirac.table(chain, nslots, mode, loop)
    out = {}
    if mode == 'F2' and not SLOW:
        return table_f2_fast(chain, nslots, maxQ, loop)
    if mode == 'F2':
        deg = nslots + 4
        xs = np.cos(np.pi*(np.arange(deg+3)+0.5)/(deg+3))    # Chebyshev nodes in [-1,1]
        ts = 2.0*xs                                           # t in [-2, 2]
        p = np.array([1, 0, 0, 0], complex)
        E = 1 - ts/2; k = np.sqrt(E*E - 1 + 0j)
        pp = np.stack([E, 0*E, 0*E, k], axis=1)
        ps = np.broadcast_to(slash(p), (len(ts), 4, 4))
        pps = np.array([slash(v) for v in pp])
        mats = {'p': ps, 'pp': pps, '1': np.broadcast_to(One, (len(ts), 4, 4))}
        P = p[None, :] + pp                                  # (T,4)
        V = np.polynomial.polynomial.polyvander(xs, deg)
        for combo, pr in combos(nslots, ['p', 'pp', '1', 'Q'], maxQ):
            M, idx, col = _batch(chain, combo, pr, mats, loop)
            M = (pps + One)[:, None] @ M @ (ps + One)[:, None]
            mu = idx[:, col['mu']]
            Tg = np.einsum('tbij,bji->t', M, GL[mu])
            TP = np.einsum('tbii,tb->t', M, eta[mu][None, :]*P[:, mu])
            N = (4-ts)*Tg - (2+ts)*TP
            assert np.max(np.abs(N.imag)) < 1e-7*max(1, np.max(np.abs(N))), combo
            c = np.linalg.lstsq(V, N.real, rcond=None)[0]
            d1 = c[1]/2.0                                    # d/dt = (1/2) d/dx
            sc = max(1, np.max(np.abs(N.real)))
            assert abs(c[0]) < 1e-7*sc, (combo, pr, c[0])
            r = round(d1)
            assert abs(d1 - r) < 1e-6*max(1, abs(d1)), (combo, pr, d1)
            if r != 0: out[(combo, pr)] = int(r)
        return out
    p = np.array([1, 0, 0, 0], complex); ps = slash(p)
    mats = {'P': ps[None], '1': One[None]}
    if mode == 'V0':
        # basis for Gamma^mu (upper index): g^mu, p^mu 1, p^mu pslash, [pslash, g^mu]
        basis = np.stack([np.stack([GU[m], p[m]*One, p[m]*ps, ps@GU[m]-GU[m]@ps]) for m in range(4)])  # (4mu,4b,4,4)
        Bm = basis.transpose(1, 0, 2, 3).reshape(4, -1).T                                              # (64,4)
        for combo, pr in combos(nslots, ['P', '1', 'Q'], maxQ):
            M, idx, col = _batch(chain, combo, pr, mats, loop)
            mu = idx[:, col['mu']]
            Gam = np.zeros((4, 4, 4), complex)
            np.add.at(Gam, mu, M[0])
            sol, res, rk, sv = np.linalg.lstsq(Bm, Gam.reshape(-1), rcond=None)
            assert np.max(np.abs(Bm@sol - Gam.reshape(-1))) < 1e-8*max(1, np.max(np.abs(Gam))), combo
            assert np.max(np.abs(sol.imag)) < 1e-9*max(1, np.max(np.abs(sol))), combo
            a, b, c, d = [_rat(x) for x in sol.real]
            if a != 0 or a+b+c != 0: out[(combo, pr)] = (a, a+b+c)
        return out
    if mode == 'SE':
        for combo, pr in combos(nslots, ['P', '1', 'Q'], maxQ):
            M, idx, col = _batch(chain, combo, pr, mats, loop)
            S = M[0].sum(axis=0)
            u = np.trace(S)/4; v = np.trace(S@ps)/4
            assert np.max(np.abs(S - u*One - v*ps)) < 1e-8*max(1, np.max(np.abs(S))), combo
            assert abs(u.imag) < 1e-9 and abs(v.imag) < 1e-9
            u, v = _rat(u.real), _rat(v.real)
            if u != 0 or v != 0: out[(combo, pr)] = (u, v)
        return out
    raise ValueError(mode)

if __name__ == '__main__':
    import time
    # one-loop vertex: chain  gamma_a (p'-k) gamma^mu (p-k) gamma_a
    ch1 = [('v', 'a'), ('s', 1), ('v', 'mu'), ('s', 0), ('v', 'a')]
    print('F2 one-loop:', table(ch1, 2, 'F2'))
    print('V0 one-loop:', table(ch1, 2, 'V0'))
    print('SE one-loop:', table([('v', 'a'), ('s', 0), ('v', 'a')], 1, 'SE'))
    t = time.time()
    ch16 = [('v','c'),('s',5),('v','a'),('s',4),('v','b'),('s',3),('v','c'),('s',2),('v','mu'),('s',1),('v','b'),('s',0),('v','a')]
    S = table(ch16, 6, 'F2')
    print('g16 structures', len(S), '%.0fs' % (time.time()-t))

# ---------------- fast F2 tables: tensor-network contraction (no index enumeration) ----------------
def _chain_open(toks, combo, wick, ps, pps, X, names):
    """Sequential product keeping unmatched Lorentz indices open: X[t, a1..ak, i, j]."""
    T = X.shape[0]
    for tok in toks:
        if tok[0] == 'v': nm = tok[1]
        else:
            ch = combo[tok[1]]
            if ch == 'p': X = X @ ps; continue
            if ch == 'pp': X = X @ pps.reshape((T,) + (1,)*len(names) + (4, 4)); continue
            if ch == '1': continue
            nm = wick[tok[1]]
        if nm in names:                       # contract with gamma_alpha (lower index)
            k = names.index(nm)
            X = np.einsum('...aij,ajk->...ik', np.moveaxis(X, 1+k, -3), GL)
            names = names[:k] + names[k+1:]
        else:                                 # open a new index (appended last)
            X = np.einsum('...ij,ajk->...aik', X, GU)
            names = names + [nm]
    return X, names

def _order(toks, tokname, mo):
    """Contraction order: list of (position, left?) and its cost (sum of 4^open indices).
    mo=False: plain left to right; mo=True: middle-out from 'mu', greedily closing indices."""
    n = len(toks)
    if not mo:
        steps = [(i, False) for i in range(n)]
    else:
        start = next((i for i, t in enumerate(toks) if t == ('v', 'mu')), 0)
        steps = [(start, False)]; l, r = start, start + 1; op = {tokname(toks[start])} - {None}
        while l > 0 or r < n:
            cand = []
            for side, pos in ((0, l-1), (1, r)):
                if 0 <= pos < n and (side == 0 or pos < n):
                    if side == 0 and l == 0: continue
                    if side == 1 and r >= n: continue
                    nm = tokname(toks[pos]); cand.append((0 if nm is None else (-1 if nm in op else 1), side, pos))
            _, side, pos = min(cand)
            nm = tokname(toks[pos])
            if nm is not None: op ^= {nm}
            steps.append((pos, side == 0))
            if side == 0: l -= 1
            else: r += 1
    op, cost = set(), 0
    for pos, _ in steps:
        nm = tokname(toks[pos])
        if nm is not None: op ^= {nm}
        cost += 4**len(op)
    return steps, cost

def _chain_mo(toks, combo, wick, ps, pps):
    """Product of the chain tokens with open Lorentz indices, contracted in the cheaper of two
    orders: left to right, or middle-out from the external vertex 'mu' (growing the segment on the
    side that closes an index first, so nested photon pairs close as soon as both ends are in).
    Returns X[t, a1..ak, i, j], names."""
    T = pps.shape[0]
    def tokname(tok):
        if tok[0] == 'v': return tok[1]
        ch = combo[tok[1]]
        return wick[tok[1]] if ch == 'Q' else None
    def apply(X, names, tok, left):
        nd = len(names)
        if tok[0] == 's':
            ch = combo[tok[1]]
            if ch == '1': return X, names
            if ch in ('p', 'pp'):
                M = ps if ch == 'p' else pps.reshape((T,) + (1,)*nd + (4, 4))
                return (M @ X if left else X @ M), names
        nm = tokname(tok)
        if nm in names:
            k = names.index(nm)
            Y = np.moveaxis(X, 1+k, -3)
            X = np.einsum('ajk,...akl->...jl', GL, Y) if left else np.einsum('...aij,ajk->...ik', Y, GL)
            return X, names[:k] + names[k+1:]
        X = np.einsum('aij,...jk->...aik', GU, X) if left else np.einsum('...ij,ajk->...aik', X, GU)
        return X, names + [nm]
    s0, c0 = _order(toks, tokname, False)
    s1, c1 = _order(toks, tokname, True)
    steps = s1 if c1 < c0 else s0
    X = np.broadcast_to(One, (T, 4, 4)).copy(); names = []
    for pos, left in steps:
        X, names = apply(X, names, toks[pos], left)
    return X, names

def _f2_network(chain, combo, pr, loop, ps, pps, Peta):
    """(Tg(t), TP(t)) for one structure, by sequential contraction (no index enumeration)."""
    wick = {}
    for n, (a, b) in enumerate(pr): wick[a] = 'w%d' % n; wick[b] = 'w%d' % n
    T = pps.shape[0]
    if loop or MO_OFF:
        X, names = _chain_open(chain, combo, wick, ps, pps, (pps + One).copy(), [])
        X = X @ (ps + One)
    else:
        X, names = _chain_mo(chain, combo, wick, ps, pps)
        X = (pps + One).reshape((T,) + (1,)*len(names) + (4, 4)) @ X @ (ps + One)
    if loop:
        Y, names2 = _chain_open(loop, combo, wick, ps, pps, np.broadcast_to(One, (T, 4, 4)).copy(), [])
        tr = np.einsum('...ii->...', Y)
        for k, nm in enumerate(names2):          # lower the loop copy of indices shared with the line
            if nm in names:
                shp = [1]*tr.ndim; shp[1+k] = 4; tr = tr*eta.reshape(shp)
        let = dict(zip(sorted(set(names) | set(names2)), 'abcdefghklmnopq'))
        outn = [n for n in names if n not in names2] + [n for n in names2 if n not in names]
        X = np.einsum('t%sij,t%s->t%sij' % (''.join(let[n] for n in names), ''.join(let[n] for n in names2),
                      ''.join(let[n] for n in outn)), X, tr)
        names = outn
    assert names == ['mu'], names
    Tg = np.einsum('tmij,mji->t', X, GL)
    TP = np.einsum('tmii,tm->t', X, Peta)
    return Tg, TP

def table_f2_fast(chain, nslots, maxQ=None, loop=None):
    out = {}
    deg = nslots + 4
    xs = np.cos(np.pi*(np.arange(deg+3)+0.5)/(deg+3)); ts = 2.0*xs
    p = np.array([1, 0, 0, 0], complex)
    E = 1 - ts/2; k = np.sqrt(E*E - 1 + 0j)
    pp = np.stack([E, 0*E, 0*E, k], axis=1)
    ps = slash(p); pps = np.array([slash(v) for v in pp])
    Peta = (p[None, :] + pp)*eta[None, :]
    V = np.polynomial.polynomial.polyvander(xs, deg)
    for combo, pr in combos(nslots, ['p', 'pp', '1', 'Q'], maxQ):
        Tg, TP = _f2_network(chain, combo, pr, loop, ps, pps, Peta)
        N = (4-ts)*Tg - (2+ts)*TP
        assert np.max(np.abs(N.imag)) < 1e-7*max(1, np.max(np.abs(N))), combo
        c = np.linalg.lstsq(V, N.real, rcond=None)[0]
        d1 = c[1]/2.0; sc = max(1, np.max(np.abs(N.real)))
        assert abs(c[0]) < 1e-7*sc, (combo, pr, c[0])
        r = round(d1)
        assert abs(d1 - r) < 1e-6*max(1, abs(d1)), (combo, pr, d1)
        if r != 0: out[(combo, pr)] = int(r)
    return out

# ---------------- vacuum-polarization loops: k^mu k^nu coefficient ----------------
def table_pi(loop, nslots, maxQ=None):
    if USE_FORM and maxQ is None:
        import formdirac
        return formdirac.table_pi(loop, nslots)
    """Closed loop (trace) with external indices 'mu','nu' and momentum k flowing through.
    Slot choices 'k' (x k-slash), '1', 'Q'.  T^{mu nu}(k) = A g^{mu nu} + C(k^2) k^mu k^nu per structure;
    returns {structure: (C_0, C_1, ...)} with C(k^2) = sum_m C_m (k^2)^m (Minkowski k^2)."""
    out = {}
    M = nslots//2 + 1
    kap = np.array([0.5 + 0.35*j for j in range(M + 3)])
    T = len(kap)
    kmat = np.array([kv*GU[0] for kv in kap])            # k = (kappa,0,0,0): k-slash = kappa gamma^0
    ps = One
    V = np.vander(kap**2, M + 2, increasing=True)          # powers (kappa^2)^0..M+1
    for combo, pr in combos(nslots, ['k', '1', 'Q'], maxQ):
        wick = {}
        for n, (a, b) in enumerate(pr): wick[a] = 'w%d' % n; wick[b] = 'w%d' % n
        cmb = tuple('pp' if c == 'k' else c for c in combo)
        X, names = _chain_open(loop, cmb, wick, ps, kmat, np.broadcast_to(One, (T, 4, 4)).copy(), [])
        assert sorted(names) == ['mu', 'nu'], names
        Tt = np.einsum('tabii->tab', X)
        if names != ['mu', 'nu']: Tt = Tt.transpose(0, 2, 1)
        y = Tt[:, 0, 0] + Tt[:, 1, 1]                      # = C(kappa^2) kappa^2
        assert np.max(np.abs(y.imag)) < 1e-8*max(1, np.max(np.abs(y))), combo
        c = np.linalg.lstsq(V, y.real, rcond=None)[0]
        assert abs(c[0]) < 1e-7*max(1, np.max(np.abs(y))), (combo, c)
        Cs = [_rat(x) for x in c[1:]]
        if any(x != 0 for x in Cs): out[(combo, pr)] = tuple(Cs)
    return out
