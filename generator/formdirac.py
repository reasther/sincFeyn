"""
Dirac algebra for the Sinc QED node tables done by FORM instead of explicit 4x4 matrices.

A chain is the same token list dirac.py uses (left = outgoing electron):
    ('v', name)  vertex gamma; an internal photon name appears twice, 'mu' = external index
    ('s', j)     numerator slot j: Q_j-slash + cp_j p-slash + cpp_j p'-slash + 1   (F2)
                                   Q_j-slash + cP_j p-slash + 1                    (V0, SE: p = p')
The Q_j are Wick-paired over all perfect matchings, each pair giving -1/2 w(i,j) g^{ab}.

Result: the node polynomial as {monomial: Fraction}, where a monomial is a sorted tuple of
factor names 'cp3', 'cpp0', 'cP2', 'w1_4' (slot indices, not slotline) -- exactly
sum over dirac.table entries of val * (-1/2)^w * prod(factors).  For V0 and SE each monomial
maps to a pair (a, a+b+c) or (u, v), as in dirac.py.
"""
import subprocess, tempfile, os, re
from fractions import Fraction

FORM = os.environ.get('FORM', 'form')

def _gamma_string(chain, mode, line=1):
    """mode F2: slots Q + cp p + cpp p' + 1;  V0/SE: Q + cP p + 1;  PI: Q + ck k + 1."""
    idx = set(); out = []
    for tok in chain:
        if tok[0] == 'v':
            nm = tok[1]
            if nm not in ('mu', 'nu'): idx.add('i' + nm)
            out.append('g_(%d,%s)' % (line, nm if nm in ('mu', 'nu') else 'i' + nm))
        else:
            j = tok[1]; q = 'g_(%d,q%d)*Q(%d,q%d)' % (line, j, j, j)
            if mode == 'F2':
                out.append('(%s+cp%d*g_(%d,p)+cpp%d*g_(%d,pp)+gi_(%d))' % (q, j, line, j, line, line))
            elif mode == 'PI':
                out.append('(%s+ck%d*g_(%d,k)+gi_(%d))' % (q, j, line, line))
            else:
                out.append('(%s+cP%d*g_(%d,p)+gi_(%d))' % (q, j, line, line))
    return '*'.join(out) if out else 'gi_(%d)' % line, idx

def program(chain, nslots, mode, loop=None):
    """loop: token list of a closed fermion loop on spin line 2; its trace multiplies every expression
    (the loop sign is not included, as in dirac.py).  mode PI: the loop alone, external indices mu, nu,
    momentum k; the result is the coefficient of k^mu k^nu, a polynomial in ksq = k.k."""
    if mode == 'PI': G, idx = 'gi_(1)', set()
    else: G, idx = _gamma_string(chain, mode)
    LP = ''
    if loop is not None:
        LP, idx2 = _gamma_string(loop, mode, 2); idx |= idx2; LP = '*' + LP
    syms = []
    for j in range(nslots):
        syms += ['cp%d' % j, 'cpp%d' % j] if mode == 'F2' else (['ck%d' % j] if mode == 'PI' else ['cP%d' % j])
    qs = ['q%d' % j for j in range(nslots)]
    L = ['#-', 'Off statistics;', 'Format nospaces;', 'Format 250;',
         'Symbols t(:1),i1,i2,ksq,KK%s;' % ''.join(',' + s for s in syms),
         'CFunction Q,Q1,w(symmetric);', 'Vectors p,pp,k;',
         'Indices mu,nu,n1,n2%s;' % ''.join(',' + s for s in sorted(idx) + qs)]
    if mode == 'F2':
        L.append('Local N = (4-t)*(g_(1,pp)+gi_(1))*%s*(g_(1,p)+gi_(1))*g_(1,mu)%s'
                 ' - (2+t)*(g_(1,pp)+gi_(1))*%s*(g_(1,p)+gi_(1))*(p(mu)+pp(mu))%s;' % (G, LP, G, LP))
    elif mode == 'V0':
        # Gamma^mu = a g^mu + b p^mu + c p^mu pslash + d [pslash, g^mu]; a+b+c = (T4+T2)/4, a = (T1-T4)/12
        L.append('Local T1 = %s*g_(1,mu)%s;' % (G, LP))
        L.append('Local T2 = %s*p(mu)%s;' % (G, LP))
        L.append('Local T4 = %s*g_(1,p)*p(mu)%s;' % (G, LP))
    elif mode == 'PI':
        L.append('Local T = %s;' % LP[1:])
    else:
        L.append('Local Tu = %s%s;' % (G, LP))
        L.append('Local Tv = %s*g_(1,p)%s;' % (G, LP))
    L += ['repeat;', '  id once Q(?a) = Q1(?a);',
          '  id all Q1(i1?,n1?)*Q(i2?,n2?) = -1/2*w(i1,i2)*d_(n1,n2);', 'endrepeat;',
          'id Q1(?a) = 0;'] + (['trace4,2;'] if loop is not None else []) + ['trace4,1;', 'id p.p = 1;']
    if mode == 'F2':
        L += ['id pp.pp = 1;', 'id p.pp = 1-t/2;', '.sort', 'if (count(t,1) != 1) discard;', 'id t = 1;']
    if mode == 'PI':
        L += ['id k(mu)*k(nu) = KK;', 'if (count(KK,1) != 1) discard;', 'id KK = 1;', 'id k.k = ksq;']
    L += ['.sort', 'Print +s;', '.end']
    return '\n'.join(L) + '\n'

_term = re.compile(r'^([+-])(?:(\d+)(?:/(\d+))?)?\*?(.*)$')

def _parse(text):
    """FORM 'Print +s' output -> {name: {monomial: Fraction}}."""
    res = {}; cur = None; buf = ''
    lines = [l.strip() for l in text.splitlines()]
    exprs = {}
    for l in lines:
        m = re.match(r'^(\w+)=$', l)
        if m: cur = m.group(1); exprs[cur] = []; continue
        if cur is None or not l: continue
        if l.endswith(';'): exprs[cur].append(l[:-1]); cur = None; continue
        exprs[cur].append(l)
    for nm, ts in exprs.items():
        d = {}
        # terms may wrap over lines: join, then split at top-level signs
        s = ''.join(ts)
        for t in re.findall(r'[+-][^+-]*(?:\([^)]*\)[^+-]*)*', s):
            m = _term.match(t)
            sign, num, den, rest = m.groups()
            c = Fraction(int(num) if num else 1, int(den) if den else 1)
            if sign == '-': c = -c
            fac = []
            for f in [x for x in rest.split('*') if x]:
                mm = re.match(r'^(\w+)(?:\((\d+),(\d+)\))?(?:\^(\d+))?$', f)
                assert mm, f
                base = mm.group(1) if mm.group(2) is None else 'w%s_%s' % tuple(sorted((int(mm.group(2)), int(mm.group(3)))))
                fac += [base]*(int(mm.group(4)) if mm.group(4) else 1)
            d[tuple(sorted(fac))] = d.get(tuple(sorted(fac)), 0) + c
        res[nm] = {k: v for k, v in d.items() if v != 0}
    return res

def run(prog, workers=1):
    with tempfile.TemporaryDirectory() as td:
        f = os.path.join(td, 'n.frm'); open(f, 'w').write(prog)
        exe = FORM if workers == 1 else 'tform'
        args = [exe, '-q'] + (['-w%d' % workers] if workers > 1 else []) + [f]
        r = subprocess.run(args, cwd=td, capture_output=True, text=True)
        if r.returncode: raise RuntimeError(r.stdout[-2000:] + r.stderr[-2000:])
        return r.stdout

def poly(chain, nslots, mode, loop=None, workers=1):
    out = _parse(run(program(chain, nslots, mode, loop), workers))
    if mode == 'F2': return out.get('N', {})
    if mode == 'PI': return out.get('T', {})
    if mode == 'SE':
        tu, tv = out.get('Tu', {}), out.get('Tv', {})
        keys = set(tu) | set(tv)
        return {k: (Fraction(tu.get(k, 0))/4, Fraction(tv.get(k, 0))/4) for k in keys}
    t1, t2, t4 = out.get('T1', {}), out.get('T2', {}), out.get('T4', {})
    # T1 = 16a + 4c, T2 = 4b, T4 = 4a + 4c  ->  a = (T1 - T4)/12, a+b+c = (T2 + T4)/4
    keys = set(t1) | set(t2) | set(t4)
    r = {k: ((t1.get(k, 0) - t4.get(k, 0))/12, (t2.get(k, 0) + t4.get(k, 0))/4) for k in keys}
    return {k: v for k, v in r.items() if v != (0, 0)}

def from_table(tab, mode):
    """dirac.table output -> same {monomial: value} form, for checking."""
    d = {}
    for (combo, pr), val in tab.items():
        fac = []
        for j, ch in enumerate(combo):
            if ch in ('p', 'pp', 'P'): fac.append('c%s%d' % (ch, j))
        for a, b in pr:
            i, j = sorted((a, b)); fac.append('w%d_%d' % (i, j))
        pw = Fraction(-1, 2)**len(pr)
        k = tuple(sorted(fac))
        v = val*pw if mode == 'F2' else tuple(Fraction(x)*pw for x in val)
        assert k not in d
        d[k] = v
    return d

def available():
    """the FORM executable, or None"""
    import shutil
    return shutil.which(FORM)

def table(chain, nslots, mode, loop=None):
    """drop-in for dirac.table (no closed loops): {(combo, pairing): value} with the Wick factor
    (-1/2)^w taken out again, so build.py treats it exactly like a dirac.py table."""
    P = poly(chain, nslots, mode, loop)
    if mode == 'PI':
        return _table_pi(P, nslots)
    tab = {}
    for mono, c in P.items():
        combo = [None]*nslots; pr = []
        for f in mono:
            m = re.match(r'c(pp|p|P)(\d+)$', f)
            if m: combo[int(m.group(2))] = m.group(1); continue
            a, b = map(int, f[1:].split('_')); combo[a] = combo[b] = 'Q'; pr.append((a, b))
        pw = Fraction(-1, 2)**len(pr)
        key = (tuple(x or '1' for x in combo), tuple(sorted(pr)))
        if mode == 'F2':
            v = c/pw; assert v.denominator == 1, (key, v); tab[key] = int(v)
        else:
            tab[key] = tuple(Fraction(x)/pw for x in c)
    return tab

def _table_pi(P, nslots):
    """PI polynomial (monomials may carry ksq) -> dirac.table_pi format {(combo, pairing): (C_0..C_M)}"""
    M = nslots//2 + 1
    tab = {}
    for mono, c in P.items():
        combo = ['1']*nslots; pr = []; m = 0
        for f in mono:
            if f == 'ksq': m += 1; continue
            mm = re.match(r'ck(\d+)$', f)
            if mm: combo[int(mm.group(1))] = 'k'; continue
            a, b = map(int, f[1:].split('_')); combo[a] = combo[b] = 'Q'; pr.append((a, b))
        key = (tuple(combo), tuple(sorted(pr)))
        v = tab.setdefault(key, [Fraction(0)]*(M + 1))
        v[m] += Fraction(c)/Fraction(-1, 2)**len(pr)
    return {k: tuple(v) for k, v in tab.items() if any(x != 0 for x in v)}

def table_pi(loop, nslots):
    return table(None, nslots, 'PI', loop)
