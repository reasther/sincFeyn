#!/usr/bin/env python3
"""compact.py KERNEL.c > OUT.c, or compact.py -i KERNEL.c ... (in place).

Rewrites the Dirac polynomials of a generated kernel (the W/Wa accumulations of every node
function and the ch_ chunk functions they call) as straight-line code over shared
subexpressions: every sub-sum of the factored (Horner) tree is normalised (integer content
pulled out, terms sorted) and computed once per node, and the |.| companion sums go in a
separate block under if(SQ_AB), so points that skip them no longer pay for them.  Works on kernels
before or after vecg.py.  Values agree with the original to rounding."""
import os, re, sys
from fractions import Fraction
from math import gcd

sys.setrecursionlimit(1000000)
TOK = re.compile(r'\s*(\(\(REAL\)-?\d+/\d+\)|\d+|[A-Za-z_]\w*(?:\[\d+\])?|[()+\-*])')

def tokenize(s):
    out, pos = [], 0
    while pos < len(s):
        m = TOK.match(s, pos)
        if not m:
            if s[pos:].strip() == '': break
            raise ValueError('cannot parse near: ' + s[pos:pos+60])
        out.append(m.group(1)); pos = m.end()
    return out

class Parser:
    """expr := ['-'] term (('+'|'-') term)* ; term := atom ['*' atom] ; atom := num | name | '(' expr ')'"""
    def __init__(self, toks): self.t, self.i = toks, 0
    def peek(self): return self.t[self.i] if self.i < len(self.t) else None
    def nxt(self): self.i += 1; return self.t[self.i-1]
    def expr(self):
        terms = []; sign = 1
        if self.peek() == '-': self.nxt(); sign = -1
        terms.append(self.term(sign))
        while self.peek() in ('+', '-'):
            sign = 1 if self.nxt() == '+' else -1
            terms.append(self.term(sign))
        return ('S', terms)
    def term(self, sign):
        f = [self.atom()]
        while self.peek() == '*': self.nxt(); f.append(self.atom())
        return (sign, f)
    def atom(self):
        x = self.nxt()
        if x == '(':
            e = self.expr(); assert self.nxt() == ')'; return e
        if x.startswith('((REAL)'):
            n, d = x[7:-1].split('/'); return ('C', Fraction(int(n), int(d)))
        if x.isdigit(): return ('C', Fraction(int(x)))
        return ('V', x)

def parse(s):
    p = Parser(tokenize(s)); e = p.expr(); assert p.i == len(p.t), s[:200]; return e

def content(cs):
    """positive rational gcd of the coefficients"""
    n = 0; d = 1
    for c in cs: n = gcd(n, abs(c.numerator)); d = d*c.denominator//gcd(d, c.denominator)
    return Fraction(n, d)

class Dag:
    """Canonical sub-sums: key = sorted tuple of (coef, var or None, child id or None)."""
    def __init__(self, chunks):
        self.chunks = chunks; self.ids = {}; self.nodes = []; self.memo = {}; self.keep = []
    def norm(self, e):
        """e -> (c, v, i): e == c * v * node_i (v, i may be None)"""
        if e[0] == 'C': return (e[1], None, None)
        if e[0] == 'V':
            m = re.fullmatch(r'c(\d+)\[0\]', e[1])
            if m: return self.norm(self.chunks[int(m.group(1))])
            return (Fraction(1), e[1], None)
        key0 = id(e); self.keep.append(e)   # keep e alive so its id is not reused
        if key0 in self.memo: return self.memo[key0]
        acc = {}
        for sign, f in e[1]:
            c = Fraction(sign); vs = []; ch = []
            for a in f:
                cc, v, i = self.norm(a)
                c *= cc
                if v is not None: vs.append(v)
                if i is not None: ch.append(i)
            if c == 0: continue
            # a product with more than one factor besides the coefficient: var * child, as the generator emits
            if len(vs) + len(ch) > 2 or len(ch) > 1 or len(vs) > 1:
                raise ValueError('unexpected product shape %r %r' % (vs, ch))
            k = (vs[0] if vs else None, ch[0] if ch else None)
            acc[k] = acc.get(k, 0) + c
        terms = [(c, v, i) for (v, i), c in acc.items() if c != 0]
        if not terms: r = (Fraction(0), None, None)
        elif len(terms) == 1 and terms[0][1] is None: r = terms[0]
        else:
            g = content([c for c, _, _ in terms])
            terms.sort(key=lambda t: (t[1] or '', -1 if t[2] is None else t[2], t[0]))
            if terms[0][0] < 0: g = -g
            key = tuple((c/g, v, i) for c, v, i in terms)
            if key not in self.ids:
                self.ids[key] = len(self.nodes); self.nodes.append(key)
            r = (g, None, self.ids[key])
        self.memo[key0] = r
        return r

def cnum(c):
    if c.denominator == 1: return '%d' % c.numerator
    return '((REAL)%d/%d)' % (c.numerator, c.denominator)

def term_code(c, v, i, tv, A):
    """c * v * t_i as C (abs version when A): returns (sign, text)"""
    f = []
    if v is not None: f.append(('A' + v) if A else v)
    if i is not None: f.append('%s%d' % (tv, i))
    s = -1 if (c < 0 and not A) else 1
    ac = abs(c)
    if ac != 1 or not f: f.append(cnum(ac))
    return s, '*'.join(f)

BLOCK = int(os.environ.get('COMPACT_BLOCK', '1500'))   # terms per block function (0: one block)

def postorder(dag, roots):
    order = []; done = set()
    for _, (c, v, i) in roots:
        if i is None: continue
        stack = [(i, False)]
        while stack:
            k, exp = stack.pop()
            if k in done: continue
            if exp: done.add(k); order.append(k)
            else:
                stack.append((k, True))
                for c2, v2, i2 in reversed(dag.nodes[k]):
                    if i2 is not None and i2 not in done: stack.append((i2, False))
    return order

def emit_split(dag, roots, order, fname, A, ind):
    """Block functions FN(<fname><j>)(F, T) computing the nodes in post-order; nodes used outside
    their block go through T.  Returns (function texts, code for the caller, size of T)."""
    tv = 'a' if A else 'd'
    blk = {}; blocks = [[]]; n = 0
    for k in order:
        if BLOCK and n >= BLOCK and blocks[-1]: blocks.append([]); n = 0
        blocks[-1].append(k); blk[k] = len(blocks) - 1; n += len(dag.nodes[k])
    exp = {}
    for k in order:
        for c, v, i in dag.nodes[k]:
            if i is not None and blk[i] != blk[k] and i not in exp: exp[i] = len(exp)
    for _, (c, v, i) in roots:
        if i is not None and i not in exp: exp[i] = len(exp)
    vars_ = sorted({v for k in order for c, v, i in dag.nodes[k] if v is not None})
    vidx = {v: j for j, v in enumerate(vars_)}
    funcs = []
    def ref(i, b):
        return '%s%d' % (tv, i) if b is not None and blk[i] == b else 'T[%d]' % exp[i]
    def code_terms(terms, b):
        parts = []
        for c, v, i in terms:
            f = []
            if v is not None: f.append(('A' if A else '') + v)
            if i is not None: f.append(ref(i, b))
            ac = abs(c)
            if ac != 1 or not f: f.append(cnum(ac))
            neg = c < 0 and not A
            parts.append(('-' if neg else ('+' if parts else '')) + '*'.join(f))
        return ''.join(parts)
    for j, ks in enumerate(blocks):
        used = sorted({v for k in ks for c, v, i in dag.nodes[k] if v is not None}, key=vidx.get)
        o = ['static void FN(%s%d)(const REAL *F, REAL *T){' % (fname, j)]
        if used: o.append('  const REAL %s;' % ', '.join('%s%s=F[%d]' % ('A' if A else '', v, vidx[v]) for v in used))
        for k in ks:
            o.append('  const REAL %s%d=%s;' % (tv, k, code_terms(dag.nodes[k], j)))
            if k in exp: o.append('  T[%d]=%s%d;' % (exp[k], tv, k))
        o.append('}')
        funcs.append('\n'.join(o))
    F = 'FA_' if A else 'F_'; TT = 'TA_' if A else 'T_'
    call = []
    if vars_: call.append('%sconst REAL %s[%d]={%s};' % (ind, F, len(vars_), ','.join(('RABS(%s)' % v) if A else v for v in vars_)))
    else: call.append('%sconst REAL %s[1]={0};' % (ind, F))
    call.append('%sREAL %s[%d];' % (ind, TT, max(1, len(exp))))
    for j in range(len(blocks)):
        if blocks[j]: call.append('%sFN(%s%d)(%s, %s);' % (ind, fname, j, F, TT))
    for tgt, (c, v, i) in roots:
        if c == 0: continue
        t = code_terms([(c, v, None)], None) if i is None else code_terms([(c, v, i)], None)
        if i is not None: t = t.replace('T[', TT + '[')
        call.append('%s%s%s=%s;' % (ind, tgt, '-' if t.startswith('-') else '+', t.lstrip('-')))
    return funcs, call


XPOLY = int(os.environ.get('COMPACT_X', '1'))   # also factor the Schwinger-parameter polynomials (U, V, cp, w)
MONO = re.compile(r'([+-]?)\s*((?:\(?(?:x\d+|\d+)(?:\*x\d+)*\)?\*?)+)')

def monomials(s):
    """'x0*x1 - 2*(x3*x3)*x4 + ...' -> list of (Fraction, [vars]) or None when s is not such a sum"""
    s = s.replace(' ', '')
    if not re.fullmatch(r'[-+x0-9*()]+', s): return None
    out = []
    for m in re.finditer(r'([+-]?)([^+-]+)', s):
        c = Fraction(-1 if m.group(1) == '-' else 1); vs = []
        for f in m.group(2).replace('(', '').replace(')', '').split('*'):
            if re.fullmatch(r'x\d+', f): vs.append(f)
            elif re.fullmatch(r'\d+', f): c *= int(f)
            else: return None
        out.append((c, vs))
    return out

def horner(terms):
    """greedy most-common-factor tree, in the parser's ('S', [(sign, [atoms])]) form"""
    from collections import Counter
    const = sum((c for c, f in terms if not f), Fraction(0))
    rest = [(c, f) for c, f in terms if f]
    parts = []
    if const != 0: parts.append((1, [('C', const)]))
    while rest:
        cnt = Counter(x for _, f in rest for x in set(f))
        best = max(cnt.items(), key=lambda kv: (kv[1], kv[0]))[0]
        inn = []; out = []
        for c, f in rest:
            if best in f: g = list(f); g.remove(best); inn.append((c, g))
            else: out.append((c, f))
        parts.append((1, [('V', best), horner(inn)])); rest = out
    return ('S', parts if parts else [(1, [('C', Fraction(0))])])

def compact_x(keep, ind, tag=''):
    """factor the polynomial definitions (U, V, cp.., cpp.., w..) of a node with shared sub-sums"""
    defs = []
    for j, ln in enumerate(keep):
        m = re.match(r'(\s*\{ const REAL U=)(.*)(, iU=1/U;)$', ln)
        if m: defs.append((j, m.group(1), m.group(2), m.group(3))); continue
        m = re.match(r'(\s*const REAL \w+=\()(.*)(\)\*iU;)$', ln)
        if m: defs.append((j, m.group(1), m.group(2), m.group(3)))
    polys = []
    for j, a, e, b in defs:
        t = monomials(e)
        if t is None: return keep
        polys.append(horner(t))
    if not polys: return keep
    dag = Dag({}); roots = [dag.norm(p) for p in polys]
    order = postorder(dag, [(None, r) for r in roots])
    def code(terms):
        parts = []
        for c, v, i in terms:
            f = []
            if v is not None: f.append(v)
            if i is not None: f.append('e%s%d' % (tag, i))
            if abs(c) != 1 or not f: f.append(cnum(abs(c)))
            parts.append(('-' if c < 0 else ('+' if parts else '')) + '*'.join(f))
        return ''.join(parts)
    pre = ['%sconst REAL e%s%d=%s;' % (ind, tag, k, code(dag.nodes[k])) for k in order]
    keep = list(keep)
    for (j, a, e, b), r in zip(defs, roots):
        keep[j] = a + code([r]) + b
    # the U line opens a block: put the shared sub-sums before it
    j0 = defs[0][0]
    return keep[:j0] + pre + keep[j0:]

def compact_node(lines, name, chunkfn, tag=''):
    """lines of one node's code (a node function, or a node block inside pt_).
    Returns (new lines, block functions, sub-sums, terms), or None when there is nothing to do."""
    acc = []; keep = []; ins = None
    for ln in lines:
        s = ln.strip()
        if s.startswith('const REAL A') and 'RABS(' in s and ins is None: ins = len(keep); continue
        if s.startswith('const REAL F_['): continue
        if re.match(r'REAL c\d+\[2\]; FN\(ch_', s): continue
        ma = re.match(r'(n_\w+_W(?:\[\d+\])+)\+=(.*); (n_\w+_Wa(?:\[\d+\])+)\+=(.*);$', s)
        if ma: acc.append((ma.group(1), ma.group(2), ma.group(3))); continue
        keep.append(ln)
    if not acc or ins is None: return None
    chunks = {j: parse(e) for j, e in chunkfn.get('ch_' + name + '_', {}).items()}
    dag = Dag(chunks)
    roots = [(w, dag.norm(parse(e))) for w, e, _ in acc]
    aroots = [(wa, r) for (w, e, wa), (_, r) in zip(acc, roots)]
    ind = '    '
    order = postorder(dag, roots)
    f1, c1 = emit_split(dag, roots, order, 'bk_%s_' % name, False, ind)
    f2, c2 = emit_split(dag, aroots, order, 'ab_%s_' % name, True, ind + '  ')
    keep[ins:ins] = c1 + [ind + 'if(SQ_AB){'] + c2 + [ind + '}']
    if XPOLY: keep = compact_x(keep, '  ', tag)
    return keep, f1 + f2, len(dag.nodes), sum(len(k) for k in dag.nodes)

def compact_half(h):
    """h: text of one REAL half (double or QREAL)."""
    chunkfn = {}
    for m in re.finditer(r'static void FN\((ch_\w+?)(\d+)\)\(const REAL \*F, REAL \*r\)\{\n(.*?)\n\}\n', h, re.S):
        body = m.group(3)
        r0 = re.search(r'r\[0\]=(.*?); r\[1\]=', body, re.S).group(1)
        chunkfn.setdefault(m.group(1), {})[int(m.group(2))] = r0
    h = re.sub(r'static void FN\(ch_\w+\)\(const REAL \*F, REAL \*r\)\{\n.*?\n\}\n', '', h, flags=re.S)
    out = []; pos = 0; nnodes = 0; nterms = 0
    # node functions (split layout) and point functions (older layout: the nodes inline in pt_)
    for m in re.finditer(r'static (?:void FN\(nd_(\w+)\)|REAL FN\(pt_(\w+)\))\(.*?\n\}\n', h, re.S):
        lines = m.group(0).split('\n')
        if m.group(1):
            r = compact_node(lines, m.group(1), chunkfn)
            if r is None: continue
            new, funcs, n, t = r
        else:
            new = []; funcs = []; n = t = 0; i = 0
            while i < len(lines):
                mm = re.match(r'\s*/\* node ', lines[i])
                if not mm: new.append(lines[i]); i += 1; continue
                j = i
                while j < len(lines) and not re.search(r'n_\w+_U=U; n_\w+_V=V; \}$', lines[j]): j += 1
                blk = lines[i:j+1]
                nm = re.search(r'REAL (n_\w+?)_W\[', '\n'.join(blk))
                r = compact_node(blk, '%s_%s' % (m.group(2), nm.group(1)), chunkfn, '%d_' % len(funcs)) if nm else None
                if r is None: new += blk
                else: new += r[0]; funcs += r[1]; n += r[2]; t += r[3]
                i = j + 1
            if not funcs: continue
        nnodes += n; nterms += t
        out.append(h[pos:m.start()]); out.append('\n'.join(funcs) + '\n'); out.append('\n'.join(new)); pos = m.end()
    out.append(h[pos:])
    return ''.join(out), nnodes, nterms

def vector_copy(body):
    """vecg.py's rule: the vector (V4) copy of the scalar double helpers"""
    v = re.sub(r'\bREAL\b', 'V4', body).replace('RABS(', 'VABS(').replace('FN(', 'FNV(')
    num = r'[-+]?\(?[-+0-9.eE*/() ]*\)?'
    def wrap(m):
        if re.fullmatch(r'\s*', m.group(2)): return m.group(0)
        if re.search(r'\bint\s+$', v[max(0, m.start() - 8):m.start()]): return m.group(0)
        return '%s=VB(%s)' % (m.group(1), m.group(2))
    v = re.sub(r'(?<![\w\]\[<>=!])(\b[A-Za-z_]\w*(?:\[\d+\])*)=(' + num + r')(?=[,;])', wrap, v)
    return '#define FNV(x) x##_v\n' + v + '#undef FNV\n'

def transform(src):
    d0 = src.index('#define REAL double'); q0 = src.index('#define REAL QREAL')
    dh = src[d0:q0]
    if '#define FNV' in dh:
        # vectorised by vecg.py: compact the scalar helpers, then rebuild their vector copy
        t0 = dh.index('#if SQ_VW == 4\n#define NGV_')
        sc, tail = dh[:t0], dh[t0:]
        sc, n1, t1 = compact_half(sc)
        a = sc.index('#define FN(x) x##_d\n') + len('#define FN(x) x##_d\n')
        v0 = tail.index('#define FNV(x) x##_v\n'); v1 = tail.index('#undef FNV\n') + len('#undef FNV\n')
        dh = sc + tail[:v0] + vector_copy(sc[a:]) + tail[v1:]
    else:
        dh, n1, t1 = compact_half(dh)
    qh, _, _ = compact_half(src[q0:])
    sys.stderr.write('  %d shared sub-sums, %d terms\n' % (n1, t1))
    return src[:d0] + dh + qh

if __name__ == '__main__':
    a = sys.argv[1:]
    if a and a[0] == '-i':
        for f in a[1:]:
            s = open(f).read()
            if 'if(SQ_AB){' in s and '/* compact */' in s: continue
            sys.stderr.write(f + '\n'); t = transform(s)
            t = t.replace(' */\n', ' */\n/* compact */\n', 1)
            open(f, 'w').write(t)
    else:
        sys.stdout.write(transform(open(a[0]).read()))
