"""Split the FORM straight-line blocks of a kernel into chunk functions, in place.

    python3 tools/splitform.py [-n STMTS] [-t THRESHOLD] FILE...

A FORM-optimised node computes its Dirac polynomials in one block
    { REAL Z1_, ..., E0, ...;  Z1_=...; ...; E0=...; ...  n_..._W[...]+=E0; ... }
(also `if(SQ_AB){ ... }`, and the V4 copy made by vecg.py).  At five loops one block is ~165k
statements in a single function, and gcc needs more than 14 GB to compile it at any -O level.
Blocks with more than THRESHOLD statements (default 20000) are cut into functions of STMTS
statements (default 4000): the block's inputs (the node's local constants) are copied into an
array I_, the temporaries and outputs live in an array Q_, and each chunk is
    static void FN(<node>_cB_K)(const REAL *I, REAL *Q)
with the statements unchanged except for the renaming, so the arithmetic (and so the value) is
the same.  The chunks go just before the node function, under the same FN/FNV macro, so the
double, vector and high-precision copies each get their own.  Idempotent.
"""
import re, sys, argparse

ASSIGN = re.compile(r'^\s+(Z\d+_|E\d+)=(.*);$')
START = re.compile(r'^(\s*)(if\(SQ_AB\)\s*)?\{ (REAL|V4) (Z1_, .*);\s*$')
FUNC = re.compile(r'^static \w+ (FNV?)\((\w+)\)\(')
IDENT = re.compile(r'\b[A-Za-z_]\w*\b')
TEMP = re.compile(r'^(Z\d+_|E\d+)$')

def split_file(fn, nst, thr):
    lines = open(fn).read().split('\n')
    out = []; pending = []          # chunk functions waiting for the next function start to be emitted before it
    fstart = None; nblk = 0; changed = 0
    i = 0
    while i < len(lines):
        l = lines[i]
        m = FUNC.match(l)
        if m: fstart = (len(out), m.group(1), m.group(2)); nblk = 0
        ms = START.match(l)
        if not ms or fstart is None:
            out.append(l); i += 1; continue
        # the block: declaration lines, then assignments
        j = i + 1; decl = []
        while j < len(lines) and lines[j].lstrip().startswith('const '): decl.append(lines[j]); j += 1
        k = j; st = []
        while k < len(lines):
            a = ASSIGN.match(lines[k])
            if not a: break
            st.append((a.group(1), a.group(2))); k += 1
        if len(st) <= thr:
            out.append(l); i += 1; continue
        ind, typ = ms.group(1), ms.group(3)
        names = sorted({lhs for lhs, _ in st}, key=lambda s: (s[0], int(re.sub(r'\D', '', s))))
        qi = {nm: q for q, nm in enumerate(names)}
        ins = sorted({t for _, r in st for t in IDENT.findall(r) if not TEMP.match(t) and t not in qi})
        ii = {nm: q for q, nm in enumerate(ins)}
        def ren(t):
            t = t.group(0)
            if t in qi: return 'Q[%d]' % qi[t]
            if t in ii: return 'I[%d]' % ii[t]
            return t
        fm, base = fstart[1], fstart[2]
        cn = []
        for c in range(0, len(st), nst):
            nm = '%s_c%d_%d' % (base, nblk, c // nst); cn.append(nm)
            body = ['static void %s(%s)(const %s *I, %s *Q){' % (fm, nm, typ, typ)]
            body += ['  Q[%d]=%s;' % (qi[lhs], IDENT.sub(ren, r)) for lhs, r in st[c:c + nst]]
            body.append('}')
            pending.extend(body)
        nblk += 1; changed += 1
        outs = sorted({lhs for lhs, _ in st if lhs.startswith('E')}, key=lambda s: int(s[1:]))
        out.append(l)                                   # keeps the Z/E declarations; the Zs go unused
        out.extend(decl)
        out.append('%s  %s I_[%d]={%s}, Q_[%d];' % (ind, typ, max(len(ins), 1), ','.join(ins) or '0', len(names)))
        out.extend('%s  %s(%s)(I_, Q_);' % (ind, fm, nm) for nm in cn)
        out.extend('%s  %s=Q_[%d];' % (ind, e, qi[e]) for e in outs)
        # insert the chunks before the enclosing function
        p = fstart[0]
        out[p:p] = pending
        fstart = (p + len(pending), fstart[1], fstart[2]); pending = []
        i = k
    if changed: open(fn, 'w').write('\n'.join(out))
    return changed

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-n', type=int, default=4000); ap.add_argument('-t', type=int, default=20000)
    ap.add_argument('files', nargs='+')
    a = ap.parse_args()
    for f in a.files: print(f, split_file(f, a.n, a.t), 'blocks split')
