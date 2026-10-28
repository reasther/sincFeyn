"""Translate the double-precision pt_<graph> of a generated two-loop kernel (gen.py, FORM O3 output with
line identifications) into a Wolfram Language function  pt<name>[x0_, ..., x5_, pw_]  with the power of V
left open: pw[V, n] stands for 1/V^n.  The |.| companion sums (SQ_AB blocks, *a arrays) are dropped."""
import re, sys

def body_of(src, g):
    i = src.index('static REAL FN(pt_%s)' % g); i = src.index('{', i) + 1
    d = 1; j = i
    while d:
        d += {'{': 1, '}': -1}.get(src[j], 0); j += 1
    return src[i:j-1]

def items(b):
    """split into statements; a for(...) with its body, and an if(SQ_AB){...} block, are single items"""
    out = []; i = 0; n = len(b); cur = ''; ini = 0
    while i < n:
        if b.startswith('/*', i):
            k = b.index('*/', i); out.append(('comment', b[i+2:k].strip())); i = k + 2; continue
        if b.startswith('for(', i) or b.startswith('if(SQ_AB)', i):
            k = b.index('{', i); d = 1; j = k + 1
            while d:
                d += {'{': 1, '}': -1}.get(b[j], 0); j += 1
            out.append(('for' if b[i] == 'f' else 'skip', b[i:j])); i = j; continue
        c = b[i]
        if c == '{' and (cur.rstrip().endswith('=') or ini):
            ini += 1; cur += c
        elif c == '}' and ini:
            ini -= 1; cur += c
        elif c in '{}':
            if cur.strip(): out.append(('stmt', cur.strip()))
            cur = ''
        elif c == ';':
            if cur.strip(): out.append(('stmt', cur.strip()))
            cur = ''
        else:
            cur += c
        i += 1
    return out

def split_top(s):
    parts = []; d = 0; cur = ''
    for c in s:
        if c in '([{': d += 1
        if c in ')]}': d -= 1
        if c == ',' and d == 0: parts.append(cur); cur = ''
        else: cur += c
    parts.append(cur); return [p.strip() for p in parts]

NAMES = {}
def nm(c):
    if c not in NAMES:
        w = c.replace('_', '')
        if c.startswith('n_') and c.endswith('_W'): w = 'W' + c[2:-2].replace('_', '')
        if c.startswith('n_') and c.endswith('_V'): w = 'V' + c[2:-2].replace('_', '')
        assert w not in NAMES.values(), (c, w)
        NAMES[c] = w
    return NAMES[c]

def idx(m):
    return '[[' + ', '.join(str(int(k) + 1) for k in re.findall(r'\[(\d+)\]', m)) + ']]'

def expr(e):
    e = e.replace('(REAL)', '').replace('+Vx', '').replace('Vx+', '')
    e = re.sub(r'\b([A-Za-z_]\w*)((?:\[\d+\])+)', lambda m: nm(m.group(1)) + idx(m.group(2)), e)
    e = re.sub(r'\b([A-Za-z_]\w*)\b(?!\[\[)', lambda m: nm(m.group(1)) if m.group(1) not in ('pw',) else m.group(1), e)
    return e

def isabs(name):
    return name in ('ra', 'tota', 'aa') or name.endswith('_Wa') or re.fullmatch(r't\d+a', name) is not None

def translate(src, g, fname):
    out = []; locs = set()
    def assign(lhs, rhs, op='='):
        base = re.match(r'[A-Za-z_]\w*', lhs).group(0)
        if isabs(base) or base.startswith('a') and base[1:2] in ('c', 'w'): return
        locs.add(nm(base))
        L = expr(lhs)
        out.append('  %s %s %s;' % (L, {'=': '=', '+=': '+=', '-=': '-=', '*=': '*='}[op], expr(rhs)))
    for kind, s in items(body_of(src, g)):
        if kind == 'comment': out.append('  (* %s *)' % s); continue
        if kind == 'skip': continue
        if kind == 'for':
            hdr, bd = s[:s.index('{')], s[s.index('{')+1:-1]
            st = [x.strip() for x in bd.split(';') if x.strip()][0]
            m = re.fullmatch(r'(n_\w+_W)\[k\]\[p\]\[w\]\*=iU2', st)
            if m: out.append('  %s = %s iU2;' % (nm(m.group(1)), nm(m.group(1)))); continue
            m = re.fullmatch(r'(t\d+)\[w\]\+=(\w+)((?:\[\d+\])*)\[w\]', st)
            if m:
                r = nm(m.group(2)) + (idx(m.group(3)) if m.group(3) else '')
                out.append('  %s += %s;' % (nm(m.group(1)), r)); continue
            raise ValueError('for: ' + s)
        if re.match(r'const REAL x\d+=x\[\d+\]$', s): continue
        if s.startswith('return') or s.startswith('*ab'): continue
        decl = s.startswith('REAL ') or s.startswith('const REAL ')
        if decl:
            for p in split_top(s.split('REAL ', 1)[1]):
                m = re.fullmatch(r'(\w+)((?:\[\d+\])*)(?:=(.*))?', p, re.S)
                name, dims, init = m.group(1), m.group(2), m.group(3)
                if isabs(name): continue
                if name == 'ivp': out.append('  nn = 0;'); locs.add('nn'); continue
                if name == 'iV':
                    out.append('  vv = %s;' % expr(re.fullmatch(r'1/\((.*)\)', init).group(1))); locs.add('vv'); continue
                if dims:
                    ds = [int(k) for k in re.findall(r'\d+', dims)]
                    out.append('  %s = ConstantArray[0, {%s}];' % (nm(name), ', '.join(map(str, ds)))); locs.add(nm(name))
                elif init is not None: assign(name, init)
            continue
        if s == 'ivp*=iV': out.append('  nn = nn + 1; ivp = pw[vv, nn];'); locs.add('ivp'); continue
        m = re.fullmatch(r'(n_\w+)_U=U', s)
        if m: continue
        m = re.fullmatch(r'([\w\[\]]+)\s*(\+=|-=|\*=|=)\s*(.*)', s, re.S)
        if m: assign(m.group(1), m.group(3), m.group(2)); continue
        raise ValueError('stmt: ' + s)
    locs |= {'tot'}
    return '%s[x0_, x1_, x2_, x3_, x4_, x5_, pw_] := Module[{%s},\n%s\n  tot]' % (
        fname, ', '.join(sorted(locs)), '\n'.join(out))

if __name__ == '__main__':
    for path, g, f in zip(sys.argv[1::3], sys.argv[2::3], sys.argv[3::3]):
        NAMES.clear(); print(translate(open(path).read(), g, f)); print()
