"""PSLQ identification of the two-loop per-graph values (optional; needs mpmath, tested with 1.3.0).

    python3 pslq/identify.py pslq_runs/values.tsv [--basis small|large]

Reads values.tsv written by run_pslq.sh (graph, value, estimated error, source).  For each graph it runs
PSLQ on {value} + basis, with tolerance 10 x the estimated error, and prints the relation it finds as a
closed form, its residual against the Sinc value, and the comparison with the known closed forms.
The basis is that of A_1^{(4)}: {1, pi^2, pi^2 log 2, zeta(3)}; --basis large adds log 2, log^2 2 and pi^4.
It then checks the sum of the graphs against A_1^{(4)} and gives the crossed graph by subtraction.
"""
import sys
from mpmath import mp, mpf, pi, log, zeta, pslq, nstr

mp.dps = 50
L2, Z3 = log(2), zeta(3)
A4 = mpf(197)/144 + pi**2/12 - pi**2*L2/2 + 3*Z3/4
BASIS = {'small': [('1', mpf(1)), ('pi^2', pi**2), ('pi^2 log2', pi**2*L2), ('zeta3', Z3)]}
BASIS['large'] = BASIS['small'] + [('log2', L2), ('log2^2', L2**2), ('pi^4', pi**4)]
# the closed forms in exact_two_loop.md (Volkov's scheme, mirror graphs included)
KNOWN = {
    'ladder':      ('11/48 + pi^2/18', mpf(11)/48 + pi**2/18),
    'self-energy': ('-1/24 - pi^2/18', -mpf(1)/24 - pi**2/18),
    'corner':      ('-55/24 + pi^2/18 + pi^2 log2/3 - zeta3/2', -mpf(55)/24 + pi**2/18 + pi**2*L2/3 - Z3/2),
    'crossed':     ('1/6 + 13 pi^2/36 - 5 pi^2 log2/6 + 5 zeta3/4', mpf(1)/6 + 13*pi**2/36 - 5*pi**2*L2/6 + 5*Z3/4),
    'Pi2':         ('119/36 - pi^2/3', mpf(119)/36 - pi**2/3),
}
KNOWN['self-energy_demo'] = KNOWN['self-energy']

def frac(p, q):
    from fractions import Fraction
    return Fraction(p, q)

def form(rel, names):
    """value = -sum(c_i b_i)/c_0 as a readable string"""
    c0 = rel[0]; out = []
    for c, n in zip(rel[1:], names):
        if c == 0: continue
        f = frac(-c, c0)
        s = ('-' if f < 0 else '+') + ' '
        a = abs(f)
        num = '' if (a.numerator == 1 and n != '1') else str(a.numerator)
        body = (num + (' ' if num and n != '1' else '') + (n if n != '1' else '')) or '1'
        out.append(s + body + ('/%d' % a.denominator if a.denominator != 1 else ''))
    s = ' '.join(out)
    return s[2:] if s.startswith('+ ') else s

def main():
    args = sys.argv[1:]
    basis = 'small'
    if '--basis' in args:
        i = args.index('--basis'); basis = args[i+1]; del args[i:i+2]
    path = args[0] if args else 'pslq_runs/values.tsv'
    rows = [l.rstrip('\n').split('\t') for l in open(path) if l.strip() and not l.startswith('graph\t')]
    B = BASIS[basis]
    vals = {}
    print('basis:', ', '.join(n for n, _ in B))
    for name, v, err, src in rows:
        val = sum(mpf(x) for x in v.split(','))      # quad-double values come as their four exact parts
        err = mpf(err); vals[name] = (val, err)
        if name == 'sum': continue
        digits = int(-mp.log10(err / abs(val)))
        print('\n%-16s %s   (%s; est. error %s, ~%d digits)' % (name, nstr(val, max(digits + 3, 17)), src, nstr(err, 2), digits))
        rel = pslq([val] + [b for _, b in B], tol=10*err/abs(val), maxcoeff=10**6, maxsteps=10**6)
        if not rel or rel[0] == 0:
            print('  PSLQ: no relation with |coefficients| <= 1e6')
        else:
            cf = sum(-mpf(c)*b for c, (_, b) in zip(rel[1:], B)) / rel[0]
            ok = abs(val - cf) <= 10*err
            print('  PSLQ: %s    relation %s, residual %s%s' % (form(rel, [n for n, _ in B]), rel, nstr(val - cf, 3),
                  '' if ok else '   <- residual > 10 x error: spurious, need more digits'))
        if name in KNOWN:
            print('  known: %s = %s, Sinc - known = %s' % (KNOWN[name][0], nstr(KNOWN[name][1], 20), nstr(val - KNOWN[name][1], 3)))
    if 'sum' in vals:
        print('\nsum of the two-loop graphs (demo) - A_1^(4) = %s' % nstr(vals['sum'][0] - A4, 3))
    need = ('ladder', 'corner', 'self-energy', 'Pi2')
    if all(n in vals for n in need):
        best = {n: KNOWN[n][1] for n in need}
        cr = A4 - sum(best.values())
        lit = mpf('-0.467645446094')
        print('crossed = A_1^(4) - (closed forms of ladder, corner, self-energy, Pi2) = %s\n'
              '  vs the literature value -0.467645446094 (12 digits): diff %s' % (nstr(cr, 20), nstr(cr - lit, 3)))

if __name__ == '__main__':
    main()
