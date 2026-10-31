"""Look up a five-loop graph without lepton loops: kernel name, D, subtraction terms, and Volkov's
2019 value when his supplemental tables are given.

    python3 scripts/g5info.py "6; 1-8, 2-9, 3-10, 4-11, 5-7" [SUPPDIR]
    python3 scripts/g5info.py g5_3086 [SUPPDIR]

The graph may be given in either orientation (Volkov's code is the lexicographically smaller of a
graph and its mirror; ours is the representative enumerate_graphs keeps).  SUPPDIR holds Volkov's
class_*_aver.txt (supplement of PRD 100, 096004): tab-separated, comma decimals, values that include
the mirror graph, as ours do.  Prints one line of tab-separated fields:
    name  volkov  D  w  terms  ref  sigma_up  sigma_down     (ref fields "-" without SUPPDIR)
"""
import sys, os, re, glob
os.environ.setdefault('_GEN_HBASIS_NOTED', '1')      # a lookup, not a generator run: no GEN_HBASIS note
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'generator'))
import graphs as gr, build

def parse(s):
    m = re.match(r'\s*(\d+)\s*;\s*(.*)$', s)
    if not m: return None
    arcs = [tuple(int(x) - 1 for x in re.split(r'[-*]', a.strip())) for a in re.split(r'[,;]', m.group(2)) if a.strip()]
    return gr.Graph(5, int(m.group(1)) - 1, arcs)

def code(g):    # Volkov's supplement notation
    return '%d;' % (g.ext + 1) + ';'.join('%d*%d' % (a + 1, b + 1) for a, b in g.arcs)

def main():
    arg = sys.argv[1]; supp = sys.argv[2] if len(sys.argv) > 2 else os.environ.get('VOLKOV_SUPP')
    reps = gr.enumerate_graphs(5)
    if re.match(r'g5_\d+$', arg): i = int(arg[3:]); g, w = reps[i]
    else:
        q = parse(arg)
        keys = {q.key(), q.mirror().key()}
        hit = [(i, gw) for i, gw in enumerate(reps) if gw[0].key() in keys]
        if not hit: sys.exit('not a 1PI five-loop graph without lepton loops: %s' % arg)
        i, (g, w) = hit[0]
    D = 1 + len(build.classes(g)[1])
    ref = ['-', '-', '-']
    if supp:
        want = {code(g), code(g.mirror())}
        for f in glob.glob(os.path.join(supp, 'class_*_aver.txt')):
            for l in open(f):
                c = l.rstrip('\n').split('\t')
                if len(c) > 4 and c[1] in want: ref = [x.replace(',', '.') for x in c[2:5]]
    print('\t'.join(['g5_%02d' % i, g.volkov(), str(D), str(w), str(len(gr.terms(g)))] + ref))

if __name__ == '__main__':
    main()
