"""
LaTeX (TikZ) drawings of every diagram the code evaluates, labelled with the program's graph
names, Volkov strings and weights.  Mirror partners (weight 2) are drawn next to their
representative, so each page shows every Feynman diagram once: 7 at two loops, 72 at three.

    python3 diagrams.py OUT.tex       then  pdflatex OUT.tex
"""
import sys, itertools
import graphs as gr

HEAD = r"""\documentclass[10pt]{article}
\usepackage[a4paper,margin=15mm]{geometry}
\usepackage{tikz}
\usetikzlibrary{decorations.pathmorphing,decorations.markings}
\tikzset{
  photon/.style={decorate,decoration={snake,amplitude=0.55mm,segment length=1.8mm,pre length=0.3mm,post length=0.3mm}},
  fermion/.style={thick,postaction={decorate},decoration={markings,mark=at position 0.55 with {\arrow{stealth}}}},
  eline/.style={thick},
  vtx/.style={circle,fill,inner sep=0.8pt},
  bubble/.style={circle,draw,thick,fill=white,minimum size=5.6mm,inner sep=0pt},
}
\newcommand{\dcell}[2]{\begin{minipage}[t]{0.24\textwidth}\centering #1\\[0.5mm]{\scriptsize #2}\end{minipage}}
\setlength{\parindent}{0pt}
\begin{document}
\section*{Diagrams evaluated by \texttt{sincqed}}
Each diagram is labelled by the program entry that evaluates it, its Volkov string
(external-photon vertex; photon ends, counting vertices along the electron line from 1) and,
in brackets, the weight with which that entry already includes it.  Entries of weight 2
include the mirror diagram, drawn next to its representative.  Electron loops are drawn as
circles (one-loop vacuum polarization), a circle with an internal photon (two-loop vacuum
polarization), or a box (light-by-light).
"""

def bez(a, b, h, t):
    """point of the cubic Bezier (a,0) (a,h) (b,h) (b,0) at parameter t"""
    x = (1-t)**3*a + 3*(1-t)**2*t*a + 3*(1-t)*t**2*b + t**3*b
    y = 3*(1-t)**2*t*h + 3*(1-t)*t**2*h
    return x, y

def draw(L, ext, arcs, vp=None, sc=None):
    """TikZ for a vertex diagram.  arcs: 0-based positions; vp: {arc index: ['1','1'] or ['V'] or
    ['S'] or ['Sm']} electron-loop insertions on that photon."""
    vp = vp or {}
    sc = sc or (0.55 if L < 3 else 0.48)
    n = 2*L + 1
    X = lambda p: p + 1.0
    o = [r"\begin{tikzpicture}[scale=%.2f,baseline=(current bounding box.south)]" % sc]
    o.append(r"\draw[fermion] (0,0) -- (%g,0);" % X(0))
    for p in range(n - 1):
        o.append(r"\draw[fermion] (%g,0) -- (%g,0);" % (X(p), X(p+1)))
    o.append(r"\draw[fermion] (%g,0) -- (%g,0);" % (X(n-1), X(n-1) + 1))
    o.append(r"\draw[photon] (%g,0) -- (%g,-1.1);" % (X(ext), X(ext)))
    for k, (a, b) in enumerate(arcs):
        h = 0.55 + 0.62*(b - a)
        o.append(r"\draw[photon] (%g,0) .. controls (%g,%g) and (%g,%g) .. (%g,0);" % (X(a), X(a), h, X(b), h, X(b)))
        ins = vp.get(k, [])
        ts = [0.5] if len(ins) == 1 else [0.2, 0.8]
        for t, kind in zip(ts, ins):
            x, y = bez(X(a), X(b), h, t)
            o.append(r"\node[bubble] (B%d%s) at (%.3f,%.3f) {};" % (k, str(t).replace('.', ''), x, y))
            nm = "B%d%s" % (k, str(t).replace('.', ''))
            if kind == 'V':      # vertex-type two-loop VP: photon across the loop
                o.append(r"\draw[photon,decoration={amplitude=0.4mm,segment length=1.1mm}] (%s.north) -- (%s.south);" % (nm, nm))
            elif kind == 'S':    # self-energy-type: photon on the upper half
                o.append(r"\draw[photon,decoration={amplitude=0.4mm,segment length=1.1mm}] (%s.north west) to[bend left=70] (%s.north east);" % (nm, nm))
            elif kind == 'Sm':   # the same on the lower half
                o.append(r"\draw[photon,decoration={amplitude=0.4mm,segment length=1.1mm}] (%s.south west) to[bend right=70] (%s.south east);" % (nm, nm))
    for p in range(n):
        o.append(r"\node[vtx] at (%g,0) {};" % X(p))
    o.append(r"\end{tikzpicture}")
    return "\n".join(o)

def draw_lbl(perm, sc=0.55):
    """light-by-light: open line vertices o1 o2 o3; loop box L0 (external photon, top) L1 L2 L3.
    perm[a] = loop vertex of photon a+1."""
    o = [r"\begin{tikzpicture}[scale=%.2f,baseline=(current bounding box.south)]" % sc]
    xs = [1, 2.5, 4]
    o.append(r"\draw[fermion] (0,0) -- (1,0);\path (0,-1.1);")
    o.append(r"\draw[fermion] (1,0) -- (2.5,0);")
    o.append(r"\draw[fermion] (2.5,0) -- (4,0);")
    o.append(r"\draw[fermion] (4,0) -- (5,0);")
    # loop vertices: L0 top, L1 right, L2 bottom-ish middle, L3 left (a square rotated)
    pos = {0: (2.5, 3.7), 1: (3.5, 2.7), 2: (2.5, 1.7), 3: (1.5, 2.7)}
    for i in range(4):
        a, b = pos[i], pos[(i+1) % 4]
        o.append(r"\draw[fermion] (%g,%g) -- (%g,%g);" % (a + b))
    o.append(r"\draw[photon] (2.5,3.7) -- (2.5,4.7);")
    for a in range(3):
        lx, ly = pos[perm[a]]
        o.append(r"\draw[photon] (%g,0) -- (%g,%g);" % (xs[a], lx, ly))
    for x in xs: o.append(r"\node[vtx] at (%g,0) {};" % x)
    for i in range(4): o.append(r"\node[vtx] at (%g,%g) {};" % pos[i])
    o.append(r"\end{tikzpicture}")
    return "\n".join(o)

def tt(s): return r"\texttt{%s}" % s.replace('_', r'\_')

def cells(items):
    """items: list of (tikz, caption); 4 per row"""
    out = []
    for i in range(0, len(items), 4):
        row = [r"\dcell{%s}{%s}" % it for it in items[i:i+4]]
        out.append("\\hfill".join(row) + "\\hfill\\mbox{}\\par\\vspace{3mm}")
    return "\n".join(out)

def mirror_arc(L, arc):
    m = lambda p: 2*L - p
    return tuple(sorted((m(arc[0]), m(arc[1]))))

def section_novp(L):
    items = []
    for i, (g, w) in enumerate(gr.enumerate_graphs(L)):
        nm = 'g%d_%02d' % (L, i)
        items.append((draw(L, g.ext, g.arcs), r"%s\ \ %s\ \ [%d]" % (tt(nm), g.volkov(), w)))
        if w == 2:
            m = g.mirror()
            items.append((draw(L, m.ext, m.arcs), r"%s mirror: %s" % (tt(nm), m.volkov())))
    return items

def vp_items():
    """mirrors vp.py's enumeration (names v2_00, v3_00..07) and vp4.py (v3_p4*)"""
    two, three = [], []
    (g1, _), = gr.enumerate_graphs(1)
    two.append((draw(1, g1.ext, g1.arcs, {0: ['1']}), r"%s\ \ one-loop VP [1]" % tt('v2_00')))
    three.append((draw(1, g1.ext, g1.arcs, {0: ['1', '1']}), r"%s\ \ two one-loop VPs [1]" % tt('v3_00')))
    i = 1
    for g, w in gr.enumerate_graphs(2):
        done = set()
        for j in range(g.L):
            if j in done: continue
            wt, mj = w, None
            if w == 1:
                ma = mirror_arc(2, g.arcs[j]); mj = g.arcs.index(ma)
                if mj != j: wt = 2; done.add(mj)
                else: mj = None
            done.add(j)
            a, b = g.arcs[j]
            nm = 'v3_%02d' % i; i += 1
            three.append((draw(2, g.ext, g.arcs, {j: ['1']}), r"%s\ \ %s, VP on %d--%d [%d]" % (tt(nm), g.volkov(), a+1, b+1, wt)))
            if wt == 2:
                if mj is not None:          # self-mirror graph, mirror photon
                    three.append((draw(2, g.ext, g.arcs, {mj: ['1']}), r"%s mirror: VP on %d--%d" % (tt(nm), g.arcs[mj][0]+1, g.arcs[mj][1]+1)))
                else:                       # mirror graph
                    m = g.mirror(); mj2 = m.arcs.index(mirror_arc(2, g.arcs[j]))
                    three.append((draw(2, m.ext, m.arcs, {mj2: ['1']}), r"%s mirror: %s" % (tt(nm), m.volkov())))
    three.append((draw(1, g1.ext, g1.arcs, {0: ['V']}), r"%s + %s\ \ two-loop VP, vertex type [1]" % (tt('v3_p4Vm'), tt('v3_p4Vc'))))
    three.append((draw(1, g1.ext, g1.arcs, {0: ['S']}), r"%s + %s\ \ two-loop VP, self-energy type [2]" % (tt('v3_p4Sm'), tt('v3_p4Sc'))))
    three.append((draw(1, g1.ext, g1.arcs, {0: ['Sm']}), r"(included in %s, %s)" % (tt('v3_p4Sm'), tt('v3_p4Sc'))))
    return two, three

def main(out):
    two_vp, three_vp = vp_items()
    t = [HEAD]
    it2 = section_novp(2); it3 = section_novp(3)
    it_lbl = [(draw_lbl(p), r"%s\ \ ordering %s" % (tt('g3_lbl'), "".join(str(x) for x in p))) for p in itertools.permutations([1, 2, 3])]
    t.append(r"\subsection*{Two loops (%d diagrams)}" % (len(it2) + len(two_vp)))
    t.append(cells(it2 + two_vp))
    t.append(r"\subsection*{Three loops without lepton loops (%d diagrams, group \texttt{q3})}" % len(it3))
    t.append(cells(it3))
    t.append(r"\subsection*{Three loops with vacuum-polarization insertions (%d diagrams, group \texttt{vp3})}" % len(three_vp))
    t.append(cells(three_vp))
    t.append(r"\subsection*{Three-loop light-by-light (%d diagrams, entry \texttt{g3\_lbl}, all orderings summed pointwise)}" % len(it_lbl))
    t.append(r"Photons 1, 2, 3 attach left to right on the open line; the ordering gives the loop vertex (1 right, 2 bottom, 3 left of the box) of each.\par\medskip")
    t.append(cells(it_lbl))
    t.append(r"\end{document}")
    open(out, 'w').write("\n".join(t) + "\n")
    print('diagrams:', len(it2) + len(two_vp), 'two-loop,', len(it3) + len(three_vp) + len(it_lbl), 'three-loop ->', out)

if __name__ == '__main__':
    main(sys.argv[1])
