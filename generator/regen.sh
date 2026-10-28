#!/bin/sh
# Regenerate kernel sources from the Python generator and sort them into OUTDIR/<class>/.
#
#   sh generator/regen.sh OUTDIR [JOBS] [CLASS ...]
#
# CLASS is any of  g2 g3 v3 v4 w4 pi6 lbl4 g4  (default: g2 g3 v3 v4 w4 pi6, everything that is
# committed to the repository).  lbl4 and g4 are long: the four-loop light-by-light Dirac tables
# take tens of core-hours, the 819 Dirac tables of g4 about 30 core-hours plus emission.
#
# The pipeline is
#   1. generate:   gen.py, vp.py, vp4.py, vp4x.py, pi6.py, lbl*.py, gen4all.sh  (node code inline,
#                  GEN_SPLIT=0, except pi6.py, whose kernels use one C function per node)
#   2. vectorise:  tools/vec4.py on the 4-node three-loop kernels, tools/vecg.py on the
#                  helper-function layout (g3_07, g3_25, every g4 file)
#   3. optimise:   tools/optimize.sh: shared sub-sums as straight-line code (compact.py) and the
#                  exact one-point rule in u (glpts.py), and the packed-point wrapper fv_ in the
#                  vectorised kernels (packpts.py).  Values are unchanged to rounding.
# Its output is byte-identical to kernels/ (tests/check_regen.sh compares the checksums in
# kernels/SHA256SUMS).  REGEN_PLAIN=1 stops after step 2 and leaves the plain factored kernels
# (slower to compile and to run; useful for reading the generated code or debugging a generator).
#
# Dirac algebra: gen.py and vp.py (g2 g3 g4, v3 v4) use FORM (https://github.com/form-dev/form) when
# the `form` executable is on the PATH (or $FORM): exact tables in seconds, node polynomials as FORM
# Format O3 code.  Without FORM they warn and fall back to Python (the output then differs from
# kernels/SHA256SUMS, which were made with FORM 4.3).  The four-loop light-by-light generators (lbl4b,
# lbl4c, lbl4d, lblvp) also use FORM when it is found, optimising each node's Wick blocks jointly
# (l4_c13 1.9x faster than the Python trees); the other lepton-loop generators (lbl, pi6, vp4, vp4x)
# emit Python factored trees, where compact.py's sub-sums shared across Wick blocks do as well.
# REGEN_G4_JOINT=1 optimises the g4 kernels jointly too: 1.27x on g4_02, but 1.02x over all 269 (0.73x
# to 1.36x per graph, single-core x86), so the released g4 kernels stay per block.
# Every generator writes the node polynomials in one symbol per distinct coefficient (line
# identifications, build.py _monomials; GEN_IDENT=0 restores one symbol per slot).
#
# Dirac tables are cached next to the generator (generator/cache_dirac.pkl, cache_lbl.pkl and
# OUTDIR/.raw/.tables, removed when g4 is done); delete the .pkl files to recompute everything from scratch.  Requires python3 with
# numpy and sympy.
set -e
OUT=${1:?usage: sh generator/regen.sh OUTDIR [JOBS] [CLASS ...]}; shift
JOBS=${1:-4}; [ $# -gt 0 ] && shift
CLASSES=${*:-g2 g3 v3 v4 w4 pi6}
HERE=$(cd "$(dirname "$0")" && pwd); TOOLS=$HERE/../tools
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd); RAW=$OUT/.raw; mkdir -p "$RAW"
export GEN_SPLIT=0 PYTHONHASHSEED=0
has(){ case " $CLASSES " in *" $1 "*) return 0;; esac; return 1; }
cd "$HERE"
log(){ echo "== $* ($(date))"; }

# two- and three-loop graphs without lepton loops (gen.py 3 also fills the cache vp*.py read)
if has g2 || has g3 || has v3 || has v4 || has w4 || has pi6; then
    log "gen.py 2, 3"; python3 gen.py 2 "$RAW" --jobs "$JOBS"; python3 gen.py 3 "$RAW" --jobs "$JOBS"
fi
if has g2 || has v3; then log "vp.py (Pi2 insertions, two and three loops)"; python3 vp.py "$RAW"; fi
if has v3; then log "vp4.py (Pi4 in the one-loop vertex)"; python3 vp4.py "$RAW"
               log "pi6.py pi4 (Pi4 check via the general loop code)"; GEN_SPLIT=1 python3 pi6.py "$RAW" pi4; fi
if has g3; then log "lbl.py (three-loop light-by-light)"; python3 lbl.py "$RAW" "$JOBS"; fi
if has v4; then log "vp.py 4 (Pi2 insertions, four loops)"; python3 vp.py "$RAW" 4; fi
if has w4; then log "vp4x.py (Pi4 insertions, four loops)"; python3 vp4x.py "$RAW"; fi
if has pi6; then log "pi6.py (proper Pi6)"; GEN_SPLIT=1 python3 pi6.py "$RAW"; fi
if has lbl4; then
    log "lbl4b.py, lbl4c.py, lbl4d.py, lblvp.py (four-loop light-by-light)"
    # GEN_CHUNK=0: unchunked, as in the kernels used for the four-loop runs
    export GEN_CHUNK=0
    python3 lbl4b.py "$RAW" "$JOBS"; python3 lbl4c.py "$RAW" "$JOBS"; python3 lbl4d.py "$RAW" "$JOBS"; python3 lblvp.py "$RAW"
    unset GEN_CHUNK
fi
if has g4; then log "gen4all.sh (four loops, no lepton loops)"; GEN_FORM=${REGEN_G4_FORM:-1} GEN_FORM_JOINT=${REGEN_G4_JOINT:-0} GEN_CHUNK=0 sh gen4all.sh "$RAW" "$JOBS" && rm -rf "$RAW/.tables"; fi

# sort into classes
class_of(){
    case $1 in
        g2_*|v2_*) echo g2;; g3_*) echo g3;; v3_*|w3_*) echo v3;; g4_*) echo g4;; v4_*) echo v4;;
        w4_p6*) echo pi6;; w4_*) echo w4;; l4_*) echo lbl4;; *) echo "";;
    esac
}
for f in "$RAW"/*.c; do
    [ -e "$f" ] || continue
    b=$(basename "$f"); c=$(class_of "$b")
    has "$c" || continue
    mkdir -p "$OUT/$c"; cp "$f" "$OUT/$c/$b"
done

# vectorise the Gauss-Legendre nodes of the double half (deterministic source rewrites)
log "vectorise"
if [ -d "$OUT/g3" ]; then
    for f in "$OUT"/g3/g3_*.c; do python3 "$TOOLS/vec4.py" "$f" > "$f.tmp" && mv "$f.tmp" "$f"; done
    python3 "$TOOLS/vecg.py" -i "$OUT/g3/g3_07.c" "$OUT/g3/g3_25.c"
fi
if [ -d "$OUT/g4" ]; then ls "$OUT"/g4/g4_*.c | xargs -n 4 -P "$JOBS" python3 "$TOOLS/vecg.py" -i; fi
if [ "${REGEN_PLAIN:-0}" = 0 ]; then log "optimise"; OPT_JOBS=$JOBS sh "$TOOLS/optimize.sh" "$OUT"; fi
log "done: $(ls "$OUT"/*/*.c | wc -l) kernel files in $OUT"
