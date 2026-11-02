#!/bin/sh
# convergence.sh -- (h, eps) convergence scan of the two- and three-loop graphs.
#
#   sh scripts/convergence.sh            run (or resume) the scan
#   sh scripts/convergence.sh summary    only rebuild the summary table from what has finished
#   sh scripts/convergence.sh plan       list the runs and which are done, run nothing
#
# One run per (graph, h, eps) on the plain grid (files keep a _plain suffix).  Each writes OUT/<graph>/h<h>_eps<eps>_<grid>.{txt,log,hist} and
# is skipped when its result exists, so the script can be killed and restarted at any time; an
# unfinished run resumes from its checkpoint (OUT/<graph>/ckpt_h<h>_eps<eps>_<grid>/).  Runs go coarse
# to fine: every h is done at every eps (and every grid) before the next, smaller h starts, so a
# partial scan is still a usable table.
#
# Settings (environment variables, defaults in brackets):
#   SQ       the executable                          [bin/sincqed]
#   OUT      output directory                        [conv]
#   GRAPHS   graph names or groups to scan            [all2 all3]  (groups are expanded to graphs;
#            -name drops one, e.g. "all3 -g3_lbl")
#   HS2      steps h for two-loop graphs              [1.5 1.3 1.2 1.0 0.9 0.8 0.7]
#   HS3      steps h for three-loop graphs            [1.5 1.3 1.2 1.0]
#   EPSS     pruning thresholds                       [1e-8 1e-9 1e-10 1e-11]
#   WINDOW   half-width of the summation box in z = ln s  [40]: N = 2*ceil(WINDOW/h)+1 for each h,
#            so h is the only discretisation parameter.
#            The hierarchical sum only visits the support, so a wide box costs little.
#            (|z| <= 24 truncates g2_02 and g2_03 at 1e-5 on the plain grid; 32: 4e-7; 40: 8e-9,
#            for 1.6x the evaluations of 24)
#   N        fixed lattice points per axis instead (33, 49, ...)  [unset]
#   THREADS  OpenMP threads                           [all cores]
#   MAXTIME  skip a (graph, h, eps) whose estimated run time exceeds this many seconds [0 = never]
#            (estimate: the same graph's last finished run, times 4 per eps decade and 2^D per halving of h)
#
# Output: OUT/totals.tsv, the sum over all two-loop (three-loop) graphs at each (h, eps, grid) against
# the exact A2 = -0.328478965579 (A3 = 1.181241456587, Laporta-Remiddi), once every graph is done;
# OUT/sets.tsv, the gauge-invariant three-loop sets with known values (Volkov's Table 4),
# OUT/compare.tsv, the grids side by side (value, diff from the reference, evals per grid),
# and OUT/summary.tsv, one row per finished run, tab-separated:
#   graph grid L D w h eps value err disc prune |S2h-Sh| evals time_s ref diff kept_terms terms_below_eps
# ref is the known value of the graph where one exists (scripts/refs.tsv: Volkov's 2018 MC values for
# all no-loop graphs at 1e-5 and exact single-graph sets; else the built-in table in driver/graphs.c: the
# graphs finite without subtraction, graph 16, the VP insertions and LbL), "-" otherwise.  Values of
# single graphs with subtractions depend on the subtraction scheme, so for most graphs there is no
# exact value to compare with; the totals (and gauge-invariant sets) are scheme independent.
# where kept_terms is the number of terms in the final (step h) sum and terms_below_eps the number
# of those smaller than eps*sum|terms| (from the histogram; a large count means the support was
# followed well past the threshold, a small one that the sum stops close to it).
# The .hist files hold, per refinement level, the count and the signed and absolute sums of the
# terms in half-decade bins of |term| (see SINCQED_HIST in README.md).

SQ=${SQ:-bin/sincqed}
OUT=${OUT:-conv}
GRAPHS=${GRAPHS:-"all2 all3"}
HS2=${HS2:-"1.5 1.3 1.2 1.0 0.9 0.8 0.7"}
HS3=${HS3:-"1.5 1.3 1.2 1.0"}
EPSS=${EPSS:-"1e-8 1e-9 1e-10 1e-11"}
WINDOW=${WINDOW:-40}
nfor(){ if [ -n "$N" ]; then echo "$N"; else awk -v h="$1" -v w="$WINDOW" 'BEGIN{ k = int(w/h); if(k*h < w - 1e-9) k++; print 2*k+1 }'; fi; }
GRIDS=plain
MAXTIME=${MAXTIME:-0}
MODE=${1:-run}

[ -x "$SQ" ] || { echo "convergence.sh: no executable $SQ (build with make, or set SQ=...)"; exit 1; }
case "$SQ" in /*) ;; *) SQ="$(pwd)/$SQ" ;; esac
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
[ -n "$THREADS" ] && export OMP_NUM_THREADS=$THREADS
export SINCQED_PROGRESS=${SINCQED_PROGRESS:-1800}

# graph list: "name L D" for every graph matched by GRAPHS (w3_* are internal checks, not graphs)
LIST="$OUT/.graphs"
"$SQ" list 2>/dev/null | awk -v sel="$GRAPHS" '
  BEGIN{ n = split(sel, s, " ") }
  /^  [a-z]/ { name = $1; L = substr($2, 3); D = ""; for(i = 3; i <= NF; i++) if($i ~ /^D=/) D = substr($i, 3);
    if(name ~ /^w3_/) next
    for(j = 1; j <= n; j++) if("-" name == s[j]) next          # "-graph" drops a graph, e.g. GRAPHS="all3 -g3_lbl"
    for(j = 1; j <= n; j++) if(s[j] == name || (s[j] == "all2" && L == 2) || (s[j] == "all3" && L == 3)) { print name, L, D; break } }' > "$LIST"
[ -s "$LIST" ] || { echo "convergence.sh: no graphs match '$GRAPHS'"; exit 1; }

# all (h, eps, graph, grid) points, coarse h first, then loose eps first
points(){
    for h in $(printf "%s\n" $HS2 $HS3 | sort -gru); do
        for e in $EPSS; do
            while read g L D; do
                if [ "$L" = 2 ]; then hs=$HS2; else hs=$HS3; fi
                case " $hs " in *" $h "*) for x in $GRIDS; do echo "$g $L $D $h $e $x"; done ;; esac
            done < "$LIST"
        done
    done
}

done_run(){ grep -q "^$1 " "$OUT/$1/h$2_eps$3_$4.txt" 2>/dev/null; }

# environment for a grid
gridenv(){
    case "$1" in
    plain)  echo "" ;;
    *) echo "convergence.sh: unknown grid '$1' (plain)" >&2; exit 1 ;;
    esac
}

# rough time estimate from the most recent finished run of the same graph and grid
estimate(){   # g D h e grid
    last=$(ls -t "$OUT/$1"/h*_eps*_$5.txt 2>/dev/null | while read f; do grep -q "^$1 " "$f" && { echo "$f"; break; }; done)
    [ -n "$last" ] || { echo 0; return; }
    awk -v g="$1" -v D="$2" -v h="$3" -v e="$4" '$1 == g { t = $(NF-2); h0 = $(NF-14); e0 = $(NF-12) }
        END { if(!t) { print 0; exit } r = t * (h0/h)^D * 4^((log(e0) - log(e))/log(10)); printf "%.0f\n", r }' "$last"
}

summary(){
    S="$OUT/summary.tsv"
    printf 'graph\tgrid\tL\tD\tw\th\teps\tvalue\terr\tdisc\tprune\tdS\tevals\ttime_s\tref\tdiff\tkept_terms\tterms_below_eps\n' > "$S"
    while read g L D; do
        for f in "$OUT/$g"/h*_eps*_*.txt; do
            [ -f "$f" ] && grep -q "^$g " "$f" || continue
            x=${f%.txt}; x=${x##*_}
            xs=$x
            hf="${f%.txt}.hist_$g.txt"
            # histogram: final level (0) term count, and count below eps * sum|terms|
            hc=$(awk 'NR == FNR { if($1 == g) eps = $(NF-12); next }
                      $1 == "0" { n += $4; a += $6; c[$2] += $4 }
                      END { cut = log(eps * a)/log(10); for(b in c) if(b + 0.5 <= cut) m += c[b]; printf "%d\t%d", n, m }' g="$g" "$f" "$hf" 2>/dev/null)
            [ -n "$hc" ] || hc="-	-"
            rf=$(awk -v g="$g" '$1 == g { print $2; exit }' "$(dirname "$0")/refs.tsv" 2>/dev/null)
            for x in $xs; do
            awk -v hc="$hc" -v g="$g" -v x="$x" -v rf="$rf" 'BEGIN{OFS="\t"} $1 == g {
                # the volkov string has spaces, so count fields from the end of the line:
                # L D w h N eps top value err disc prune trunc round S_2h evals time_s ref diff
                v = $(NF-10); s2 = $(NF-4); d = s2 - v; if(d < 0) d = -d; ds = s2 == "-" ? "-" : sprintf("%.3e", d)   # "-": no 2h level (top = 0)
                r = $(NF-1); df = $NF
                if(rf != "" && (r == "-" || g ~ /^g3_/)){ r = sprintf("%.12e", rf); df = sprintf("%+.3e", v - rf) }   # scripts/refs.tsv
                print g, x, $(NF-17), $(NF-16), $(NF-15), $(NF-14), $(NF-12), v, $(NF-9), $(NF-8), $(NF-7), ds, $(NF-3), $(NF-2), r, df, hc }' "$f"
            done
        done
    done < "$LIST" >> "$S"
    { head -1 "$S"; tail -n +2 "$S" | sort -t'	' -k1,1 -k6,6gr -k7,7gr -k2,2; } > "$S.tmp" && mv "$S.tmp" "$S"
    # grids side by side: one row per (graph, h, eps), value / diff / evals for each grid in GRIDS order
    C="$OUT/compare.tsv"
    awk -F'\t' -v grids="$GRIDS" 'BEGIN{ OFS = "\t"; ng = split(grids, G, " "); printf "graph\th\teps"; for(i = 1; i <= ng; i++) printf "\tvalue_%s\tdiff_%s\tevals_%s", G[i], G[i], G[i]; print "" }
        NR > 1 { k = $1 OFS $6 OFS $7; if(!(k in seen)){ seen[k] = 1; ord[++n] = k } V[k, $2] = $8; Df[k, $2] = $16; E[k, $2] = $13 }
        END { for(j = 1; j <= n; j++){ k = ord[j]; line = k; for(i = 1; i <= ng; i++){ g = G[i]
                  line = line OFS ((k, g) in V ? V[k, g] OFS Df[k, g] OFS E[k, g] : "-" OFS "-" OFS "-") } print line } }' "$S" > "$C"
    # totals: sum over the graphs of each loop order, against the exact coefficient, for every
    # (h, eps, grid) at which all of that order's graphs have finished
    # (only when GRAPHS covers the whole order: a subset is not compared with A2 or A3)
    T="$OUT/totals.tsv"
    nall(){ "$SQ" list 2>/dev/null | awk -v L="L=$1" '/^  [a-z]/ && $2 == L && $1 !~ /^w3_/' | wc -l; }
    n2=$(awk '$2 == 2' "$LIST" | wc -l); [ "$n2" -eq "$(nall 2)" ] || n2=-1
    n3=$(awk '$2 == 3' "$LIST" | wc -l); [ "$n3" -eq "$(nall 3)" ] || n3=-1
    awk -F'\t' -v n2="$n2" -v n3="$n3" 'BEGIN{ OFS = "\t"
            ex[2] = -0.328478965579194; ex[3] = 1.18124145658720; ng[2] = n2; ng[3] = n3 }
        NR > 1 { k = $3 OFS $2 OFS $6 OFS $7; if(!(k in n)) ord[++m] = k; n[k]++; v[k] += $8; e[k] += $9; ev[k] += $13 }
        END { for(j = 1; j <= m; j++){ k = ord[j]; split(k, q, OFS); L = q[1]
                  if(n[k] == ng[L] && (L in ex)) printf "%s\t%d\t%+.10f\t%+.10f\t%+.3e\t%.2e\t%.0f\n", k, n[k], v[k], ex[L], v[k] - ex[L], e[k], ev[k] } }' "$S" | sort -t'	' -k1,1n -k3,3gr -k4,4gr -k2,2 > "$T.tmp"
    { printf 'L\tgrid\th\teps\tgraphs\tsum\texact\tdiff\terr_sum\tevals\n'; cat "$T.tmp"; } > "$T"; rm -f "$T.tmp"
    # gauge-invariant sets with known values (Volkov, arXiv:1507.06435 Table 4, Fig. 3 numbering; the
    # analytic values are from the papers cited there): the sum of the members at each (h, eps, grid)
    # at which all of them have finished.  Our graphs include the mirror graph (w = 2), so one of ours
    # is a pair of his, e.g. g3_15 = his 39-40.
    V="$OUT/sets.tsv"
    awk -F'\t' 'BEGIN{ OFS = "\t"; ns = 0
        # name; members; exact
        d["V39-40"]    = "g3_15;-0.334695103723"
        d["V23-24"]    = "g3_23;1.861907872591"
        d["V25"]       = "g3_22;-0.026799490"
        d["V26-27"]    = "g3_24 g3_26;-3.17668477"
        d["V28"]       = "g3_27;1.79027778"
        d["V29-30"]    = "g3_12;-1.757936342"
        d["V31-32,37-38"] = "g3_10 g3_16;5.35763265"
        d["V33-34,37-38"] = "g3_13 g3_16;0.45545185"
        d["V31-32,35-36"] = "g3_10 g3_11;1.541649"
        d["V33-36"]    = "g3_13 g3_11;-3.360532"
        d["V41-48"]    = "g3_14 g3_17 g3_18 g3_19;-0.4029"
        d["V49-72"]    = "g3_00 g3_01 g3_02 g3_03 g3_04 g3_05 g3_06 g3_07 g3_08 g3_09 g3_20 g3_21 g3_25;0.9541"
        # all 28 no-lepton-loop graphs (Volkov arXiv:1807.05281: 0.90485(10), 5.3 days on one K80 GPU);
        # exact value = A3 minus the lepton-loop sets below (LbL, Pi4, Pi2Pi2, Pi2 in 4th-order vertices),
        # good to ~1e-5 because the 7-10 set is known to 4 digits
        d["V2018_noloop"] = "g3_00 g3_01 g3_02 g3_03 g3_04 g3_05 g3_06 g3_07 g3_08 g3_09 g3_10 g3_11 g3_12 g3_13 g3_14 g3_15 g3_16 g3_17 g3_18 g3_19 g3_20 g3_21 g3_22 g3_23 g3_24 g3_25 g3_26 g3_27;0.904978"
        d["V1-6_LbL"]  = "g3_lbl;0.3710052921"
        d["V7-10"]     = "v3_05 v3_06 v3_07;0.05015"
        d["V11-12,15-16"] = "v3_02 v3_04;-0.0879847"
        d["V13-14,17-18"] = "v3_01 v3_03;-0.112336"
        d["V19-21_Pi4"] = "v3_p4Sc v3_p4Sm v3_p4Vc v3_p4Vm;0.0528706524757670"
        d["V22_Pi2Pi2"] = "v3_00;0.00255852493651735"
        for(k in d){ split(d[k], p, ";"); nm[k] = split(p[1], M, " "); for(i = 1; i <= nm[k]; i++) mem[M[i]] = mem[M[i]] " " k; ex[k] = p[2] } }
        NR > 1 && ($1 in mem) { n2 = split(substr(mem[$1], 2), S, " ")
            for(i = 1; i <= n2; i++){ k = S[i] OFS $2 OFS $6 OFS $7; if(!(k in c)) ord[++m] = k; c[k]++; v[k] += $8; e[k] += $9; ev[k] += $13 } }
        END { for(j = 1; j <= m; j++){ k = ord[j]; split(k, q, OFS)
                  if(c[k] == nm[q[1]]) printf "%s\t%s\t%+.10f\t%+.10f\t%+.3e\t%.2e\t%.0f\n", k, d[q[1]] ~ /;/ ? substr(d[q[1]], 1, index(d[q[1]], ";") - 1) : "", v[k], ex[q[1]], v[k] - ex[q[1]], e[k], ev[k] } }' "$S" \
        | sort -t'	' -k1,1 -k3,3gr -k4,4gr -k2,2 > "$V.tmp"
    { printf 'set\tgrid\th\teps\tgraphs\tsum\texact\tdiff\terr_sum\tevals\n'; cat "$V.tmp"; } > "$V"; rm -f "$V.tmp"
    echo "summary: $S ($(($(wc -l < "$S") - 1)) runs); grids side by side: $C; totals per loop order: $T; known sets: $V"
}

case "$MODE" in
summary) summary; exit 0 ;;
plan)
    points | while read g L D h e x; do
        if done_run "$g" "$h" "$e" "$x"; then st=done; else st="todo (est. $(estimate "$g" "$D" "$h" "$e" "$x") s)"; fi
        printf '%-8s h=%-4s eps=%-6s %-7s %s\n' "$g" "$h" "$e" "$x" "$st"
    done; exit 0 ;;
run) ;;
*) echo "usage: sh scripts/convergence.sh [run|summary|plan]"; exit 1 ;;
esac

for x in $GRIDS; do gridenv "$x" > /dev/null || exit 1; done
echo "# convergence scan $(date '+%Y-%m-%d %H:%M:%S')  SQ=$SQ  N=${N:-auto (window $WINDOW)}  grids: $GRIDS  graphs: $(awk '{printf "%s ", $1}' "$LIST")" >> "$OUT/scan.log"
points | while read g L D h e x; do
    done_run "$g" "$h" "$e" "$x" && continue
    if [ "$MAXTIME" -gt 0 ]; then
        t=$(estimate "$g" "$D" "$h" "$e" "$x")
        if [ "$t" -gt "$MAXTIME" ]; then echo "skip $g h=$h eps=$e $x (est. ${t}s > MAXTIME)" | tee -a "$OUT/scan.log"; continue; fi
    fi
    d="$OUT/$g"; b="$d/h${h}_eps${e}_$x"; c="$d/ckpt_h${h}_eps${e}_$x"; mkdir -p "$d" "$c"
    echo "$(date '+%H:%M:%S') start $g h=$h eps=$e $x" | tee -a "$OUT/scan.log"
    env $(gridenv "$x") SINCQED_RESULTS="$b.txt" SINCQED_LOG="$b.log" SINCQED_HIST="$b.hist" SINCQED_CKPT_DIR="$c" \
        "$SQ" run "$h" "$(nfor "$h")" "$e" "$g" > /dev/null 2>&1 < /dev/null
    if done_run "$g" "$h" "$e" "$x"; then
        rmdir "$c" 2>/dev/null
        echo "$(date '+%H:%M:%S') done  $(grep "^$g " "$b.txt" | awk -v x="$x" '{printf "%s h=%s eps=%s %s value %s err %s evals %s time %ss", $1, $(NF-14), $(NF-12), x, $(NF-10), $(NF-9), $(NF-3), $(NF-2)}')" | tee -a "$OUT/scan.log"
    else
        echo "$(date '+%H:%M:%S') FAILED $g h=$h eps=$e $x (see $b.log)" | tee -a "$OUT/scan.log"
    fi
    summary > /dev/null     # tables stay current while the scan runs
done
summary
