#!/bin/sh
# Pointwise check of rewritten kernels against reference kernels of the same graphs, e.g. the
# one-point u rule (tools/glpts.py) against the 2-point Gauss-Legendre rule:
#   cp -r kernels /tmp/ref; python3 tools/glpts.py -n 2 /tmp/ref/*/*.c
#   sh tests/check_onept.sh /tmp/ref kernels [JOBS]
# For every kernel in NEWDIR that differs from REFDIR it evaluates f (double and the high-precision
# type) at 400 random nodes, sigma_i = e^z, |z| <= 3, and prints the largest difference relative to
# the reference |.| sum; where the kernel has the packed-point wrapper fv_, also the largest
# difference between fv_ and f in double (it should be 0).  Built with gcc as C (__float128).
REF=$1 NEW=$2 JOBS=${3:-4}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
W=${TMPDIR:-/tmp}/check_onept.$$; mkdir -p "$W"
one(){   # ref.c new.c
    d=$W/$(basename "$2" .c); mkdir -p "$d"
    name=$(grep -o '^const graph_t G_[A-Za-z0-9_]*' "$2" | sed 's/const graph_t G_//')
    fv=$(grep -c "^static void fv_$name(" "$2")
    sed "s/^const graph_t G_$name=/const graph_t Gref_$name=/" "$1" > "$d/ref.c"
    cat > "$d/main.c" <<EOT
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
__thread int SQ_AB = 1;
extern const graph_t G_$name, Gref_$name;
int main(void){
    const graph_t *g = &G_$name, *r = &Gref_$name; int D = g->D;
    double md = 0, mq = 0, mv = 0; srand(12345);
    for(int t = 0; t < 400; t++){
        double sf[7][16], ab, abr, out[7], abv[7]; QREAL sq[16], aq, aqr;
        for(int j = 0; j < 7; j++) for(int i = 0; i < D; i++) sf[j][i] = exp(6.0*rand()/RAND_MAX - 3);
        double a = g->fd(sf[0], &ab), b = r->fd(sf[0], &abr);
        for(int i = 0; i < D; i++) sq[i] = sf[0][i];
        QREAL qa = g->fq(sq, &aq), qb = r->fq(sq, &aqr);
        double e = fabs(a - b)/abr, eq = fabs((double)(qa - qb))/(double)aqr;
        if(e > md) md = e; if(eq > mq) mq = eq;
        if($fv){ g->fv(7, sf[0], 16, out, abv);
            for(int j = 0; j < 7; j++){ double aj, vj = g->fd(sf[j], &aj); e = fabs(out[j] - vj) + fabs(abv[j] - aj); if(e > mv) mv = e; } }
    }
    printf("%-10s D=%-2d  f %.1e  f(QREAL) %.1e", "$name", D, md, mq);
    if($fv) printf("  fv-f %.1e", mv);
    printf("%s\n", md < 1e-13 && mq < 1e-13 && mv == 0 ? "" : "  FAIL");
    return 0;
}
EOT
    CC="gcc -O1 -w -I$ROOT/driver"
    $CC -c -o "$d/ref.o" "$d/ref.c" && $CC -c -o "$d/new.o" "$2" && $CC -c -o "$d/main.o" "$d/main.c" &&
        gcc -o "$d/t" "$d/main.o" "$d/ref.o" "$d/new.o" -lquadmath -lm && "$d/t" || echo "$name: build failed  FAIL"
    rm -rf "$d"
}
for f in "$NEW"/*/*.c; do
    r=$REF/${f#$NEW/}; [ -f "$r" ] && ! cmp -s "$f" "$r" && echo "$r $f"
done | { n=0; while read r f; do one "$r" "$f" & n=$((n+1)); [ $((n % JOBS)) = 0 ] && wait; done; wait; } | tee "$W/out.txt"
bad=$(grep -c FAIL "$W/out.txt"); echo "$(wc -l < "$W/out.txt") kernels checked, $bad failed"; rm -rf "$W"
[ "$bad" = 0 ]
