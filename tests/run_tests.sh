#!/bin/sh
# Quick regression runs against known values: sh tests/run_tests.sh [bin/sincqed]
# Each check runs one selection and compares the value with a reference within a tolerance set
# for the (coarse) step used here.  A few minutes on 4 cores.  Box |z| <= 40.
SQ=${1:-bin/sincqed}
OUT=${TMPDIR:-/tmp}/sincqed_tests.$$
mkdir -p "$OUT"
export SINCQED_LOG=- SINCQED_PROGRESS=0 SINCQED_CKPT=0
fail=0; n=0; skip=0
check(){   # name  expected  tolerance  run-args...
    name=$1 ref=$2 tol=$3; shift 3
    SINCQED_RESULTS="$OUT/$name.txt" "$SQ" run "$@" > "$OUT/$name.log" 2>&1
    if grep -q "graphs=0 " "$OUT/$name.log"; then   # its class is not in this build (make CLASSES=...)
        skip=$((skip+1)); printf "%-5s %-28s (no matching graph in this build)\n" SKIP "$name"; return; fi
    v=$(grep -v '^#' "$OUT/$name.txt" | awk '{for(i=2;i<=NF;i++) if($i ~ /^[+-][0-9]\.[0-9]+e[+-][0-9]+$/ && length($i) >= 20){s+=$i; break}} END{printf "%.12f", s}')
    ok=$(awk -v v="$v" -v r="$ref" -v t="$tol" 'BEGIN{d=v-r; if(d<0)d=-d; print (d<=t) ? "PASS" : "FAIL"}')
    n=$((n+1)); [ "$ok" = PASS ] || fail=$((fail+1))
    printf "%-5s %-28s %+.10f  ref %+.10f  diff %+.1e  tol %.0e\n" "$ok" "$name" "$v" "$ref" "$(awk -v v="$v" -v r="$ref" 'BEGIN{print v-r}')" "$tol"
}
# two-loop A2 = 197/144 + pi^2/12 - pi^2/2 ln2 + 3/4 zeta3 (all seven diagrams)
check two-loop-sum      -0.328478965579  3e-5  1.0 81 all2
# one-loop vertex with the Pi2 insertion: 119/36 - pi^2/3
check pi2-insertion      0.015687421859  3e-6  1.0 81 "VP: 1-loop vertex + Pi2"
# three-loop graph 16 (Volkov 3; 1-6, 2-5, 4-7), analytic to 6 digits
check graph16           -0.334695        5e-4  1.2 69 1e-8 "3; 1-6, 2-5, 4-7"
# three-loop light-by-light (6 orderings), analytic
check lbl3               0.3710052921    5e-4  1.5 55 1e-8 g3_lbl
sk=""; [ $skip -gt 0 ] && sk=", $skip skipped (classes not in this build)"
echo "$((n-fail))/$n passed$sk  (outputs in $OUT)"
[ $fail -eq 0 ] && [ $n -gt 0 ]
