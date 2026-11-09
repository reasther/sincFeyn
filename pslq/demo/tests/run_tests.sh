#!/bin/sh
# Quick regression runs against known values: sh tests/run_tests.sh [bin/sincqed]
# Each check runs one selection and compares the value with a reference within a tolerance set
# for the (coarse) step used here.  A few minutes on 4 cores.  Box |z| <= 40.
SQ=${1:-bin/sincqed}
OUT=${TMPDIR:-/tmp}/sincqed_tests.$$
mkdir -p "$OUT"
export SINCQED_LOG=- SINCQED_PROGRESS=0 SINCQED_CKPT=0
fail=0; n=0
check(){   # name  expected  tolerance  run-args...
    name=$1 ref=$2 tol=$3; shift 3
    SINCQED_RESULTS="$OUT/$name.txt" "$SQ" run "$@" > "$OUT/$name.log" 2>&1
    v=$(grep -v '^#' "$OUT/$name.txt" | awk '{for(i=2;i<=NF;i++) if($i ~ /^[+-][0-9]\.[0-9]+e[+-][0-9]+$/ && length($i) >= 20){s+=$i; break}} END{printf "%.12f", s}')
    ok=$(awk -v v="$v" -v r="$ref" -v t="$tol" 'BEGIN{d=v-r; if(d<0)d=-d; print (d<=t) ? "PASS" : "FAIL"}')
    n=$((n+1)); [ "$ok" = PASS ] || fail=$((fail+1))
    printf "%-5s %-28s %+.10f  ref %+.10f  diff %+.1e  tol %.0e\n" "$ok" "$name" "$v" "$ref" "$(awk -v v="$v" -v r="$ref" 'BEGIN{print v-r}')" "$tol"
}
# two-loop A2 = 197/144 + pi^2/12 - pi^2/2 ln2 + 3/4 zeta3 (all seven diagrams)
check two-loop-sum      -0.328478965579  3e-5  1.0 81 all2
# one-loop vertex with the Pi2 insertion: 119/36 - pi^2/3
check pi2-insertion      0.015687421859  3e-6  1.0 81 "VP: 1-loop vertex + Pi2"
# demo: unreduced Sinc sums over all Schwinger parameters (one loop: 1/2; two loops: A2)
democheck(){   # name  L  expected  tolerance  demo-args...
    name=$1 L=$2 ref=$3 tol=$4; shift 4
    SINCQED_DEMO_RESULTS="$OUT/$name.txt" "$SQ" demo "$@" > "$OUT/$name.log" 2>&1
    v=$(awk -v L="L=$L" '$2 == "SUM" && $3 == L {print $6}' "$OUT/$name.txt" | tail -1)
    ok=$(awk -v v="$v" -v r="$ref" -v t="$tol" 'BEGIN{d=v-r; if(d<0)d=-d; print (v != "" && d<=t) ? "PASS" : "FAIL"}')
    n=$((n+1)); [ "$ok" = PASS ] || fail=$((fail+1))
    printf "%-5s %-28s %+.10f  ref %+.10f  diff %+.1e  tol %.0e\n" "$ok" "$name" "$v" "$ref" "$(awk -v v="$v" -v r="$ref" 'BEGIN{print v-r}')" "$tol"
}
democheck demo-one-loop  1  0.5            3e-7  1.0 1e-12 all1
democheck demo-two-loop  2 -0.328478965579 1e-6  1.0 1e-12 all2
echo "$((n-fail))/$n passed  (outputs in $OUT)"
[ $fail -eq 0 ]
