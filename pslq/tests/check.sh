#!/bin/sh
# Quick check of the whole chain: run_pslq.sh at coarse steps, values against the closed forms
# (to ~1e-9, the h error at h = 0.7).  No Python.
DIR=$(cd "$(dirname "$0")/.." && pwd)
OUT=${TMPDIR:-/tmp}/sincqed_pslq_test.$$
QUICK=1 OUT="$OUT" sh "$DIR/run_pslq.sh" > /dev/null || { echo "FAIL: run_pslq.sh"; exit 1; }
# closed forms (pi, log 2, zeta 3 to 20 digits)
awk -F'\t' 'BEGIN {
    pi2 = 9.8696044010893586188; l2 = 0.69314718055994530942; z3 = 1.2020569031595942854
    ref["ladder"] = 11/48 + pi2/18; ref["self-energy"] = -1/24 - pi2/18; ref["self-energy_demo"] = ref["self-energy"]
    ref["corner"] = -55/24 + pi2/18 + pi2*l2/3 - z3/2; ref["crossed"] = 1/6 + 13*pi2/36 - 5*pi2*l2/6 + 5*z3/4
    ref["Pi2"] = 119/36 - pi2/3; ref["sum"] = 197/144 + pi2/12 - pi2*l2/2 + 3*z3/4; tol = 5e-9; fail = 0 }
    NR > 1 { split($2, p, ","); v = p[1] + p[2]; d = v - ref[$1]; if (d < 0) d = -d
             ok = d <= tol ? "PASS" : "FAIL"; if (ok == "FAIL") fail++
             printf "%-5s %-17s %+.12f  closed form %+.12f  diff %.1e\n", ok, $1, v, ref[$1], d }
    END { print fail ? fail " check(s) failed" : "all checks passed"; exit fail > 0 }' "$OUT/values.tsv"
st=$?; rm -rf "$OUT"; exit $st
