#!/bin/sh
# Generate all four-loop graphs into OUTDIR.  Usage: sh gen4all.sh OUTDIR [JOBS] [BATCH]
#  1. all Dirac tables (819) in one parallel pool, longest first, one file each in OUTDIR/.tables
#     (so no table is computed twice and cores stay busy until the last few tables);
#  2. C emission in batches of BATCH graphs (default JOBS), reading the tables from OUTDIR/.tables.
# Both steps are resumable: finished tables and non-empty g4_*.c files are skipped.
# OUTDIR/.tables can be deleted once table_4.h has been written.
OUT=${1:?outdir}; JOBS=${2:-4}; B=${3:-$JOBS}
N=$(python3 -c "import graphs; print(len(graphs.enumerate_graphs(4)))") || exit 1
mkdir -p "$OUT"
echo "== Dirac tables  $(date)"
python3 gen.py 4 "$OUT" --jobs $JOBS --cache-dir "$OUT/.tables" --dirac-only --missing || exit 1
i=0
while [ $i -lt $N ]; do
  j=$((i + B)); [ $j -gt $N ] && j=$N
  miss=0; k=$i
  while [ $k -lt $j ]; do f=$(printf "%s/g4_%02d.c" "$OUT" $k); [ -s "$f" ] || miss=1; k=$((k + 1)); done
  if [ $miss = 1 ]; then
    echo "== emit $i:$j  $(date)"
    python3 gen.py 4 "$OUT" --range $i:$j --jobs $JOBS --cache-dir "$OUT/.tables" --missing || exit 1
  fi
  i=$j
done
python3 gen.py 4 "$OUT" --table-only
echo "== done $(date)"
