#!/bin/sh
cd "$(dirname "$0")/.."
for ch in A B; do
  for i in 1 2 3 4 5 6 7 8; do
    timeout 590 python agg_zz_3d.py $ch >> out/agg_zz_3d_$ch.log 2>&1
    if grep -q "hierarchy complete" out/agg_zz_3d_$ch.log; then break; fi
  done
done
echo ALLDONE >> out/agg_zz_3d_B.log
