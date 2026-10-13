#!/bin/sh
# third driver: after aggfull_drive2.sh, the old-restriction diagnostics (cell average) left
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while ! grep -q DRIVER2_END "$1" 2>/dev/null; do sleep 30; done
for j in ca8_2s_N192 ca8_1s ca8_2s; do
  for t in 1 2 3 4 5; do
    out=$(timeout 590 python aggfull_run.py $j 2>&1 | grep -v Warn | grep -v "X = np")
    echo "$out" | tail -6
    echo "$out" | grep -q "converged True\|done already\|dV eigenvalues" && break
  done
done
echo DRIVER3_END
