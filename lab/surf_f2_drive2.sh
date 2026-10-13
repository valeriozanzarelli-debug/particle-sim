#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while pgrep -f surf_f2_drive.sh > /dev/null; do sleep 20; done
for w in ritz smooth; do
  for k in $(seq 1 10); do
    timeout 590 python surf_f2_valli.py $w; rc=$?
    [ $rc -ne 3 ] && [ $rc -ne 124 ] && break
  done
done
timeout 590 python surf_f2_valli.py analyze
echo VALLIDONE
