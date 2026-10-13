#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while pgrep -f surf_drive_all.sh > /dev/null; do sleep 20; done
for k in $(seq 1 30); do
  [ -f out/surf_beat.json ] && break
  timeout 590 python surf_run.py beat; rc=$?
  [ $rc -ne 3 ] && [ $rc -ne 124 ] && break
done
echo BEATDONE
