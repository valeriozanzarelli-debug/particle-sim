#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while pgrep -f surf_f2_drive2.sh > /dev/null; do sleep 20; done
for i in 0 1 2 3 4 5; do timeout 590 python surf_f2_posthoc.py $i; done
echo POSTDONE
