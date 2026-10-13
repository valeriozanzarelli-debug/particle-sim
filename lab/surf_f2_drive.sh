#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for n in 48 96; do for M in 1 4 16 25; do
  [ -f out/surf_f2_val_${n}_${M}.json ] || timeout 590 python surf_f2_valentini.py $n $M
done; done
echo VALDONE
