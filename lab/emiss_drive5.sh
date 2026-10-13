#!/bin/bash
# phase 3: 96^3 levels (pieces), then the 96^3 runs
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for s in "$@"; do
  for i in $(seq 1 20); do
    timeout 590 python emiss_ritz96.py $s 96; rc=$?
    [ $rc -eq 3 ] || break
  done
  echo "=== ritz $s rc=$rc"
done
