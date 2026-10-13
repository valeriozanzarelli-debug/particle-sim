#!/bin/bash
# runs the given configurations one after the other, each in pieces of <= 590 s (exit 3 = run again)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for c in "$@"; do
  for i in $(seq 1 30); do
    timeout 590 python emiss_run.py $c 540; rc=$?
    [ $rc -eq 3 ] || break
  done
  echo "=== $c finished rc=$rc"
done
