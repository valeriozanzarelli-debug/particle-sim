#!/bin/bash
# listening of the 96^3 levels in pieces (exit 3 = run again)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for sg in 1s+:0.89955 2p3/2+3:0.97340; do
  s=${sg%%:*}; export GUESS=${sg##*:}
  for i in 1 2 3 4 5; do
    timeout 590 python emiss_levels.py 96 $s 1000; rc=$?
    [ $rc -eq 3 ] || break
  done
done
