#!/bin/bash
# phase 7 driver: glue and chain runs, each call <= 590 s, repeated until 'all done'
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for s in aggfull_p7_glue.py aggfull_p7_chain.py; do
  for i in 1 2 3 4 5 6; do
    timeout 590 python -u $s run | tee /dev/stderr | grep -q "all done" && break
  done
done
