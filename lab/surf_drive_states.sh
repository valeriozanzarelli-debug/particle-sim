#!/bin/bash
# runs surf_states.py in <= 10 min pieces until both states are done
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for s in 1s 2p; do
  for k in $(seq 1 12); do
    timeout 590 python surf_states.py $s && break
  done
done
