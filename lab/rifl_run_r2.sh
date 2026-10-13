#!/bin/bash
# corse E2 ed E3 dello studio RIFLUSSO; ogni pezzo < 10 min
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")/.."
for U in 10 14.142 20 28.284 40; do timeout 590 python3 lab/rifl_r2.py E2 $U; done
timeout 590 python3 lab/rifl_r2.py E2sum
for e in 10 5 2.5 1 0.5; do timeout 590 python3 lab/rifl_r2.py E3 $e; done
