#!/bin/bash
# corse E0, E1 (una eps alla volta), E4 dello studio RIFLUSSO; ogni pezzo < 10 min
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")/.."
for e in 2.5 2.0 1.6 1.0 0.8 0.5 0.1; do timeout 590 python3 lab/rifl_r1.py E1 $e; done
timeout 590 python3 lab/rifl_r1.py E1sum
timeout 590 python3 lab/rifl_r1.py E0
timeout 590 python3 lab/rifl_r1.py E4
