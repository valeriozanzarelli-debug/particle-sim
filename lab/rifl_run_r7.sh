#!/bin/bash
# Q2osc: scansione fitta in N della banda intera a theta = 1.5 (e 1.3); ogni pezzo < 10 min
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")/.."
for N0 in 1000 2000 5000 10000; do timeout 590 python3 lab/rifl_r3.py Q2osc 1.5 $N0 2; done
for rep in 1 2 3 4 5; do timeout 590 python3 lab/rifl_r3.py Q2osc 1.5 20000 4; done
for N0 in 1000 2000 4000; do timeout 590 python3 lab/rifl_r3.py Q2osc 1.3 $N0 1; done
