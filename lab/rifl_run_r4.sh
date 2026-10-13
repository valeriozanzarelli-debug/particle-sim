#!/bin/bash
# corse Q1 con griglia geometrica (pre-registrazione 2) e Q2long (N grandi a theta grande); ogni pezzo < 10 min
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")/.."
for e in 2.5 5 10 20 50 100 200 500 1000 3000 10000; do timeout 590 python3 lab/rifl_r3.py Q1 $e; done
timeout 590 python3 lab/rifl_r3.py Q1sum
for th in 1.5 1.3 1.0 0.6; do for w in fisica intera; do for rep in 1 2; do timeout 590 python3 lab/rifl_r3.py Q2long $th $w; done; done; done
