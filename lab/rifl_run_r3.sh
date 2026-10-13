#!/bin/bash
# corse Q1 (Dirac ultrarelativistico) e Q2 (regola a theta finito) dello studio RIFLUSSO; ogni pezzo < 10 min
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")/.."
for e in 2.5 3.5 5 7 10 14 20 28 35 50 70 100; do timeout 590 python3 lab/rifl_r3.py Q1 $e; done
timeout 590 python3 lab/rifl_r3.py Q1sum
for th in 0.01 0.03 0.1 0.3 0.6 1.0 1.3 1.5; do
  for rep in 1 2 3; do timeout 590 python3 lab/rifl_r3.py Q2 $th; done
done
timeout 590 python3 lab/rifl_r3.py Q2sum
