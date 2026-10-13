#!/bin/bash
# amendments A1, A2 (lab/zstar_amend.py) and the film export
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")"
for a in 0.003 0.001 0.0003 0.0001; do
  for v in S0.25:fine R+i:fine; do timeout 590 python3 zstar_run.py $v $a; done
done
for a in 0.001 0.0003 0.0001; do timeout 590 python3 zstar_run.py K+2:fine $a; done
for a in 0.003 0.001 0.0003 0.0001; do timeout 590 python3 zstar_run.py K-2:fine $a; done
timeout 590 python3 zstar_film.py
