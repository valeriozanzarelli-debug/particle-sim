#!/bin/bash
# driver for the ZSTAR production runs (each python call < 10 min; checkpoint in lab/out/zstar_runs.json)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")"
timeout 590 python3 zstar_run.py order
for a in 0.003 0.001 0.0003 0.0001; do
  for v in R+1 S0.25 K+1 K+1R+1 R-i R+i S1 S0.866 S0.75 split; do
    timeout 590 python3 zstar_run.py $v $a
  done
done
timeout 590 python3 zstar_run.py refl1 0.001
for a in 0.003 0.001 0.0003 0.0001; do timeout 590 python3 zstar_run.py K-2 $a; done
for a in 0.001 0.0003 0.0001; do timeout 590 python3 zstar_run.py K+2 $a; done
