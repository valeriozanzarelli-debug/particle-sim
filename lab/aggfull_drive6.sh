#!/bin/bash
# phase 6 driver: the T-level s states (3-D, for the coarse-level alpha^4 exam), each call <= 590 s, repeated until converged
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for job in eT8_1s eT8_2s eT16_1s; do
  for i in 1 2 3 4; do
    timeout 590 python -u aggfull_run.py $job
    python -c "import json,sys; d=json.load(open('out/aggfull_atomT.json')); sys.exit(0 if d.get('$job',{}).get('converged') else 1)" && break
  done
done
