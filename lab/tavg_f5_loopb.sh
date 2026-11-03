#!/bin/bash
# ciclo non sorvegliato della fase 5: una corsa per processo (<= 590 s) finche' tutto e' fatto
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for i in $(seq 1 100); do
  out=$(timeout 590 python -W ignore tavg_f5_drive.py --rev 2>&1); echo "$out" >> out/tavg_f5b.log
  echo "$out" | grep -q TUTTO_FATTO && { echo "== FINITO $(date)" >> out/tavg_f5b.log; exit 0; }
done
