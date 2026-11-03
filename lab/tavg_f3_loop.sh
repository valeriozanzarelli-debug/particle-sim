#!/bin/bash
# ciclo non sorvegliato di un gruppo della fase 3: spezzoni da <= 590 s finche' il gruppo e' finito
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
G=$1; LOG=out/tavg_f3_$G.log
for i in $(seq 1 400); do
  out=$(timeout 590 python tavg_f3_drive.py $G 2>&1); echo "$out" >> $LOG
  echo "$out" | grep -q TUTTO_FATTO && { echo "== $G FINITO $(date)" >> $LOG; exit 0; }
done
