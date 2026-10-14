#!/bin/bash
# corse a pezzi: rilancia "python <script> <args>" finche' non esce con codice diverso da 3 (ogni pezzo <= 590 s)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")"
while true; do
  timeout 590 python "$@"
  rc=$?
  if [ $rc -ne 3 ] && [ $rc -ne 124 ]; then echo "DRIVE FINE rc=$rc"; break; fi
done
