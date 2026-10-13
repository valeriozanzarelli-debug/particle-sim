#!/bin/sh
# rerun lab/vuoto_run3d.py in <= 10 min pieces until it finishes (exit 3 = checkpoint)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")"
while true; do
  timeout 590 python3 vuoto_run3d.py; rc=$?
  [ $rc -eq 3 ] || [ $rc -eq 124 ] || break
done
echo "driver exit $rc"
