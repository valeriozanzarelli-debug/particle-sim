#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for a in 0.5 0.4 0.3 0.25; do
  until timeout 590 python surf_f4.py prep $a; do [ $? -eq 3 ] || { echo FAILED prep $a; break; }; done
  until timeout 590 python surf_f4.py run $a; do [ $? -eq 3 ] || { echo FAILED run $a; break; }; done
done
echo F4DONE
