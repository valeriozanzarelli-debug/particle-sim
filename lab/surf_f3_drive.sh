#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for k in learned bare; do
  until timeout 590 python surf_f3.py prep $k; do [ $? -eq 3 ] || { echo FAILED prep $k; break; }; done
done
for k in learned bare; do
  until timeout 590 python surf_f3.py run $k; do [ $? -eq 3 ] || { echo FAILED run $k; break; }; done
done
echo F3DONE
