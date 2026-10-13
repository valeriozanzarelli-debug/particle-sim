#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
until timeout 590 python surf_f4.py prep 0.2; do [ $? -eq 3 ] || { echo FAILED prep; break; }; done
until timeout 590 python surf_f4.py run 0.2 600; do [ $? -eq 3 ] || { echo FAILED run; break; }; done
echo F5DONE
