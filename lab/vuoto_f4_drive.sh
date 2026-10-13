#!/bin/sh
# phase-4 driver: reruns lab/$1 in pieces of <= 10 min until it says ALL DONE (exit 3/124 = rerun)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")"
i=0; while [ $i -lt 60 ]; do i=$((i+1)); timeout 590 python3 -u "$1"; rc=$?; [ $rc -eq 3 ] || [ $rc -eq 124 ] || break; done
echo "driver exit $rc"
