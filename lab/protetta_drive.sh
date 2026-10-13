#!/bin/sh
# reruns lab/protetta_vestito.py in pieces of <= 10 min until ALL DONE (exit 3/124 = rerun)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")"
i=0; while [ $i -lt 40 ]; do i=$((i+1)); timeout 590 python3 -u protetta_vestito.py; rc=$?; [ $rc -eq 3 ] || [ $rc -eq 124 ] || break; done
echo "driver exit $rc"
