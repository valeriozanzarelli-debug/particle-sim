#!/bin/sh
# runs lab/fw1vuoto_run.py with two single-thread workers in <= 10 min pieces until both finish
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")"
run() {
  while true; do
    FW1V_WORKER=$1 timeout 590 python3 fw1vuoto_run.py > /dev/null 2>&1; rc=$?
    [ $rc -eq 3 ] || [ $rc -eq 124 ] || break
  done
  echo "worker $1 exit $rc"
}
run 0 & run 1 & wait
