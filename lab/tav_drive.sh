#!/bin/bash
# tav_drive.sh <job> ... : rilancia lab/tav_run.py a pezzi da <= 10 minuti finche' ogni job stampa DONE (log in lab/out/tav_<job>.log, non tracciato)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")/.."
for job in "$@"; do
  for i in $(seq 1 60); do
    timeout 590 python3 lab/tav_run.py "$job" 540 >> "lab/out/tav_$job.log" 2>&1
    tail -n 1 "lab/out/tav_$job.log" | grep -q DONE && break
  done
done
