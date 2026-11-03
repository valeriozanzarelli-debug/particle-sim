#!/bin/bash
# esegue i punti elencati in un file (una riga = argomenti di acqua_run.py), pezzi da <= 590 s, riprende dai checkpoint
cd /home/claude/particle-sim/lab; export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while read -r line; do
  [ -z "$line" ] && continue
  for i in $(seq 1 30); do
    timeout 590 python acqua_run.py $line > /dev/null 2>> out/acqua_pts_err.txt; rc=$?
    [ $rc -eq 0 ] && break
    echo "$(date +%T) rc=$rc $line" >> out/acqua_pts_err.txt
  done
done < "$1"
echo "FINE $1" >> out/acqua_pts_err.txt
