#!/bin/bash
# produzione H2O: pianifica -> 3 lavoratori in parallelo -> ripianifica, finche' la fase e' "fatto"
cd /home/claude/particle-sim/lab; export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
NW=${1:-3}
for it in $(seq 1 30); do
  f=$(python acqua_plan.py next | awk '{print $1}')
  echo "$(date +%T) fase $f" >> out/acqua_prod_log.txt
  [ "$f" = "fatto" ] && break
  rm -f out/acqua_coda_w*.txt
  awk -v n=$NW '{print > ("out/acqua_coda_w" (NR % n) ".txt")}' out/acqua_coda.txt
  for w in out/acqua_coda_w*.txt; do ./acqua_pts.sh $w & done
  wait
done
