#!/bin/bash
# produzione trimero: pianifica -> lavoratori in parallelo -> ripianifica, finche' "fatto"; poi giudizio
cd /home/claude/particle-sim/lab; export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PROG=acqua_trimero.py TOUT=900
for it in $(seq 1 40); do
  f=$(python acqua_trimero_plan.py next | awk '{print $1}')
  echo "$(date +%T) fase $f" >> out/acqua_tprod_log.txt
  [ "$f" = "fatto" ] && break
  NW=2; [ "$f" = "griglia2" ] && NW=4
  rm -f out/acqua_tcoda_w*.txt
  awk -v n=$NW '{print > ("out/acqua_tcoda_w" (NR % n) ".txt")}' out/acqua_tcoda.txt
  for w in out/acqua_tcoda_w*.txt; do ./acqua_dpts.sh $w & done
  wait
done
python acqua_trimero_plan.py judge > out/acqua_trimero_judge.out 2>&1
