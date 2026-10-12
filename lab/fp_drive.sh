#!/bin/bash
# runs the three polarizations in chunks (each call <= 10 min)
cd "$(dirname "$0")"
for pol in circ+ circ- lin_x; do
  while true; do
    timeout 590 python3 fp_run.py $pol 520; rc=$?
    [ $rc -eq 3 ] && continue
    [ $rc -ne 0 ] && echo "ERRORE $pol rc=$rc" && exit $rc
    break
  done
done
echo FINITO
