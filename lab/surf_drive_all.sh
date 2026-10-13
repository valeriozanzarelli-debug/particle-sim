#!/bin/bash
# production driver: 1-D/radial runs, then (after the states) the 3-D runs, each in <= 10 min pieces
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while pgrep -f surf_drive_states.sh > /dev/null; do sleep 20; done
for w in f6 cont h1s beat; do
  [ -f out/surf_$w.json ] || timeout 590 python surf_run.py $w
done
for k in $(seq 1 30); do
  [ -f out/surf_orb.json ] && break
  timeout 590 python surf_run.py orb; rc=$?
  [ $rc -ne 3 ] && [ $rc -ne 124 ] && break
done
for s in free 1s 1s2p; do
  for k in $(seq 1 30); do
    [ -f out/surf_a3d_$s.json ] && break
    timeout 590 python surf_run.py a3d $s; rc=$?
    [ $rc -ne 3 ] && [ $rc -ne 124 ] && break
  done
done
echo ALLDONE
