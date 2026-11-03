#!/bin/bash
# rilancia acqua_opt.py a pezzi da <= 590 s finche' la fase non scrive DONE (cache delle valutazioni = checkpoint)
cd /home/claude/particle-sim/lab; export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1
for f in "$@"; do
  set -- $f
  for i in $(seq 1 40); do
    timeout 590 python acqua_opt.py $(echo $f | tr ':' ' ') && break
  done
done
