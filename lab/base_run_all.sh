#!/bin/bash
# runs every experiment of the base rule, one after the other (each < 5 min on 4 cores); log in out/base_run_all_<ver>.log
cd "$(dirname "$0")"
run() { echo "### $*"; t0=$SECONDS; timeout 570 python "$@" 2>&1 | grep -v "^    \.\.\." ; echo "wall $((SECONDS-t0)) s"; }
{
run base_test.py
run base_conv.py wave
run base_conv.py photons 0.5
run base_conv.py photons 1 2 4
run base_conv.py eps
for e in free double which single grating talbot edge; do run base_exp2d.py $e; done
for e in free gouy slit3d airy arago malus three; do run base_exp3d.py $e; done
run base_mz.py phase 0 90
run base_mz.py phase 180 270
run base_mz.py bomb
run base_mz.py wheeler
} 2>&1
