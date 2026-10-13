#!/bin/bash
# VIDEO SUPERCRITICO: draw frames IV0..IV1 in calls of <= 8 min each (resumable: existing frames are skipped)
L=$(cd "$(dirname "$0")" && pwd)
SCR=${VZ_SCR:-/tmp/claude-0/-home-claude-particle-sim/a13a264e-04a0-55c0-8c14-0e293e6cf6fd/scratchpad/vz}
A=$1; B=$2; S=${3:-1.0}
W=$(python3 -c "print(int(1920*$S))")
while true; do
  n=0; for i in $(seq $A $((B-1))); do [ -f "$SCR/frames_$W/$(printf %05d $i).png" ] || n=$((n+1)); done
  [ $n -eq 0 ] && break
  VZ_BUDGET=480 timeout 600 nice -n 10 python3 "$L/vz_render.py" batch $A $B --scala $S || exit 1
done
