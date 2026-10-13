#!/bin/bash
# VIDEO ATOMO: render frames IV0..IV1 at scale S in resumable pieces (each python call <= ~8 min)
SCR=${VL_SCR:-/tmp/claude-0/-home-claude-particle-sim/a13a264e-04a0-55c0-8c14-0e293e6cf6fd/scratchpad/vl}
a=$1; b=$2; s=${3:-1.0}
w=$(python3 -c "print(int(1920*$s))")
need=$((b-a))
while true; do
  have=$(python3 -c "import os; print(sum(os.path.exists('$SCR/frames_$w/%05d.png'%i) for i in range($a,$b)))")
  [ "$have" -ge "$need" ] && break
  VL_BUDGET=470 python3 "$(dirname "$0")/vl_render.py" batch $a $b --scala $s
  have2=$(python3 -c "import os; print(sum(os.path.exists('$SCR/frames_$w/%05d.png'%i) for i in range($a,$b)))")
  [ "$have2" -eq "$have" ] && sleep 60
done
echo done $a $b
