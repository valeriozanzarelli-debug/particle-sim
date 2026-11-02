#!/bin/bash
# VIDEO ATOMO driver: resumable pieces (each python call <= ~8 min) until DONE
SCR=${VL_SCR:-/tmp/claude-0/-home-claude-particle-sim/a13a264e-04a0-55c0-8c14-0e293e6cf6fd/scratchpad/vl}
kind=$1
while [ ! -f "$SCR/${kind}_DONE.json" ]; do
  python3 "$(dirname "$0")/vl_palline.py" "$kind" || exit 1
done
