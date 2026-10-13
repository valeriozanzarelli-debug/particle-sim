#!/bin/bash
# VIDEO IDROGENO: render all frames in resumable pieces, then the mp4s and stills
L=$(cd "$(dirname "$0")" && pwd)
NV=$(python3 -c "import sys; sys.path.insert(0,'$L'); import vl_idro as V; print(V.NV)")
D=${VL_SCR:-/tmp/claude-0/-home-claude-particle-sim/a13a264e-04a0-55c0-8c14-0e293e6cf6fd/scratchpad/vl}/idro/frames
while [ "$(ls $D 2>/dev/null | grep -c '^[0-9]*\.png$')" -lt "$NV" ]; do python3 "$L/vl_idro.py" render ${1:-0} ${2:-$NV} || exit 1; done
