#!/bin/bash
# VIDEO IDROGENO (fase 2): whole pipeline, resumable - fields, balls, exams, frames, mp4s and stills
L=$(cd "$(dirname "$0")" && pwd)
S=${VL_SCR:-/tmp/claude-0/-home-claude-particle-sim/a13a264e-04a0-55c0-8c14-0e293e6cf6fd/scratchpad/vl}/idro2
until [ -f $S/campi.json ]; do python3 "$L/vl_idro.py" campi || exit 1; done
[ -f $S/fermo.npz ] || python3 "$L/vl_idro.py" palline || exit 1
[ -f $S/esami.done ] || { python3 "$L/vl_idro.py" esami && touch $S/esami.done; } || exit 1
NV=$(python3 -c "import sys; sys.path.insert(0,'$L'); import vl_idro as V; print(V.NV)")
D=$S/frames
while [ "$(ls $D 2>/dev/null | grep -c '^[0-9]*\.png$')" -lt "$NV" ]; do python3 "$L/vl_idro.py" render ${1:-0} ${2:-$NV} || exit 1; done
python3 "$L/vl_idro.py" video
