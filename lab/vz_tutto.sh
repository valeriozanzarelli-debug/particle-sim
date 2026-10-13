#!/bin/bash
# VIDEO SUPERCRITICO, tutto con un comando (dati: lab/out/zstar_film.json al commit fc8a08c, letto con git show)
#   VZ_ETICHETTA="" lab/vz_tutto.sh      (versione finale, senza "BOZZA")
set -e
L=$(cd "$(dirname "$0")" && pwd)
SCR=${VZ_SCR:-/tmp/claude-0/-home-claude-particle-sim/a13a264e-04a0-55c0-8c14-0e293e6cf6fd/scratchpad/vz}
mkdir -p "$SCR"
python3 "$L/vz_esami.py"
[ -f "$SCR/vz_refs.json" ] || python3 "$L/vz_render.py" refs 60
N=$(python3 -c "import sys; sys.path.insert(0,'$L'); from vz_comune import NV; print(NV)")
H=$((N/2))
"$L/vz_render_drive.sh" 0 $H & "$L/vz_render_drive.sh" $H $N & wait
"$L/vz_video.sh"
python3 "$L/vz_figura.py"
