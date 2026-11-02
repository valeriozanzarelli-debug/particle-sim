#!/bin/bash
# VIDEO ATOMO, tutto con un comando: rifà ogni passo dall'esportazione in /mnt/project-files/ripresa/dati/atomo_accoppiato.
# Se l'esportazione è cambiata (impronta sha256 di manifest + fotogrammi.csv), cancella i prodotti intermedi e riparte.
#   lab/vl_tutto.sh            (VL_ETICHETTA=BOZZA di default; VL_ETICHETTA="" per la versione finale)
set -e
export VL_SIG=${VL_SIG:-0} VL_FLUSSO=${VL_FLUSSO:-corrente}   # emendamento 3 (lab/out/vl_amend3.json)
L=$(cd "$(dirname "$0")" && pwd)
SCR=${VL_SCR:-/tmp/claude-0/-home-claude-particle-sim/a13a264e-04a0-55c0-8c14-0e293e6cf6fd/scratchpad/vl}
mkdir -p "$SCR"
IMP=$(python3 -c "import sys; sys.path.insert(0,'$L'); from vl_comune import impronta; print(impronta()[0])")
if [ "$(cat "$SCR/impronta.txt" 2>/dev/null)" != "$IMP" ]; then
  echo "esportazione nuova ($IMP): riparto da zero"
  rm -rf "$SCR"/el_* "$SCR"/lu_* "$SCR"/*_ckpt.npz "$SCR"/*_part.npz "$SCR"/*_DONE.json "$SCR"/frames_1920 "$SCR"/vl_refs.json
  python3 "$L/vl_campi.py" > /dev/null
  echo "$IMP" > "$SCR/impronta.txt"
fi
"$L/vl_drive.sh" elettrone & "$L/vl_drive.sh" luce & wait
python3 -c "
import numpy as np; d=np.load('$SCR/el_dati.npy',mmap_mode='r')
np.save('$SCR/el_dipolo.npy', np.array([float(np.asarray(d[k,2],np.float64).mean()) for k in range(201)]))"
python3 "$L/vl_esami.py" > /dev/null
[ -f "$SCR/vl_refs.json" ] || python3 "$L/vl_render.py" refs 0 50
for r in "0 374" "374 748" "748 1121"; do "$L/vl_render_drive.sh" $r 1.0 & done; wait
python3 "$L/vl_render.py" ritocca 0 1121      # info panel from the current code (resumable)
"$L/vl_video.sh"
python3 "$L/vl_figura.py"
