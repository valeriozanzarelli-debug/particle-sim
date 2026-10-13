#!/usr/bin/env bash
# Studio "bagliori": una foto del Laboratorio senza schermo a un tick preciso, poi gli esami I1-I3.
# uso: lab/gl_foto.sh <nome> <tick> -- <argomenti di atomo>     (env ATOMO_BAGLIORE per provare il bagliore)
set -euo pipefail
nome=$1; tick=$2; shift 3
dir=lab/out/gl_frames
mkdir -p "$dir"
timeout 300 xvfb-run -a -s "-screen 0 1280x720x24" target/release/atomo "$@" --foto "$dir/$nome.png" \
    --sequenza 1 --tick "$tick" --frames 999999 2>&1 | grep -E "sequenza|panic" || true
python3 lab/gl_immagini.py "$dir/${nome}_0000.png"
