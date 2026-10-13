#!/usr/bin/env bash
# Studio "bagliori": una clip al rallentatore dal Laboratorio senza schermo (xvfb + llvmpipe), un fotogramma per
# tick della regola (--sequenza), poi ffmpeg. Uso:
#   lab/gl_clip.sh <nome> <fotogrammi> <tick_iniziale> <passo> <fps> -- <argomenti di atomo>
# es. lab/gl_clip.sh fotone_davanti 48 6 1 12 -- --scena fotone --yaw 1.5708 --pitch 0
# Scrive lab/out/gl_frames/<nome>_NNNN.png e lab/out/gl_<nome>.mp4 / .gif (i PNG dei fotogrammi non si committano).
set -euo pipefail
nome=$1; k=$2; t0=$3; passo=$4; fps=$5; shift 6
dir=lab/out/gl_frames
mkdir -p "$dir"
timeout 590 xvfb-run -a -s "-screen 0 1280x720x24" target/release/atomo "$@" --foto "$dir/$nome.png" \
    --sequenza "$k" --tick "$t0" --passo "$passo" --frames 999999 2>&1 | grep -E "sequenza|error|panic" || true
n=$(ls "$dir"/"$nome"_*.png 2>/dev/null | wc -l)
echo "fotogrammi: $n"
ffmpeg -loglevel error -y -framerate "$fps" -i "$dir/${nome}_%04d.png" -vf "crop=1600:760:0:70,format=yuv420p" -c:v libx264 -crf 20 "lab/out/gl_$nome.mp4"
ffmpeg -loglevel error -y -framerate "$fps" -i "$dir/${nome}_%04d.png" -vf "crop=1600:760:0:70,scale=640:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128[p];[b][p]paletteuse" "lab/out/gl_$nome.gif"
ls -la "lab/out/gl_$nome.mp4" "lab/out/gl_$nome.gif"
