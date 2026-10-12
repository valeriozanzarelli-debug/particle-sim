#!/bin/sh
# Laboratorio Atomo: the long validation run (several minutes on 4 cores).
# Numbers -> assets/atomo/risultati.txt (read by the viewer and compiled into it), log -> lab/out/atomo_misura.log
cd "$(dirname "$0")/.." || exit 1
mkdir -p lab/out
cargo build --bin atomo -j 2 || exit 1
RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-3} ./target/debug/atomo --misura tutto 2>&1 | tee lab/out/atomo_misura.log
# the long test too
cargo test --bin atomo -j 2 -- --ignored --nocapture 2>&1 | tail -20 | tee -a lab/out/atomo_misura.log
