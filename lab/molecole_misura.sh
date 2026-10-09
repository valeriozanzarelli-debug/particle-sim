#!/bin/sh
# (dev profile = opt-level 3, no LTO: same speed, much faster to build)
# Runs every measurement of the molecule level (each invocation < 10 min; molecule scans resume from the
# points already stored in assets/molecole/risultati.txt). Log: lab/out/molecole_misura.log
cd "$(dirname "$0")/.." || exit 1
cargo build --bin molecole || exit 1
for t in idrogeno eccitati nodo2s nodo3p campo riga h2piu h2piu h2 h2 h2 h2 elio elio1p elio1p eliocampo eliocampo; do
  timeout 590 ./target/debug/molecole --misura "$t" 2>&1 | tee -a lab/out/molecole_misura.log
done
