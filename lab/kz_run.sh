#!/bin/bash
# K1 production queue (sequential, 2 threads).  log: out/kz_k1.log
cd "$(dirname "$0")"
run() { python kz_k1.py "$@" >> out/kz_k1.log 2>&1; }
for t in 64 128 256 512 1024; do run 48 $t 0 8; done
run 48 2048 0 4; run 48 2048 4 4
run 48 4096 0 2; run 48 4096 2 2; run 48 4096 4 2; run 48 4096 6 2
for t in 64 128 256 512 1024 2048 4096; do run 32 $t 0 8; done
for t in 64 128 256 512 1024 2048; do run 48 $t 0 6 0.5; done
run 48 4096 0 3 0.5; run 48 4096 3 3 0.5
for t in 256 512; do run 64 $t 0 4; done
run 64 1024 0 2; run 64 1024 2 2
for s in 0 1 2 3; do run 64 2048 $s 1; done
for s in 0 1 2 3; do run 64 4096 $s 1; done
echo DONE >> out/kz_k1.log
