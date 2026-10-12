#!/bin/bash
# K2 and K3 queue, starts after the K1 queue has written DONE.  logs: out/kz_k2.log, out/kz_k3.log
cd "$(dirname "$0")"
until grep -q DONE out/kz_k1.log 2>/dev/null; do sleep 30; done
for a in "8 20000" "12 20000" "16 20000" "24 12000" "32 8000"; do python kz_k2.py $a >> out/kz_k2.log 2>&1; done
python kz_k2.py fit >> out/kz_k2.log 2>&1
python kz_k3.py > out/kz_k3.log 2>&1
echo DONE2 >> out/kz_k3.log
