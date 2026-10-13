#!/bin/bash
# emendamento A9: rifare E2, E3 e Q2 (N <= 400) con la griglia corretta, poi r5; ogni pezzo < 10 min
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")/.."
python3 - <<'PY'
import json; f='lab/out/rifl_r2.json'; d=json.load(open(f)); d['prima_di_A9']={k: d.pop(k) for k in list(d) if k!='tempi'}; json.dump(d,open(f,'w'),indent=1)
PY
for U in 10 14.142 20 28.284 40; do timeout 590 python3 lab/rifl_r2.py E2 $U; done
timeout 590 python3 lab/rifl_r2.py E2sum
for e in 10 5 2.5 1 0.5; do timeout 590 python3 lab/rifl_r2.py E3 $e; done
timeout 590 python3 lab/rifl_r2.py E3sum
for th in 0.01 0.03 0.1 0.3 0.6 1.0 1.3 1.5; do timeout 590 python3 lab/rifl_r3.py Q2redo $th; done
timeout 590 python3 lab/rifl_r3.py Q2sum
timeout 590 python3 lab/rifl_r5.py
