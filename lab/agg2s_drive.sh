#!/bin/bash
# agg2s_drive.sh - runs the agg2s_ production jobs in order, each call <= 590 s, resuming until converged.
# usage: ./agg2s_drive.sh JOB [JOB ...]   (special: "contact" = agg2s_contact.py until all 6 states are converged)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for job in "$@"; do
  for try in 1 2 3 4 5 6 7 8; do
    if [ "$job" = "contact" ]; then
      n=$(python3 -c "import json,os; f='out/agg2s_contact.json'; d=json.load(open(f)) if os.path.exists(f) else {}; print(sum(1 for v in d.values() if v.get('converged')))")
      [ "$n" = "6" ] && break
      timeout 590 python3 agg2s_contact.py 2>&1 | grep -v -i warn
    else
      ok=$(python3 -c "import json,os; f='out/agg2s_atom.json'; d=json.load(open(f)) if os.path.exists(f) else {}; print(d.get('$job',{}).get('converged',False))")
      [ "$ok" = "True" ] && break
      timeout 590 python3 agg2s_run.py "$job" 2>&1 | grep -v -i warn
    fi
  done
done
echo DRIVE DONE "$@"
