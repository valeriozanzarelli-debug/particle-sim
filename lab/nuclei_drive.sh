#!/bin/bash
# runs the nuclei jobs in chunks of <= 480 s until each one writes its final json
cd "$(dirname "$0")"
export OMP_NUM_THREADS=2
for job in "$@"; do
  n=0
  while [ ! -f "out/nuclei_$job.json" ] && [ $n -lt 12 ]; do
    timeout 590 python nuclei_run.py "$job" --wall 480 || echo "chunk exit $? for $job"
    n=$((n+1))
  done
  echo "JOBDONE $job"
done
