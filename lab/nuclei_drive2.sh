#!/bin/bash
# like nuclei_drive.sh but with a chunk cap per job (budget, declared): after NMAX chunks the job is finalized unconverged
cd "$(dirname "$0")"
export OMP_NUM_THREADS=2
NMAX=$1; shift
for job in "$@"; do
  n=$(grep -c "chunk end" "out/nuclei_logs/$job.log" 2>/dev/null || echo 0)
  while [ ! -f "out/nuclei_$job.json" ] && [ $n -lt $NMAX ]; do
    timeout 590 python nuclei_run.py "$job" --wall 480 || echo "chunk exit $? for $job"
    n=$((n+1))
  done
  [ -f "out/nuclei_$job.json" ] || python nuclei_run.py "$job" --finalize
  echo "JOBDONE $job"
done
