#!/bin/bash
# usage: bash zoomknot_drive.sh "CONFIG MODE" ...   (sequential; re-invokes each run until its json exists)
cd "$(dirname "$0")"
mkdir -p out/zoomknot_logs
for job in "$@"; do
  set -- $job
  tag=$1; [ "$2" = "learned" ] && tag=${1}_learned
  for i in 1 2 3 4 5; do
    [ -f out/zoomknot_${tag}.json ] && break
    python -u zoomknot_run.py $job >> out/zoomknot_logs/${tag}.log 2>&1
  done
done
