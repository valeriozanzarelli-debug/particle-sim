#!/bin/bash
# colla_drive.sh JOB1 JOB2 ... : run each COLLA job in <= 8-minute resumable chunks until its result json exists.
cd "$(dirname "$0")"
for j in "$@"; do
  n=0
  while [ ! -f "out/colla_$j.json" ] && [ $n -lt 12 ]; do
    python -W ignore colla_run.py "$j" --wall 480 >> out/colla_logs/drive.log 2>&1 || { echo "ERROR $j" >> out/colla_logs/drive.log; break; }
    n=$((n+1))
  done
done
echo "DRIVER DONE $*" >> out/colla_logs/drive.log
