#!/bin/bash
# driver A (second version): AN relaxations, then normal modes (C2b), then sextic scan.
cd "$(dirname "$0")"
while tasklist 2>/dev/null | grep -q " 4036 "; do sleep 10; done
COLLA_THREADS=1 ./colla_drive.sh an_B2 an_B4 an_B3
for j in an_B1:12 an_B2:20; do
  job=${j%%:*}; k=${j##*:}; n=0
  while [ ! -f "out/colla_modes_$job.json" ] && [ $n -lt 10 ]; do
    COLLA_THREADS=1 python -W ignore colla_modes.py $job $k --wall 540 >> out/colla_logs/drive.log 2>&1 || { echo "ERROR modes $job" >> out/colla_logs/drive.log; break; }
    n=$((n+1))
  done
done
COLLA_THREADS=1 ./colla_drive.sh sx50_B1 sx50_B2 sx50_B4 sx200_B1 sx200_B2 sx200_B4
