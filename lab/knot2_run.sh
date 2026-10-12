#!/bin/bash
# knot2_run.sh SCRIPT ARGS... : rerun a resumable knot2 script (exit 3 = checkpoint) until it finishes.
export KNOT_SCRATCH="${KNOT_SCRATCH:-C:/Users/valer/AppData/Local/Temp/claude/C--Users-valer-Desktop-Programma-di-simulazione-particle-sim/2978b0e8-1bc3-443b-89d8-2e1ba7531c9b/scratchpad}"
mkdir -p lab/out/knot2_logs
s=$1; shift
log="lab/out/knot2_logs/$(basename $s .py)_$(echo "$@" | tr ' ' '_').log"
[ -n "$KNOT_APPEND" ] || : > "$log"
while true; do
  timeout 400 python -u lab/$s "$@" >> "$log" 2>&1; rc=$?
  [ $rc -eq 3 ] || break
done
echo "FINISHED $s $@ rc=$rc"
