#!/bin/bash
# run jobs in order, repeating each while it checkpoints (exit 3)
cd "$(dirname "$0")"
for j in "$@"; do
  while true; do python3 emiss6_run.py "$j" 540; c=$?; [ $c -eq 3 ] || break; done
  echo "== $j exit $c"
done
