#!/bin/sh
cd "$(dirname "$0")/.."
for i in $(seq 1 12); do
  timeout 590 python agg_zz_atom.py >> out/agg_zz_atom.log 2>&1
  python -c "
import json,sys
d=json.load(open('out/agg_zz_atom.json'))
sys.exit(0 if len([k for k,v in d.items() if v.get('converged')])>=4 else 1)" && break
done
echo ATOMDONE >> out/agg_zz_atom.log
