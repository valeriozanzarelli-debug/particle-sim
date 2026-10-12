#!/bin/sh
cd "$(dirname "$0")/.."
for i in 1 2 3 4 5; do
  timeout 590 python agg_zz_zoom.py >> out/agg_zz_zoom.log 2>&1
  grep -q "hybrid 1s:" out/agg_zz_zoom.log && break
done
echo ZOOMDONE >> out/agg_zz_zoom.log
