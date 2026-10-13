#!/bin/sh
# second driver: after aggfull_drive.sh, the 2p jobs left and the amendment (K) jobs; each call <= 590 s
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while ! grep -q DRIVER_END "$1" 2>/dev/null; do sleep 30; done
for j in e8_2p q8_2p eK16_1s eK8_1s eK8_2s eK8_2p qK16_1s qK8_1s qK8_2s qK8_2p pK16_1s pK8_1s pK8_2s pK8_2p; do
  for t in 1 2 3 4 5; do
    out=$(timeout 590 python aggfull_run.py $j 2>&1 | grep -v Warn | grep -v "X = np")
    echo "$out" | tail -6
    echo "$out" | grep -q "converged True\|done already\|dV eigenvalues" && break
  done
done
echo DRIVER2_END
