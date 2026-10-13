#!/bin/sh
# phase 4 driver: 2p at a0/8 with the float64 field (N = 192) and three kinetics; each call <= 590 s
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for j in e8_2p_N192 q8_2p_N192 eK8_2p_N192 qK8_2p_N192 eT8_2p_N192 qT8_2p_N192; do
  for t in 1 2 3 4 5; do
    out=$(timeout 590 python aggfull_run.py $j 2>&1 | grep -v Warn | grep -v "X = np")
    echo "$out" | tail -3
    echo "$out" | grep -q "converged True\|done already\|dV eigenvalues" && break
  done
done
echo DRIVER4_END
