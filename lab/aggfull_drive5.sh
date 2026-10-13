#!/bin/sh
# phase 5 driver: clean-solver 2p at a0/8 (N = 256) with three kinetics, and their first-order jobs; each call <= 590 s
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for j in e8_2p_c eK8_2p_c eT8_2p_c q8_2p_c qK8_2p_c qT8_2p_c; do
  for t in 1 2 3 4; do
    out=$(timeout 590 python aggfull_run.py $j 2>&1 | grep -v Warn | grep -v "X = np")
    echo "$out" | tail -3
    echo "$out" | grep -q "converged True\|done already\|dV eigenvalues\|needs" && break
  done
done
echo DRIVER5_END
