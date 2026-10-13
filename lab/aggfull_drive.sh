#!/bin/sh
# driver of the production jobs (each call <= 590 s; repeats until the job reports converged/done)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for j in e8_1s e8_2s e8_2p p8_1s p8_2s p8_2p e16_1s p16_1s q8_1s q8_2s q8_2p q16_1s m8_1s e8_1s_N128 e8_1s_N192 e8_2s_N192 ca8_1s_N192 ca8_2s_N192 ca8_1s ca8_2s; do
  for t in 1 2 3 4 5; do
    out=$(timeout 590 python aggfull_run.py $j 2>&1 | grep -v Warn | grep -v "X = np")
    echo "$out" | tail -6
    echo "$out" | grep -q "converged True\|done already\|M_atom\|dV eigenvalues" && break
  done
done
echo DRIVER_END
