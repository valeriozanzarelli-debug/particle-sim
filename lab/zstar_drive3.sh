#!/bin/bash
# amendment A3: Z_cr by the sign method for every variant (validation + the narrow-window ones), a down to 1e-5
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$(dirname "$0")"
A="3e-3 1e-3 3e-4 1e-4 3e-5 1e-5"
timeout 590 python3 zstar_run_sign.py base 131.0 134.0 0.01 3e-5 1e-5
timeout 590 python3 zstar_run_sign.py R+1 120 130 0.01 $A
timeout 590 python3 zstar_run_sign.py K+1 120 130 0.01 $A
timeout 590 python3 zstar_run_sign.py K+1R+1 128 136 0.01 $A
timeout 590 python3 zstar_run_sign.py R-i 125 133 0.01 $A
timeout 590 python3 zstar_run_sign.py S0.25 102.7 106 0.002 $A
timeout 590 python3 zstar_run_sign.py R+i 94.1 101 0.002 $A
timeout 590 python3 zstar_run_sign.py S1 137 160 0.01 $A
timeout 590 python3 zstar_run_sign.py S0.866 137 160 0.01 $A
timeout 590 python3 zstar_run_sign.py S0.75 137 160 0.01 $A
timeout 590 python3 zstar_run_sign.py split 131 134 0.01 $A
timeout 590 python3 zstar_run_sign.py K-2 195 290 0.01 $A
timeout 590 python3 zstar_run_sign.py K+2 41 44 0.0005 1e-3 3e-4 1e-4 3e-5 1e-5
ZSTAR_OMEGA=0 timeout 590 python3 zstar_run_sign.py K+2 41 44 0.0005 1e-3 3e-4 1e-4 3e-5 1e-5
