#!/bin/bash
# secondo tentativo: 3 lavoratori, poi giudizio
cd /home/claude/particle-sim/lab
printf "calc AB 3 2.88 4.8 -55.3 1e-06\ncalc A 3 2.88 4.8 -55.3 1e-06\ncalc B 3 2.88 4.8 -55.3 1e-06\n" > out/acqua_d2coda_w0.txt
printf "calc AB 3 2.92 4.8 -55.3 1e-06\ncalc A 3 2.92 4.8 -55.3 1e-06\ncalc B 3 2.92 4.8 -55.3 1e-06\n" > out/acqua_d2coda_w1.txt
: > out/acqua_d2coda_w2.txt
for R in 2.88 2.92; do for k in AB A B; do echo "calc $k 2 $R 4.8 -55.3 1e-06" >> out/acqua_d2coda_w2.txt; done; done
for w in out/acqua_d2coda_w*.txt; do ./acqua_dpts.sh $w & done; wait
python acqua_dimero2.py > out/acqua_dimero2_judge.out 2>&1
