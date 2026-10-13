"""aggfull_dev_pot.py - DEVELOPMENT run (not production, declared): band-limited restriction of -alpha/r vs the old
cell average, with the OLD learned free generator of chain A and gamma0 = 1. Usage: python lab/aggfull_dev_pot.py LEVEL_B N STATE POT"""
import sys, time, json, os
import numpy as np
import aggfull_core as ac
import aggfull_atom as aa
B, N, kind, pot = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
lv = [l for l in ac.load_chain_A() if l["B"] == B][0]
h = B * ac.TH_A
Gk, spin = aa.kinetic_grid(lv, N)
if pot == "bl":
    V, x = ac.coulomb_bl(N, h)
else:
    V, x = ac.coulomb_cellavg(N, h)
at = aa.Atom(Gk, h, V, x)
t0 = time.time()
Y = None
if kind == "2s":
    w, v, r, it = aa.solve(at, at.guess("1s"), maxiter=40)
    Y = v
w, v, r, it = aa.solve(at, at.guess(kind), Y=Y, maxiter=80)
E, res = at.energy(v[:, 0])
n = 1 if kind == "1s" else 2
print(f"B={B} N={N} {kind} {pot}: E = {E*ac.MC2:.6f} eV, rel {E*ac.MC2/(-ac.RY/n**2)-1:+.5e}, resid {res:.1e}, lobpcg {r:.1e} it {it}, {time.time()-t0:.0f}s", flush=True)
