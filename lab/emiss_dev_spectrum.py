"""DEVELOPMENT (declared, not a test): spectrum of <psi0|psi(t)> for start states built from the radial rule,
in 48^3 and 96^3 boxes (a = 0.5), to see which levels the 3-D rule has near 1s and 2p."""
import sys, time, json
import numpy as np
from emiss_core import *

n = int(sys.argv[1]); which = sys.argv[2]; ticks = int(sys.argv[3]); ar = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
a = 0.5
Za = 60 * ALPHA
D = Dirac3D(n, a)
kap, nr = dict(s=(-1, 1), p12=(1, 2), p32=(-2, 2), s2=(-1, 2))[which]
om_r, rg, g, f = radial_state(kap, nr, ar, Za)
psi = dirac_3d_from_radial(rg, g, f, kap, 0.5, D.x)
p0 = psi.copy()
C = []
t0 = time.time()
for t in range(ticks):
    if t % 2 == 0:
        C.append(ovl(p0, psi))
    D.tick(psi, t)
C = np.array(C)
T = len(C)
tt = np.arange(T) * 2 * a
h = np.sin(np.pi * (np.arange(T) + 0.5) / T)**2
om = np.linspace(0.80, 1.0, 40001)
S = np.abs(np.exp(1j * np.outer(om, tt)) @ (h * C))**2
S /= S.max()
pk = [(om[i], S[i]) for i in range(1, len(om) - 1) if S[i] > S[i - 1] and S[i] > S[i + 1] and S[i] > 1e-3]
pk.sort(key=lambda p: -p[1])
print(f"n={n} {which} (radial a={ar}: omega {om_r:.5f}, Dirac {dirac_formula(nr, kap, Za):.5f}) {ticks} ticks {time.time()-t0:.0f}s")
print("peaks:", ", ".join(f"{p[0]:.5f} ({p[1]:.3f})" for p in pk[:8]))
