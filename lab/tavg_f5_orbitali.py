"""tavg_f5_orbitali.py - the computed orbitals for the periodic table (src/atomo/orbatomo.rs, lab/out/orbatomo_prereg.json).

    python lab/tavg_f5_orbitali.py

assets/esportazioni/<O|He>/orbitali.json: {r, orbitali: [{nome, l, P}], spin_orbitali: [{orb, m, s}]}, P(r) = r R(r), bohr.
O: fine-level ROHF, the Hund determinant the energies use. He: P = r sqrt(4 pi rho1) from the ECG's exact density, one per spin.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("OMP_NUM_THREADS", "1")
import tavg_blocks as tb  # noqa: E402
import tavg_radial as tr  # noqa: E402

ASSETS = os.path.join(HERE, "..", "assets")
L = "spdf"


def norma(r, P):
    return float(np.trapezoid(P * P, r))


a = tb.get_ae(8, "1s2 2s2 2p4")
r0 = a.g.r
print("determinante", tr.hund_determinant(a.config))
r = np.geomspace(max(r0[0], 1e-4), 14.0, 400)
orb = []
for i, c in enumerate(a.config):
    P = np.interp(r, r0, a.P[i])
    if P[np.argmax(np.abs(P))] < 0:
        P = -P
    orb.append(dict(nome=f"{c[0]}{L[c[1]]}", l=int(c[1]), P=P.tolist()))
    print(orb[-1]["nome"], "norma", norma(r, P))
so = []
for x in tr.hund_determinant(a.config):
    so.append(x)
print(so[:3])
json.dump(dict(r=r.tolist(), orbitali=orb, spin_orbitali=[dict(orb=int(x[0]), m=int(x[1]), s=int(np.sign(x[2]))) for x in so]),
          open(os.path.join(ASSETS, "esportazioni", "O", "orbitali.json"), "w"))

m = json.load(open(os.path.join(ASSETS, "scene_calcolate", "elio", "manifest.json")))
nr, nz, rmax = m["rho"]["nr"], m["rho"]["nz"], m["rho"]["rmax"]
T = np.fromfile(os.path.join(ASSETS, "scene_calcolate", "elio", "rho1.f32"), "<f4").reshape(nr, nz)
z = np.linspace(-rmax, rmax, nz)
keep = z >= 0
rh, dens = z[keep], T[0, keep].astype(float)
P = rh * np.sqrt(4 * np.pi * dens)
print("He norma", norma(rh, P))
json.dump(dict(r=rh.tolist(), orbitali=[dict(nome="1s", l=0, P=P.tolist())], spin_orbitali=[dict(orb=0, m=0, s=1), dict(orb=0, m=0, s=-1)]),
          open(os.path.join(ASSETS, "esportazioni", "He", "orbitali.json"), "w"))
print("fatto")
