"""emiss_x5d.py - declared extra diagnosis for X5: Ritz from the P4_E1pi final state in [0.86, 0.92] (learned nucleus,
48^3, 2400 ticks, filters every 0.0025): every level with its weight in the state (down to 1e-8) and its dipole to the
2p states.  Run: python lab/emiss_x5d.py -> lab/out/emiss_x5d.json"""
import json, os
import numpy as np
from emiss_core import CK, OUT, norm2
import emiss_ritz as RZ
from emiss_v4_dev import make_atom, S_LEARNED
os.environ["EMISS_CFG"] = "none"
import emiss_configs as EC

d = np.load(os.path.join(CK, "run_P4_E1pi.npz"))
psi = d["psi"].copy()
psi /= np.sqrt(norm2(psi))
D = make_atom(n=48, s=S_LEARNED)
RZ.NT = 2400
out = RZ.ritz(D, psi, 0.86, 0.92, 0.0025)
p12, p32 = EC.level4("2p1/2+"), EC.level4("2p3/2+")
rows = []
for w, res, wt, v in sorted(out, key=lambda r: r[0]):
    if not (0.85 < w < 0.93):
        continue
    v = v.reshape(psi.shape)
    d12 = float(np.sum(np.abs(D.dip_matrix(v, p12))**2))
    d32 = float(np.sum(np.abs(D.dip_matrix(v, p32))**2))
    rows.append(dict(omega=float(w), residual=float(res), weight=float(wt), d2_to_2p12=d12, d2_to_2p32=d32))
    print(f"omega {w:.5f} residual {res:.1e} weight {wt:.2e}  |d|^2 to 2p1/2 {d12:.3f}, to 2p3/2 {d32:.3f}")
json.dump(rows, open(os.path.join(OUT, "emiss_x5d.json"), "w"), indent=1)
