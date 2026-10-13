"""aggfull_fineref.py - the FINE LEVEL reference for the exams: the zig-zag rule itself (radial form, lab/zz_core.py
RadialWalk, frozen) at the SAME fine theta as chain A (theta = a0/4096 = 0.033456), infinitely heavy nucleus.
Levels 1s1/2, 2s1/2, 2p1/2, 2p3/2 by inverse iteration on the walk (= listening for an infinite time), as zz_hydrogen.py.
Also the walk at theta_mu (the walk run with the reduced mass, for the recoil exam: Dirac-Coulomb levels scale with the
mass exactly). Writes lab/out/aggfull_fineref.json.  Run: timeout 590 python lab/aggfull_fineref.py"""
import json, os, time
import numpy as np
import scipy.sparse.linalg as spl
from zz_core import RadialWalk, dirac_formula, ALPHA, MC2_EV
import aggfull_core as ac

OUT = ac.OUT
TH = ac.TH_A
A0 = 1 / ALPHA
NJ = {"1s1/2": (1, -1), "2s1/2": (2, -1), "2p1/2": (2, 1), "2p3/2": (2, -2)}
res = dict(theta=TH, levels={})
t0 = time.time()
for nm, (n, kappa) in NJ.items():
    rw = RadialWalk(kappa, TH, (40 if n == 1 else 70) * A0)
    Wm = rw.matrix()
    guess = 1 - 0.5 * 1.002 * ALPHA ** 2 / n ** 2
    w, v = spl.eigs(Wm, k=3, sigma=np.exp(-1j * guess * TH), which="LM", tol=0)
    om = np.sort(-np.angle(w) / TH)
    E = om[np.argmin(np.abs(om - guess))]
    D = dirac_formula(n, kappa)
    res["levels"][nm] = dict(E_bind_eV=(1 - E) * MC2_EV, dirac_bind_eV=(1 - D) * MC2_EV, rel_to_dirac=(1 - E) / (1 - D) - 1,
                             N=rw.N, candidates=((1 - om) * MC2_EV).tolist())
    print(f"{nm}: walk theta={TH:.5f}: binding {(1 - E) * MC2_EV:.7f} eV, Dirac {(1 - D) * MC2_EV:.7f}, rel {(1 - E) / (1 - D) - 1:+.2e} "
          f"(N={rw.N}, {time.time() - t0:.0f} s)", flush=True)
L = res["levels"]
res["fs_2p_eV"] = L["2p1/2"]["E_bind_eV"] - L["2p3/2"]["E_bind_eV"]
res["fs_2p_dirac_eV"] = L["2p1/2"]["dirac_bind_eV"] - L["2p3/2"]["dirac_bind_eV"]
print(f"fine structure 2p: walk {res['fs_2p_eV'] / 4.135667696e-15 / 1e9:.5f} GHz, Dirac {res['fs_2p_dirac_eV'] / 4.135667696e-15 / 1e9:.5f} GHz")
json.dump(res, open(os.path.join(OUT, "aggfull_fineref.json"), "w"), indent=1)
