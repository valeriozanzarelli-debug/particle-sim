"""fluid_f2.py - F2 production: GP vortex energy vs ln R, core constant, swirl energy vs missing balls x c^2.
Writes out/fluid_f2.json, out/fluid_f2.npz. Protocol in out/fluid_prereg.json (F2)."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fluid_core import gp_vortex, gp_energy_parts

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
r, f = gp_vortex(400.0, 0.02)
Rs = np.array([10, 20, 40, 80, 160.0])
P = [gp_energy_parts(r, f, R) for R in Rs]
E = np.array([p["total"] for p in P]); Es = np.array([p["swirl"] for p in P]); Nm = np.array([p["nmiss"] for p in P])
lnR = np.log(Rs)
a, b = np.polyfit(lnR, E, 1)
A = float(np.exp(b / a))
A_fetter = 1.464 * np.sqrt(2.0)
d1, d0 = np.polyfit(lnR, Es - Nm, 1)  # m = c = 1
# core-only missing balls (r < 3 xi) and other candidate masses, at each R
core = gp_energy_parts(r, f, 3.0)
res = dict(R=Rs.tolist(), E=E.tolist(), E_swirl=Es.tolist(), N_miss=Nm.tolist(),
           slope=float(a), slope_over_pi=float(a / np.pi), A=A, A_fetter_converted=float(A_fetter),
           A_rel_dev=float(A / A_fetter - 1), swirl_minus_missing_slope=float(d1),
           swirl_minus_missing_const=float(d0), N_miss_core_3xi=float(core["nmiss"]),
           E_over_c2=E.tolist(),
           exams=dict(E2a=bool(abs(a / np.pi - 1) < 0.01), E2b=bool(abs(A / A_fetter - 1) < 0.02),
                      E2c=bool(abs(d1) < 0.02 * np.pi)))
json.dump(res, open(os.path.join(OUT, "fluid_f2.json"), "w"), indent=1)
np.savez_compressed(os.path.join(OUT, "fluid_f2.npz"), r=r[r < 30], f=f[r < 30])
print(json.dumps(res, indent=1))
