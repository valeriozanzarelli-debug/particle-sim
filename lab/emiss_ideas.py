"""emiss_ideas.py - numerical check of the most promising "geometric/visual" idea of docs/EMISSIONE.md:
"the light of the atom is a DONUT (energy density ~ sin^2 theta around the dipole axis), and on screen it should be
round, not cubic".  If the lattice light rule were visibly anisotropic (cubic), balls drawn from |F|^2 would show
a square donut and the idea "density -> balls" would lie.  INFORMATIVE test (not in the pre-registration v3).

Uses the light field of a driven-dipole run (P_R0 if present, else dev_light): energy density 1/2|F|^2 in the shell
1.5 lambda < r < (n/2 - sponge - 2) cells, binned in cos(theta) and phi.
Run: python lab/emiss_ideas.py -> lab/out/emiss_ideas.json
"""
import json, os
import numpy as np
from emiss_core import CK, OUT

name = "P_R0" if os.path.exists(os.path.join(CK, "run_P_R0.npz")) else "dev_light"
d = np.load(os.path.join(CK, f"run_{name}.npz"), allow_pickle=True)
F = d["F"][0]
n = F.shape[-1]
u = 0.5 * np.sum(np.abs(F)**2, 0)
q = np.arange(n) - n // 2
X, Y, Z = np.meshgrid(q, q, q, indexing="ij")
r = np.sqrt(X**2 + Y**2 + Z**2) + 1e-12
lam_cells = 20.0
m = (r > 1.5 * lam_cells) & (r < n // 2 - 12 - 2)
ct = (Z / r)[m]
ph = np.arctan2(Y, X)[m]
w = (u * r**2)[m]                       # r^2 u: constant along a ray for a 1/r wave
# polar profile
cb = np.linspace(-1, 1, 21)
idx = np.clip(np.digitize(ct, cb) - 1, 0, 19)
prof = np.array([w[idx == i].mean() for i in range(20)])
cc = 0.5 * (cb[1:] + cb[:-1])
M = np.vstack([1 - cc**2, np.ones_like(cc)]).T
(A, B), *_ = np.linalg.lstsq(M, prof, rcond=None)
fit = M @ np.array([A, B])
# azimuthal profile on the equator band |cos theta| < 0.3
eq = np.abs(ct) < 0.3
pb = np.linspace(-np.pi, np.pi, 33)
ip = np.clip(np.digitize(ph[eq], pb) - 1, 0, 31)
az = np.array([w[eq][ip == i].mean() for i in range(32)])
pc = 0.5 * (pb[1:] + pb[:-1])
c4 = 2 * np.mean(az * np.cos(4 * pc)) / np.mean(az)
res = dict(source=name, n=n, shell_cells=[1.5 * lam_cells, n // 2 - 14],
           polar_cos=cc.tolist(), polar=prof.tolist(), fit_sin2=float(A), fit_const=float(B),
           const_over_sin2=float(B / A), polar_rms_resid=float(np.sqrt(np.mean((prof - fit)**2)) / A),
           az_phi=pc.tolist(), az=az.tolist(), az_cos4_rel=float(c4), az_minmax=float(az.max() / az.min() - 1))
json.dump(res, open(os.path.join(OUT, "emiss_ideas.json"), "w"), indent=1)
print(f"{name}: r^2 u(theta) = {A:.4e} sin^2 + {B:.4e}  (const/sin^2 = {B / A:.4f}, rms residual {res['polar_rms_resid']:.4f}); "
      f"equator: cos4phi amplitude {c4:+.4f}, max/min - 1 = {res['az_minmax']:.4f}")
