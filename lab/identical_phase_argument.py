"""Why a density-only rule (HDW / Madelung without phase) cannot give the right fringes.

Two double slits that differ ONLY by an extra phase π on the right slit:
  * their screen patterns are opposite (central bright fringe vs central dark fringe);
  * but the balls start (z0 = 1) with the same density j_z and the same velocities j_x/(k|U|²)
    wherever there are balls, so a rule that only uses positions and velocities must give the same
    pattern for both -> it can be right for at most one of them.
  python lab/identical_phase_argument.py
"""
import numpy as np
from optics import sources, rs_field, flux, screen_reference, K, COARSE

z0 = 1.0
xg = np.linspace(-8, 8, 16001)
out = {}
for name, jones in (("same phase", ((1.0,), (1.0,))), ("right slit +π", ((1.0,), (-1.0,)))):
    s, w = sources(400, jones)
    U, Ux, Uz = rs_field(xg, z0, s, w)
    jx, jz = flux(U, Ux, Uz)
    out[name] = (jz, jx / (K * (np.abs(U) ** 2).sum(-1)), screen_reference("flux", jones=jones))
(j1, v1, r1), (j2, v2, r2) = out.values()
lit = j1 > 1e-3 * j1.max()                       # where the balls actually start
print(f"start density: max relative difference where j_z > 1e-3·max: {np.max(np.abs(j2[lit]-j1[lit])/j1[lit]):.2e}")
print(f"start velocity: max difference there: {np.max(np.abs(v2[lit]-v1[lit])):.2e}  (velocities up to {np.abs(v1[lit]).max():.2f})")
rel = np.abs(j2 - j1) / np.maximum(j1, 1e-300)
order = np.argsort(rel); cw = np.cumsum(j1[order]) / j1.sum()
for q in (0.5, 0.9, 0.99):
    print(f"  {100*q:.0f}% of the light starts where the two densities differ by less than {100*rel[order][np.searchsorted(cw, q)]:.2f}%")
dv = np.abs(v2 - v1); order = np.argsort(dv); cw = np.cumsum(j1[order]) / j1.sum()
for q in (0.5, 0.9, 0.99):
    print(f"  {100*q:.0f}% of the light starts where the two velocities differ by less than {dv[order][np.searchsorted(cw, q)]:.4f}")
c = np.corrcoef(r1.reshape(-1, COARSE).sum(1), r2.reshape(-1, COARSE).sum(1))[0, 1]
print(f"screen patterns: coarse correlation between the two = {c:.3f}; centre bin probabilities "
      f"{r1[118:122].sum():.4f} vs {r2[118:122].sum():.4f}")
