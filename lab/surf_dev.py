"""surf_dev.py - DEVELOPMENT checks (declared, not production): correctness of the movers on small random
problems, timing. No exam is judged here."""
import time
import numpy as np
from surf_core import *

rng = np.random.default_rng(1)
t0 = time.time()
n = 12
W = Walk3D(n, 0.3, V=0.3 * rng.random((n, n, n)))
psi = rng.normal(size=(4, n, n, n)) + 1j * rng.normal(size=(4, n, n, n))
psi /= np.linalg.norm(psi)
# 1) chi2 after many ticks of S4 from a random (very rough) state
N = 400000
r = density(psi)
cell, xi = born_sample(rng, r, N)
hist = []
Fs_store, rho_store = [], []
cell_s, xi_s = cell.copy(), xi.copy()
for t in range(30):
    r0 = density(psi)
    psi, subs = W.tick(psi)
    r1 = density(psi)
    F = [None] * 3
    for (j, A, B) in subs:
        F[j] = net_flux(j, A, B)
    ncr, st = flow_tick(cell, xi, r0, r1, F)
    Fs_store.append(F); rho_store.append((r0, r1))
    z = chi2_z(cell, r1)
    hist.append(z[0])
print("S4 random 12^3, 30 ticks: z =", np.round(hist[::3], 2), "max crossings", ncr.max(), "stuck", st)
# 2) reversibility: run back with -F
c2, x2 = cell.copy(), xi.copy()
for F, (r0, r1) in zip(Fs_store[::-1], rho_store[::-1]):
    flow_tick(c2, x2, r1, r0, [-f for f in F])
dpos = np.abs((c2 + x2) - (cell_s + xi_s))
dpos = np.minimum(dpos, n - dpos)
print("reversibility: max |x_back - x0| = %.2e, median %.2e" % (dpos.max(), np.median(dpos)))
# 3) Liouville: rho J conserved (finite differences)
M = 2000
h = 1e-7
c0, x0 = born_sample(rng, rho_store[0][0], M)
x0 = np.clip(x0, 0.01, 0.99)
pts = [(c0.copy(), x0.copy())]
for a_ in range(3):
    xx = x0.copy(); xx[:, a_] += h
    pts.append((c0.copy(), xx))
for F, (r0, r1) in zip(Fs_store, rho_store):
    for (c, x) in pts:
        flow_tick(c, x, r0, r1, F)
P = [(c + x) for (c, x) in pts]
Jm = np.stack([((P[k + 1] - P[0] + n / 2) % n - n / 2) / h for k in range(3)], axis=2)
J = np.abs(np.linalg.det(Jm))
rho_start = rho_store[0][0][tuple(c0.T)]
rho_end = rho_store[-1][1][tuple(pts[0][0].T)]
err = np.abs(rho_end * J / rho_start - 1)
print("Liouville rho J / rho0 - 1: median %.2e, 90%% %.2e, 99%% %.2e" % (np.median(err), np.quantile(err, 0.9), np.quantile(err, 0.99)))
print("time %.1f s" % (time.time() - t0))
