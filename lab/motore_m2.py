"""motore_m2.py - local time: exams C1 (causal-cone freeze, bit-exact), C2 (lazy clock: catch-up with W(k)^m),
C3 (tail of the learned 1-D agglomerated rules outside the cone, descriptive) + the Berger-Oliger saving computed on
the R3 world of RETICOLO_UNIVERSO (calculation). Pre-registration: lab/out/motore_prereg.json.
Usage: python lab/motore_m2.py  -> lab/out/motore_m2.json"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import motore_core as mc
import zz_core as zz

rng = np.random.default_rng(20261013)
res = {}


def frozen_run(field, tick, n, T, box, shrink=0):
    D = mc.chebyshev_dist_to_box(n, box)
    A = field.copy()
    Bf = field.copy()
    active = 0
    for t in range(T):
        A = tick(A, t)
        new = tick(Bf, t)
        m = D <= (T - 1 - t - shrink)
        Bf = np.where(m[..., None], new, Bf)
        active += int(m.sum())
    v = (slice(box[0], box[1]), slice(box[2], box[3]))
    err = float(np.abs(A[v] - Bf[v]).max() / np.abs(A[v]).max())
    return err, active / (T * n * n)


# ---------------- C1 ----------------
n, T, W = 192, 64, 48
box = (n // 2 - W // 2, n // 2 + W // 2, n // 2 - W // 2, n // 2 + W // 2)
x = np.arange(n) - n / 2 + 0.5
X, Y = np.meshgrid(x, x, indexing="ij")
# light: three packets + a smooth random background filling the box (so a too-narrow cone must fail)
F = mc.light_packet(n, 12, 30, 10, 10, center=(50, 60)) + mc.light_packet(n, 8, 200, 6, 6, center=(150, 120)) \
    + mc.light_packet(n, 16, 100, 12, 12, center=(96, 30))
noise = np.fft.ifft2(np.fft.fft2(rng.normal(size=(n, n, 3)) + 1j * rng.normal(size=(n, n, 3)), axes=(0, 1)) *
                     (np.exp(-(np.fft.fftfreq(n)[:, None] ** 2 + np.fft.fftfreq(n)[None, :] ** 2) / 0.01))[..., None], axes=(0, 1))
F = F + 0.3 * noise / np.abs(noise).max()
t0 = time.time()
eL, fracL = frozen_run(F, mc.light_tick_real, n, T, box)
eL_s, _ = frozen_run(F, mc.light_tick_real, n, T, box, shrink=1)
# Dirac, Z = 60 nucleus in the window, cell 0.3 lambda_C
th = 0.3
r = np.sqrt(X ** 2 + Y ** 2)
V = -60 * zz.ALPHA / np.maximum(r, 0.5) / th          # potential per unit time in cell units: V(r = x th) ; phase V*dt with dt = th
phase = np.exp(-1j * (-60 * zz.ALPHA / (np.maximum(r, 0.5) * th)) * th)
P = (rng.normal(size=(n, n, 4)) + 1j * rng.normal(size=(n, n, 4))) * np.exp(-r ** 2 / (2 * 40 ** 2))[..., None]
tickD = lambda psi, t: mc.dirac_tick_real(psi, t, th, phase)
eD, fracD = frozen_run(P, tickD, n, T, box)
eD_s, _ = frozen_run(P, tickD, n, T, box, shrink=1)
frac_an = sum(min(n, W + 2 * (T - 1 - t)) ** 2 for t in range(T)) / (T * n * n)
res["C1"] = dict(n=n, T=T, window=W, light_err=eL, light_err_cone_minus1=eL_s, dirac_err=eD, dirac_err_cone_minus1=eD_s,
                 active_fraction_measured=fracL, active_fraction_analytic=frac_an,
                 verdict="PASS" if (eL <= 1e-13 and eD <= 1e-13 and eL_s > 1e-6 and eD_s > 1e-6 and abs(fracL - frac_an) < 1e-12) else "FAIL",
                 wall=time.time() - t0)
print("C1", res["C1"], flush=True)

# saving as a function of the horizon (analytic, the same counting)
res["C1_saving_table"] = [dict(n=n, window=W, T=TT, active_fraction=sum(min(n, W + 2 * (TT - 1 - t)) ** 2 for t in range(TT)) / (TT * n * n))
                          for TT in (8, 16, 32, 64, 128)]

# ---------------- C2 ----------------
t0 = time.time()
n2, m = 128, 1024
F2 = mc.light_packet(n2, 10, 30, 12, 12) + 0.2 * mc.light_packet(n2, 7, 120, 8, 8, center=(30, 90))
A = F2.copy()
for t in range(m):
    A = mc.light_tick_real(A, t)
t_tick = time.time() - t0
t1 = time.time()
Bk = mc.apply_k(mc.matpow(mc.lattice_W(n2), m // 2), F2)
t_jump = time.time() - t1
e2 = float(np.abs(A - Bk).max() / np.abs(A).max())
res["C2"] = dict(n=n2, ticks=m, rel_err=e2, wall_tick_by_tick=t_tick, wall_jump=t_jump,
                 verdict="PASS" if e2 <= 1e-11 else "FAIL")
print("C2", res["C2"], flush=True)

# ---------------- C3 ----------------
th1, lv = mc.load_levels_1d()
c3 = []
for B in (16, 64, 192):
    rule = lv[B]
    w = int(np.abs(rule.offs).max())
    nn = 1024
    psi = np.zeros((1, 2, nn), complex)
    psi[0, 0, nn // 2] = 1.0
    psi[0, 1, nn // 2] = 1.0
    psi /= np.sqrt(2)
    xs = np.arange(nn) - nn // 2
    for s in (1, 4, 16):
        q = rule.step(psi, s)[0]
        rho = np.sum(np.abs(q) ** 2, axis=0)
        out_phys = float(rho[np.abs(xs) > s].sum() / rho.sum())
        out_sten = float(rho[np.abs(xs) > w * s].sum() / rho.sum())
        amp_far = float(np.sqrt(rho[np.abs(xs) > 2 * w * s].max() / rho.max())) if (2 * w * s) < nn // 2 else None
        c3.append(dict(B=B, w=w, steps=s, outside_physical_cone=out_phys, outside_stencil_cone=out_sten, amp_beyond_2ws=amp_far))
        print("C3", c3[-1], flush=True)
res["C3"] = c3

# ---------------- Berger-Oliger saving on the R3 world (calculation from RETICOLO_UNIVERSO.md section 5) ----------------
cells = [86272, 50176, 35840, 30720] + [25600] * 15 + [30464, 100096]
glob = sum(cells)
local = sum(c / 2 ** k for k, c in enumerate(cells))
res["berger_oliger_R3"] = dict(cells_per_level=cells, global_step_cost=glob, local_step_cost=local, saving=glob / local)
print("BO", res["berger_oliger_R3"]["saving"])
json.dump(res, open(os.path.join(mc.OUT, "motore_m2.json"), "w"), indent=1)
