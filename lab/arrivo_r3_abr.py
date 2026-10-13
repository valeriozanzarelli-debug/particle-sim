"""arrivo_r3_abr.py - Tumulka's absorbing boundary rule (NR, exam/comparison), the rule's detector reflection
for stationary NR waves (E13), the flux CDF, and the prediction numbers for the proposed experiment (P).
Writes lab/out/arrivo_r3.json."""
import json
import os
import time
import numpy as np
import arrivo_abr as B
import arrivo_nr as A
from tunnel_core import _Msite, _modes

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
t0 = time.time()
res = {}

# E13: stationary reflection of the rule's cell detector (lattice, exact) vs (1 - v)/(1 + v)
e13 = []
a = 0.025
for v in (0.05, 0.1, 0.3):
    E = 1 / np.sqrt(1 - v**2)
    # lattice quasi-energy with the same group velocity is not needed: use omega = E (continuum energy)
    Mf = _Msite(E, a, a, 0.0)
    lam, vv, f = _modes(Mf)
    iin = int(np.argmax(f))
    vin = vv[:, iin] / np.sqrt(abs(f[iin]))
    vref = vv[:, 1 - iin] / np.sqrt(abs(f[1 - iin]))
    R_lat = abs(vin[1] / vref[1])**2
    # lattice group velocity at this quasi-energy: cos(w a) = cos th cos(k a)
    k = np.arccos(np.cos(E * a) / np.cos(a)) / a
    vg = np.cos(a) * np.sin(k * a) / np.sin(E * a)
    kNR = v * E   # momentum
    e13.append(dict(v=v, R_lattice=R_lat, R_dirac=(1 - v) / (1 + v), R_dirac_vg=(1 - vg) / (1 + vg),
                    R_schrod_kappa2m=((2 - kNR) / (2 + kNR))**2, rel=abs(R_lat / ((1 - vg) / (1 + vg)) - 1)))
res["E13"] = e13
print("E13", e13, flush=True)

# flux CDF at L = 10 and L = 100
tg = np.linspace(1e-4, 200, 400001)
for L in (10.0, 100.0):
    P = A.flux(L, tg)
    C = np.concatenate([[0], np.cumsum(0.5 * (P[1:] + P[:-1]) * np.diff(tg))])
    sel = np.arange(0, tg.size, 200)
    res["flux_L%g" % L] = dict(tau=tg[sel].tolist(), Pi=P[sel].tolist(), F=C[sel].tolist(),
                               min_Pi=float(P.min()))

# ABR, L = 10: kappa scan and kappa = 2 D (the rule's detector)
abr = {}
for kap, name in ((np.pi / 2, "pi/2"), (np.pi, "pi"), (2 * np.pi, "2pi"), (10 * np.pi, "10pi"),
                  (40.0, "2D_20"), (80.0, "2D_40"), (160.0, "2D_80")):
    t, r, n = B.cn_abr(10.0, kap, 100.0, h=0.005, dt=5e-4, rec_every=20)
    F = np.cumsum(r) * 5e-4 * 20
    abr[name] = dict(kappa=kap, tau=t.tolist(), F=F.tolist(), rho=r.tolist(), F20=float(np.interp(20, t, F)),
                     F100=float(F[-1]), norm_loss=float(1 - n[-1]), mean_trunc=float(np.sum(t * r) / np.sum(r)))
    print("ABR L=10 kappa", name, "F(20) %.4f F(100) %.4f" % (abr[name]["F20"], abr[name]["F100"]), flush=True)
# grid check for kappa = pi and 2D_80
for kap in (np.pi, 160.0):
    t, r, n = B.cn_abr(10.0, kap, 100.0, h=0.01, dt=1e-3, rec_every=10)
    F = np.cumsum(r) * 1e-3 * 10
    abr["check_h2_%g" % kap] = dict(F20=float(np.interp(20, t, F)), F100=float(F[-1]))
res["ABR_L10"] = abr
# L = 100, kappa = pi (Jozani-Tumulka 90 %)
t, r, n = B.cn_abr(100.0, np.pi, 200.0, h=0.01, dt=1e-3, rec_every=20)
F = np.cumsum(r) * 1e-3 * 20
res["ABR_L100_pi"] = dict(tau=t[::5].tolist(), F=F[::5].tolist(), F50=float(np.interp(50, t, F)),
                          F100=float(np.interp(100, t, F)), F200=float(F[-1]))
print("ABR L=100 kappa=pi F(50,100,200)", res["ABR_L100_pi"]["F50"], res["ABR_L100_pi"]["F100"], F[-1], flush=True)

# modal formula check against CN (kappa = 160, L = 10)
mod = B.modal_large_kappa(10.0, 160.0, 20000)
res["modal_check_160"] = dict(modal_F=B.modal_cdf(mod, [20, 50, 100]).tolist(),
                              CN_F=[float(np.interp(x, abr["2D_80"]["tau"], abr["2D_80"]["F"])) for x in (20, 50, 100)])

# PREDICTION for the proposed experiment: electron, d = 50 um, L = 100 d = 5 mm
HBAR, ME, C = 1.054571817e-34, 9.1093837015e-31, 299792458.0
lamC = HBAR / (ME * C)
d = 50e-6
Dphys = d / lamC
tu = ME * d**2 / HBAR
kap = 2 * Dphys
pred = dict(d_m=d, L_m=100 * d, D_over_lambdaC=Dphys, time_unit_s=tu, kappa_DD=kap,
            v_typ_m_s=np.pi * HBAR / (ME * d), v_over_c=np.pi * HBAR / (ME * d) / C,
            absorption_per_hit=4 * np.pi / kap)
mod = B.modal_large_kappa(100.0, kap, 400000)
pred["mean_s"] = mod["mean"] * tu
pred["median_s"] = B.modal_quantile(mod, 0.5) * tu
pred["q10_s"] = B.modal_quantile(mod, 0.1) * tu
pred["q90_s"] = B.modal_quantile(mod, 0.9) * tu
pred["F_at_1ms"] = float(B.modal_cdf(mod, [1e-3 / tu])[0])
pred["F_at_taumax_DD"] = float(B.modal_cdf(mod, [42.9])[0])
pred["F_at_1s"] = float(B.modal_cdf(mod, [1.0 / tu])[0])
pred["wsum_modes"] = mod["wsum"]
pred["flux_mean_s"] = None
# flux (ideal non-reflecting detector) and DD Bohm in SI for the same set-up
Pf = A.flux(100.0, tg)
pred["flux_mean_s_trunc200"] = float(np.trapezoid(tg * Pf, tg) / np.trapezoid(Pf, tg) * tu)
pred["DD_taumax_s_paper"] = 42.9 * tu
res["PRED"] = pred
print("PRED", pred, flush=True)
res["secs"] = time.time() - t0
json.dump(res, open(os.path.join(OUT, "arrivo_r3.json"), "w"), indent=0)
print("done %.0f s" % res["secs"])
