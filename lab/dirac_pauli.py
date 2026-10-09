"""dirac_pauli.py - low-speed limit of the Dirac surfer rule (2+1 D).

Same initial positions, three rules:
  D  Dirac surfer           v = psi^dag alpha psi / psi^dag psi, positive-energy packet, spinor (1,0)
  P  Pauli surfer           v = [Im(phi* grad phi) + (1/2) curl(rho z)] / (m rho), phi Schroedinger-evolved
                            (the non-relativistic surfer rule + the spin term of Gordon / Holland)
  S  Schroedinger surfer    v = Im(phi* grad phi) / (m rho)  (the old non-relativistic rule, no spin term)
Small parameter eps = 1/sigma = k0 (speed ~ eps c). Prediction: |x_D - x_P| = O(eps^2) relative to sigma
|x_D - x_S| = O(1) relative to sigma at times T = m sigma^2 / 2 (spin rotation).
Version 2 (after version 1 - eps = 0.2, 0.1, 0.05 - turned out to need ~40 min; only its eps = 0.2 point
was obtained: D-P 3.27e-2 sigma, D-S 0.348 sigma): eps = 0.2, 0.2/sqrt2, 0.1, expected ratio 2 per step.
Tolerances (v2, fixed before its run): each ratio in [1.5, 2.5]; eps = 0.1: D-P < 2e-2 sigma; D-S > 0.1 sigma.
Run: timeout 590 python lab/dirac_pauli.py
"""
import json, os, time
import numpy as np
from scipy.ndimage import map_coordinates, spline_filter
from dirac_core import Dirac, gaussian_packet, SEED, DT, INTERP_ORDER

t_start = time.time()
OUT = os.path.join(os.path.dirname(__file__), "out")
rng = np.random.default_rng(SEED)
verdict, res = [], {}


def check(name, ok, msg):
    verdict.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def interp(D, coef, pos):
    idx = [(p / D.dx + D.n // 2) for p in pos]
    return map_coordinates(coef, idx, order=INTERP_ORDER, mode="grid-wrap", prefilter=False)


def nr_fields(D, phik, t, spin):
    KX, KY = D.k
    ft = phik * np.exp(-0.5j * (KX**2 + KY**2) * t)
    ph = np.fft.ifft2(ft)
    dxp = np.fft.ifft2(1j * KX * ft)
    dyp = np.fft.ifft2(1j * KY * ft)
    rho = np.abs(ph)**2
    jx = np.imag(np.conj(ph) * dxp)
    jy = np.imag(np.conj(ph) * dyp)
    if spin:
        drx = 2 * np.real(np.conj(ph) * dxp); dry = 2 * np.real(np.conj(ph) * dyp)
        jx = jx + 0.5 * dry; jy = jy - 0.5 * drx
    return [spline_filter(f, order=INTERP_ORDER, mode="grid-wrap") for f in (rho, jx, jy)]


def run_nr(D, phi0, pos0, T, spin):
    phik = np.fft.fft2(phi0)
    pos = pos0.copy()

    def v(c, p):
        r = interp(D, c[0], p)
        return np.array([interp(D, c[1], p) / r, interp(D, c[2], p) / r])

    n = int(round(T / DT))
    c_t = nr_fields(D, phik, 0.0, spin)
    for s in range(n):
        t = s * DT
        c_h = nr_fields(D, phik, t + DT / 2, spin); c_1 = nr_fields(D, phik, t + DT, spin)
        k1 = v(c_t, pos); k2 = v(c_h, pos + DT / 2 * k1); k3 = v(c_h, pos + DT / 2 * k2); k4 = v(c_1, pos + DT * k3)
        pos = pos + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        c_t = c_1
    return pos


rows = []
for sig in (5.0, 5.0 * np.sqrt(2), 10.0):
    t0 = time.time()
    k0 = 1 / sig
    D = Dirac(96, 12 * sig, 2)
    psi = gaussian_packet(D, sig, [k0, 0.0], [1, 0])
    phi0 = psi[0].copy()                                  # non-relativistic wave = upper component
    psiD = D.project(psi, +1); psiD /= np.sqrt(D.norm(psiD))
    pos0 = D.birth(psiD, 300, rng)
    T = 0.5 * sig**2
    _, trD, vmax = D.run_surfers(psiD, pos0, T, record_every=int(round(T / DT)))
    xD = trD[-1]
    xP = run_nr(D, phi0, pos0, T, spin=True)
    xS = run_nr(D, phi0, pos0, T, spin=False)
    dDP = np.sqrt(np.mean(np.sum((xD - xP)**2, axis=0))) / sig
    dDS = np.sqrt(np.mean(np.sum((xD - xS)**2, axis=0))) / sig
    rows.append(dict(sigma=sig, k0=k0, T=T, rms_D_minus_P_over_sigma=dDP, rms_D_minus_S_over_sigma=dDS, vmax=vmax,
                     runtime=time.time() - t0))
    print(f"eps = 1/sigma = k0 = {k0:.3f}, T = {T:.1f}: RMS|x_D - x_Pauli|/sigma = {dDP:.3e}, "
          f"RMS|x_D - x_Schroedinger(no spin)|/sigma = {dDS:.3e}, max|v_D| = {vmax:.4f} c  [{time.time() - t0:.1f} s]", flush=True)
res["rows"] = rows
r1 = rows[0]["rms_D_minus_P_over_sigma"] / rows[1]["rms_D_minus_P_over_sigma"]
r2 = rows[1]["rms_D_minus_P_over_sigma"] / rows[2]["rms_D_minus_P_over_sigma"]
res["ratios"] = [r1, r2]
check("P1 Dirac -> Pauli at second order in v/c", 1.5 <= r1 <= 2.5 and 1.5 <= r2 <= 2.5, f"ratios {r1:.2f}, {r2:.2f} (expected 2 per factor sqrt2 in eps)")
check("P2 smallest eps: Dirac and Pauli surfers coincide", rows[2]["rms_D_minus_P_over_sigma"] < 2e-2,
      f"{rows[2]['rms_D_minus_P_over_sigma']:.2e} sigma")
check("P3 the spin term is needed (Schroedinger-only rule differs at O(1))", min(r["rms_D_minus_S_over_sigma"] for r in rows) > 0.1,
      ", ".join(f"{r['rms_D_minus_S_over_sigma']:.3f}" for r in rows) + " sigma")
res["verdict"] = verdict
res["runtime_s"] = time.time() - t_start
with open(os.path.join(OUT, "dirac_pauli.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
print(f"=== {sum(v for _, v in verdict)}/{len(verdict)} PASS, runtime {time.time() - t_start:.1f} s ===", flush=True)
