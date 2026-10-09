"""kk_speed_kk.py - TASK 1: pure Kaluza-Klein check of a separate hop speed c_w along the hidden circle.

(a) the walk cannot tell c_w from the circle size: spectrum of (N, c_w) == spectrum of (N / c_w, 1);
    a charged packet in a twisted geometry (N, c_w, slide D) evolves exactly like (N / c_w, 1, D / c_w).
(b) rest energy of n = 1 vs c_w on the walk: sqrt(theta^2 + (c_w n / R)^2).
(c) KK numbers: with the gauge coupling fixed by gravity (KK reduction) the charge/mass ratio of
    every KK mode stays q^2 / (G m^2) <= 4 (Gaussian units) for ANY c_w; electron: 4.17e42.
Seeded, deterministic. Writes out/kk_speed_kk.json. Run from lab/.
"""
import json, time
import numpy as np
import kk_speed_core as S
import kk_core as K

t0 = time.time()
out = {}
th = 0.02

# ---- (a1) spectrum equivalence ----
rows = []
for (N, cw, Nref) in ((32, 0.5, 64), (16, 0.25, 64), (32, 0.25, 128)):
    for n in (0, 1, 2, 3):
        E = S.rest_energy(N, th, n, cw)
        Eref = S.rest_energy(Nref, th, n, 1.0)
        rows.append(dict(N=N, cw=cw, N_equiv=Nref, n=n, E=E, E_equiv=Eref, diff=abs(E - Eref)))
        print(f"(a1) N={N:3d} c_w={cw:4.2f}  n={n}: E = {E:.12f}   vs  N={Nref:3d} c_w=1: {Eref:.12f}   |diff| = {abs(E-Eref):.1e}")
out["a1_spectrum_equivalence"] = rows

# ---- (a2) a charged packet in a twisted geometry: (N, cw, D) vs (N/cw, 1, D/cw) ----
rng = np.random.default_rng(K.SEED)
nx = ny = 48
D0 = rng.normal(size=(nx, ny)); D1 = rng.normal(size=(nx, ny))
res = []
for (N, cw) in ((16, 0.5), (12, 0.25)):
    Nref = int(round(N / cw))
    for n in (1, -1, 2):
        pa = S.packet_cw(nx, ny, N, 0.2, n, cw=cw, k0=(0.4, 0.1), sigma=4)
        pb = S.packet_cw(nx, ny, Nref, 0.2, n, cw=1.0, k0=(0.4, 0.1), sigma=4)
        Wa = S.KKWalkCw(nx, ny, N, 0.2, cw=cw, D={0: 0.7 * D0, 1: 0.7 * D1})
        Wb = S.KKWalkCw(nx, ny, Nref, 0.2, cw=1.0, D={0: 0.7 * D0 / cw, 1: 0.7 * D1 / cw})
        for t in range(120):
            pa = Wa.tick(pa); pb = Wb.tick(pb)
        da, db = K.density(pa), K.density(pb)
        err = float(np.max(np.abs(da - db)) / np.max(db))
        res.append(dict(N=N, cw=cw, N_equiv=Nref, n=n, max_rel_density_diff=err))
        print(f"(a2) twisted geometry, n={n:+d}: (N={N}, c_w={cw}, D) vs (N={Nref}, c_w=1, D/c_w): max density diff / max = {err:.1e}")
out["a2_twisted_equivalence"] = res

# ---- (b) rest energy of n = 1 vs c_w ----
N = 32; R = N / (2 * np.pi)
rb = []
for cw in (1.0, 0.5, 0.25, 0.1, 0.05, 0.02):
    E = S.rest_energy(N, th, 1, cw)
    pred = np.sqrt(th**2 + (cw / R)**2)
    rb.append(dict(cw=cw, E1=E, pred=pred, rel_err=(E - pred) / pred))
    print(f"(b) N={N} theta={th}: c_w={cw:5.2f}  E(n=1) = {E:.6f}   sqrt(theta^2+(c_w/R)^2) = {pred:.6f}  rel.err {(E-pred)/pred:+.1e}")
out["b_rest_energy_vs_cw"] = rb

# ---- (c) KK numbers in SI ----
hbar, c, G, me, QE = 1.054571817e-34, 299792458.0, 6.67430e-11, 9.1093837015e-31, 1.602176634e-19
eps0 = 8.8541878128e-12
alpha = 7.2973525693e-3
lP = np.sqrt(hbar * G / c**3)
R_KK = 2 * lP / np.sqrt(alpha)                 # coordinate radius giving alpha with c_w = c
electron_ratio = QE**2 / (4 * np.pi * eps0 * G * me**2)
print("(c) electron: e^2 / (4 pi eps0 G m_e^2) = %.3e" % electron_ratio)
rc = []
for cw in (1.0, 1e-5, 1e-10, 1e-15, 1e-20, R_KK * me * c / hbar):
    Rp = R_KK / cw                              # proper radius (c_w is the radion: g_ww = (c / c_w)^2)
    m1 = hbar / (Rp * c)                        # kg  (= hbar c_w / (R c^2))
    alpha_KK = 4 * lP**2 / Rp**2                # KK coupling fixed by the proper radius
    q2_over_Gm2 = 4 * np.pi * alpha_KK * hbar * c / (G * m1**2) / (4 * np.pi)   # Gaussian: q^2/(G m^2)
    rc.append(dict(cw=cw, proper_radius_m=Rp, m1_over_me=m1 / me, alpha_KK=alpha_KK, q2_over_G_m2=q2_over_Gm2))
    print(f"(c) c_w/c = {cw:8.1e}: proper R = {Rp:.2e} m, m(n=1) = {m1/me:.2e} m_e, alpha_KK = {alpha_KK:.2e}, q^2/(G m^2) = {q2_over_Gm2:.3f}")
out["c_kk_numbers"] = dict(electron_q2_over_Gm2=electron_ratio, R_KK_m=R_KK, rows=rc,
                           note="Gaussian units; Heaviside-Lorentz q^2/(G m^2) = 16 pi = 50.3")
out["time_s"] = time.time() - t0
json.dump(out, open("out/kk_speed_kk.json", "w"), indent=1)
print("time %.1f s" % out["time_s"])
