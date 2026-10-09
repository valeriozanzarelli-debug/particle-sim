"""trap_dirac.py - diagnostica (NON la regola del surfista relativistico): Zitterbewegung e velocita'
di Bohm-Dirac per un pacchetto d'onda dell'equazione di Dirac libera.

Unita': hbar = c = m = 1  (lunghezza = lambdabar_C = hbar/mc = 3.86e-13 m, tempo = hbar/mc^2 = 1.29e-21 s).
1D:  H = sigma_x p + sigma_z m        (alpha = sigma_x, beta = sigma_z)
2D:  H = sigma_x p_x + sigma_y p_y + sigma_z m
Evoluzione esatta nello spazio di Fourier: exp(-i H(k) t) = cos(E t) - i sin(E t) H(k)/E, E = sqrt(k^2+m^2).
Velocita' di Bohm-Dirac: v = psi^dag alpha psi / psi^dag psi   (|v| <= 1 per Cauchy-Schwarz)

Misure:
 (1) 1D, pacchetto con spinore (1, i)/sqrt2 (mescola energie + e -): <x>(t) oscilla a frequenza ~2m
     (= 2 m c^2/hbar) con ampiezza ~ 1/(2m) (= hbar/2mc); il pacchetto proiettato sulle sole energie
     positive NON oscilla.
 (2) 2D, pacchetto a energia positiva, spin su, fermo: la velocita' di Bohm-Dirac GIRA attorno al
     centro (corrente di spin di Gordon); traiettorie, velocita' (<= c), raggio dei giri.

uso: timeout 590 python lab/trap_dirac.py
"""
import json, os, time
import numpy as np

OUT = os.path.join(os.path.dirname(__file__), "out")
T0 = time.time()
m = 1.0
res = {}


def log(*a):
    print(f"[{time.time()-T0:6.1f}s]", *a, flush=True)


# =============================================================== (1) 1D Zitterbewegung
n, L = 4096, 400.0
x = (np.arange(n) - n // 2) * L / n
k = 2 * np.pi * np.fft.fftfreq(n, L / n)
E = np.sqrt(k**2 + m**2)


def evolve1d(phik, t):
    c, s = np.cos(E * t), np.sin(E * t) / E
    a, b = phik
    # H = [[m, k],[k, -m]]
    return np.array([c * a - 1j * s * (m * a + k * b), c * b - 1j * s * (k * a - m * b)])


def pos_proj(phik):
    a, b = phik
    # proiettore (1 + H/E)/2
    return 0.5 * np.array([a + (m * a + k * b) / E, b + (k * a - m * b) / E])


ts = np.arange(0, 40.0, 0.05)
out1 = {}
for sig in (5.0, 1.0):
    g = np.exp(-x**2 / (4 * sig**2))
    psi0 = np.array([g, 1j * g]) / np.sqrt(2)
    psi0 /= np.sqrt(np.sum(np.abs(psi0) ** 2) * (L / n))
    for label, phik in (("misto", np.fft.fft(psi0, axis=1)), ("solo_E_pos", pos_proj(np.fft.fft(psi0, axis=1)))):
        X, VB = [], 0.0
        norm = None
        for t in ts:
            ps = np.fft.ifft(evolve1d(phik, t), axis=1)
            rho = np.sum(np.abs(ps) ** 2, 0)
            norm = rho.sum() if norm is None else norm
            X.append(float(np.sum(x * rho) / rho.sum()))
            jb = 2 * np.real(np.conj(ps[0]) * ps[1])            # psi^dag sigma_x psi
            msk = rho > 1e-12 * rho.max()
            VB = max(VB, float(np.max(np.abs(jb[msk] / rho[msk]))))
        X = np.array(X)
        Xc = X - X.mean()
        spec = np.abs(np.fft.rfft(Xc * np.hanning(len(Xc))))
        fr = 2 * np.pi * np.fft.rfftfreq(len(Xc), ts[1] - ts[0])
        wpk = float(fr[np.argmax(spec[1:]) + 1])
        amp_early = float(0.5 * (X[: int(5 / 0.05)].max() - X[: int(5 / 0.05)].min()))
        amp_late = float(0.5 * (X[-int(5 / 0.05):].max() - X[-int(5 / 0.05):].min()))
        out1[f"sigma{sig}_{label}"] = dict(sigma=sig, omega_peak=wpk, amp_t0_5=amp_early, amp_t35_40=amp_late,
                                           max_bohm_speed=VB, x=X[::4].tolist(), t=ts[::4].tolist())
        log(f"  1D sigma = {sig} lambdabar_C, {label:10s}: frequenza dominante di <x> = {wpk:.4f} (2m = 2),"
            f" ampiezza 0<t<5: {amp_early:.4f} (hbar/2mc = 0.5), 35<t<40: {amp_late:.4f}; max |v_Bohm| = {VB:.6f}")
res["zb_1d"] = out1

# =============================================================== (2) 2D velocita' di Bohm-Dirac che gira
n2, L2 = 256, 48.0
g1 = (np.arange(n2) - n2 // 2) * L2 / n2
X2, Y2 = np.meshgrid(g1, g1, indexing="ij")
kk = 2 * np.pi * np.fft.fftfreq(n2, L2 / n2)
KX, KY = np.meshgrid(kk, kk, indexing="ij")
E2 = np.sqrt(KX**2 + KY**2 + m**2)


def evolve2d(phik, t):
    c, s = np.cos(E2 * t), np.sin(E2 * t) / E2
    a, b = phik
    # H = [[m, kx - i ky],[kx + i ky, -m]]
    Ha = m * a + (KX - 1j * KY) * b
    Hb = (KX + 1j * KY) * a - m * b
    return np.array([c * a - 1j * s * Ha, c * b - 1j * s * Hb])


def pos2d(phik):
    a, b = phik
    Ha = m * a + (KX - 1j * KY) * b
    Hb = (KX + 1j * KY) * a - m * b
    return 0.5 * np.array([a + Ha / E2, b + Hb / E2])


def bohm_v(ps):
    rho = np.sum(np.abs(ps) ** 2, 0)
    vx = 2 * np.real(np.conj(ps[0]) * ps[1]) / rho          # sigma_x
    vy = 2 * np.imag(np.conj(ps[0]) * ps[1]) / rho          # sigma_y
    return vx, vy, rho


out2 = {}
for sig in (3.0, 1.0, 0.3):
    g = np.exp(-(X2**2 + Y2**2) / (4 * sig**2))
    phik = pos2d(np.fft.fft2(np.array([g, 0 * g]), axes=(1, 2)))     # spin su, energia positiva
    ps = np.fft.ifft2(phik, axes=(1, 2))
    vx, vy, rho = bohm_v(ps)
    rr = np.hypot(X2, Y2)
    vphi = (-Y2 * vx + X2 * vy) / np.maximum(rr, 1e-12)
    vr = (X2 * vx + Y2 * vy) / np.maximum(rr, 1e-12)
    msk = rho > 1e-8 * rho.max()
    sp = np.hypot(vx, vy)
    # velocita' azimutale media pesata entro r < 2 sigma e a r = sigma
    ring = msk & (np.abs(rr - sig) < L2 / n2)
    vphi_sigma = float(np.mean(vphi[ring]))
    weighted = float(np.sum(rho * sp) / np.sum(rho))
    # traiettorie di Bohm (RK2, campo ricalcolato a ogni passo)
    from scipy.interpolate import RegularGridInterpolator as RGI
    pts = np.array([[sig * np.cos(a), sig * np.sin(a)] for a in np.linspace(0, 2 * np.pi, 6, endpoint=False)]
                   + [[0.5 * sig, 0.0], [2 * sig, 0.0]])
    dt, tmax = 0.02, 12.0
    ang0 = np.arctan2(pts[:, 1], pts[:, 0])
    ang = ang0.copy()
    prev = ang0.copy()
    rmax_traj = np.hypot(pts[:, 0], pts[:, 1])
    speeds = []
    for it in range(int(tmax / dt)):
        def vel(p, t):
            q = np.fft.ifft2(evolve2d(phik, t), axes=(1, 2))
            vx_, vy_, _ = bohm_v(q)
            fx = RGI((g1, g1), vx_, bounds_error=False, fill_value=0.0)
            fy = RGI((g1, g1), vy_, bounds_error=False, fill_value=0.0)
            return np.stack([fx(p), fy(p)], 1)
        t = it * dt
        v1 = vel(pts, t)
        v2 = vel(pts + dt * v1, t + dt)
        pts = pts + 0.5 * dt * (v1 + v2)
        speeds.append(np.hypot(*(0.5 * (v1 + v2)).T))
        a = np.arctan2(pts[:, 1], pts[:, 0])
        ang += np.angle(np.exp(1j * (a - prev)))
        prev = a
        rmax_traj = np.maximum(rmax_traj, np.hypot(pts[:, 0], pts[:, 1]))
        if sig < 1 and it > 150:   # pacchetto stretto: si allarga in fretta, basta un tratto
            break
    speeds = np.array(speeds)
    turns = (ang - ang0) / (2 * np.pi)
    out2[f"sigma{sig}"] = dict(sigma=sig, max_speed_field=float(sp[msk].max()), mean_speed_weighted=weighted,
                               vphi_at_sigma=vphi_sigma, vphi_pred_1_over_2m_sigma=1 / (2 * m * sig) * 1.0,
                               mean_vr_abs=float(np.sum(rho * np.abs(vr)) / np.sum(rho)),
                               traj_turns=turns.tolist(), traj_time=float((it + 1) * dt),
                               traj_max_speed=float(speeds.max()), traj_mean_speed=float(speeds.mean()),
                               traj_rmax=rmax_traj.tolist())
    log(f"  2D sigma = {sig}: max|v| nel campo = {sp[msk].max():.4f} (<= 1 = c), |v| medio pesato = {weighted:.4f},"
        f" v_phi a r = sigma: {vphi_sigma:.4f} (spin-corrente ~ 1/(2 m sigma) = {1/(2*sig):.4f})")
    log(f"     traiettorie per t = {(it+1)*dt:.1f}: giri = {np.round(turns, 3).tolist()}, velocita' max {speeds.max():.4f},"
        f" media {speeds.mean():.4f}; raggio max {np.round(rmax_traj, 2).tolist()}")
res["bohm_2d"] = out2

with open(os.path.join(OUT, "trap_dirac.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
log("scritto lab/out/trap_dirac.json")
