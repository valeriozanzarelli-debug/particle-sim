"""trap_null.py - Parte B1: luce chiusa nel vuoto LINEARE (= regola delle palline-onda = Maxwell).

Configurazioni:
  2D  anello che gira (fascio circolante): E_z = exp(-(r-R)^2/w^2) e^{i m phi}, m = round(2 pi R / lambda)
  3D  "fotone toroidale" (tipo Williamson-van der Mark): stesso anello in un tubo gaussiano, campo
      polarizzato lungo z e reso trasverso (div E = 0)
  3D  hopfione di Ranada (campo di luce annodato, ogni coppia di linee di campo e' concatenata)
Misure: frazione di energia rimasta nella regione iniziale in funzione del tempo; tempo di fuga
t_1/2 (meta' energia uscita) in funzione di R/lambda e w/lambda; per l'hopfione l'elicita'
(numero topologico) nel tempo.  Controllo "a palline": formula di Poisson (somma di ondine
ritardate) contro il propagatore spettrale.

uso:  timeout 590 python lab/trap_null.py 2d   (poi: torus, hopf)
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from trap_core import (grid2d, grid3d, tm2d_positive_frequency, tm2d_evolve,
                       transverse_project, maxwell3d_evolve, cross_khat, helicity3d,
                       kirchhoff_poisson_2d, TWO_PI)

OUT = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261009)
T0 = time.time()


def log(*a):
    print(f"[{time.time()-T0:6.1f}s]", *a, flush=True)


def t_half(ts, fr, level=0.5):
    """primo tempo in cui la frazione scende sotto 'level' (interpolazione lineare)"""
    for i in range(1, len(ts)):
        if fr[i] < level <= fr[i - 1]:
            return ts[i - 1] + (fr[i - 1] - level) * (ts[i] - ts[i - 1]) / (fr[i - 1] - fr[i])
    return float("nan")


PART = sys.argv[1] if len(sys.argv) > 1 else "2d"
res = {"2d_ring": [], "3d_torus": [], "hopfion": {}}

if PART == "2d":
    # ============================================================ 2D: anello circolante
    log("2D anello circolante (modo TM), propagatore esatto")
    cases = [(1, 0.5), (2, 0.5), (2, 1.0), (4, 0.5), (4, 1.0), (4, 2.0), (8, 1.0), (8, 2.0), (16, 1.0)]
    for (R, w) in cases:
        L = max(32.0, 5 * R + 6 * w + 8)              # nessun rientro dal bordo periodico prima di t_max
        dx = 1 / 10
        n = int(2 ** np.ceil(np.log2(L / dx)))
        L = n * dx
        X, Y, KX, KY = grid2d(n, L)
        r = np.hypot(X, Y)
        m = int(round(TWO_PI * R))
        psi = np.exp(-((r - R) / w) ** 2) * np.exp(1j * m * np.arctan2(Y, X))
        Ek, Bxk, Byk, K = tm2d_positive_frequency(psi, KX, KY)
        tmax = 2.0 * R + 4 * w
        ts = np.linspace(0, tmax, 49)
        fa, fd = [], []
        annulus = np.abs(r - R) < 2 * w
        disk = r < 2 * R
        E0 = None
        for t in ts:
            Ez, Bx, By = tm2d_evolve(Ek, Bxk, Byk, K, t)
            e = np.abs(Ez) ** 2 + np.abs(Bx) ** 2 + np.abs(By) ** 2   # media sul ciclo (x2)
            tot = e.sum()
            E0 = tot if E0 is None else E0
            fa.append(e[annulus].sum() / tot)
            fd.append(e[disk].sum() / tot)
        cons = tot / E0 - 1
        th_a, th_d = t_half(ts, fa), t_half(ts, fd)
        pred_a = np.sqrt(2 * R * 2 * w + (2 * w) ** 2)   # raggio tangente rettilineo esce da |r-R|<2w
        pred_d = np.sqrt(3) * R                            # ... ed esce dal disco r<2R
        res["2d_ring"].append(dict(R=R, w=w, m=m, n=n, L=L, t_half_annulus=th_a, t_half_disk=th_d,
                                   ray_pred_annulus=pred_a, ray_pred_disk=pred_d, energy_drift=cons,
                                   frac_annulus_end=fa[-1], frac_disk_end=fd[-1]))
        log(f"  R={R:5.1f} w={w:3.1f} m={m:4d}  t1/2(anello)={th_a:6.2f} [raggio {pred_a:5.2f}]  "
            f"t1/2(disco 2R)={th_d:6.2f} [raggio {pred_d:5.2f}]  frazione finale disco={fd[-1]:.3f}  "
            f"deriva energia={cons:.1e}")
        if (R, w) == (4, 1.0):
            res["2d_ring_curve"] = dict(t=ts.tolist(), annulus=fa, disk=fd)
            # ---- controllo a palline: somma di ondine ritardate (Poisson) contro lo spettrale
            u0 = np.fft.ifft2(Ek).real
            u1 = np.fft.ifft2(-1j * K * Ek).real
            tchk = 3.0
            uspec = np.fft.ifft2(Ek * np.exp(-1j * K * tchk)).real
            probes = [(4.0, 0.0), (0.0, 5.0), (-6.0, 1.0), (2.5, -2.5), (7.0, 3.0)]
            ub = kirchhoff_poisson_2d(u0, u1, X, Y, probes, tchk, nalpha=400, ntheta=1080)
            xs = X[:, 0]
            us = np.array([uspec[np.argmin(abs(xs - px)), np.argmin(abs(xs - py))] for px, py in probes])
            err = np.max(np.abs(ub - us)) / np.max(np.abs(u0))
            res["ball_check"] = dict(t=tchk, probes=probes, balls=ub.tolist(), spectral=us.tolist(),
                                     max_err_rel_to_max_u0=float(err))
            log(f"  controllo palline (Poisson) vs spettrale a t={tchk}: errore max {err:.2e} (rel. al max del campo)")
            np.savez_compressed(os.path.join(OUT, "trap_null_ring2d.npz"),
                                e0=np.abs(np.fft.ifft2(Ek))[::4, ::4] ** 2,
                                e_mid=np.abs(tm2d_evolve(Ek, Bxk, Byk, K, 4.0)[0])[::4, ::4] ** 2,
                                e_end=np.abs(tm2d_evolve(Ek, Bxk, Byk, K, 10.0)[0])[::4, ::4] ** 2, L=L)


if PART == "torus":
    # ============================================================ 3D: fotone toroidale
    log("3D fotone toroidale (anello in un tubo), Maxwell vettoriale esatto")
    n3, L3 = 128, 32.0
    (X, Y, Z), K = grid3d(n3, L3)
    dV = (L3 / n3) ** 3
    rho = np.hypot(X, Y)
    rr = np.sqrt(X**2 + Y**2 + Z**2)
    phi = np.arctan2(Y, X)
    for (R, w) in [(2, 0.75), (4, 0.75), (6, 0.75), (4, 1.5)]:
        m = int(round(TWO_PI * R))
        d2 = (rho - R) ** 2 + Z**2
        E = np.zeros((3,) + X.shape, complex)
        E[2] = np.exp(-d2 / w**2) * np.exp(1j * m * phi)       # polarizzazione z: trasversa alla corsa
        Ek = transverse_project(np.array([np.fft.fftn(c) for c in E]), K)
        Bk = cross_khat(K, Ek)                                   # frequenza positiva: B = khat x E
        kabs = np.sqrt(K[0]**2 + K[1]**2 + K[2]**2)
        tube = d2 < (2 * w) ** 2
        ball = rr < 2 * R
        ts = np.linspace(0, min(2.0 * R + 4 * w, 12.0), 21)
        ft, fb = [], []
        for t in ts:
            ph = np.exp(-1j * kabs * t)
            Et = np.array([np.fft.ifftn(c * ph) for c in Ek])
            Bt = np.array([np.fft.ifftn(c * ph) for c in Bk])
            e = (np.abs(Et) ** 2).sum(0) + (np.abs(Bt) ** 2).sum(0)
            tot = e.sum()
            ft.append(e[tube].sum() / tot)
            fb.append(e[ball].sum() / tot)
        th_t, th_b = t_half(ts, ft), t_half(ts, fb)
        pred_t = np.sqrt(2 * R * 2 * w + (2 * w) ** 2)
        res["3d_torus"].append(dict(R=R, w=w, m=m, t_half_tube=th_t, t_half_ball=th_b,
                                    ray_pred_tube=pred_t, ray_pred_ball=np.sqrt(3) * R,
                                    frac_tube_end=ft[-1], frac_ball_end=fb[-1], t_end=ts[-1]))
        log(f"  R={R} w={w} m={m}: t1/2(tubo)={th_t:5.2f} [raggio {pred_t:5.2f}]  "
            f"t1/2(sfera 2R)={th_b:5.2f} [raggio {np.sqrt(3)*R:5.2f}]  frazione nel tubo a t={ts[-1]:.1f}: {ft[-1]:.3f}")


if PART == "hopf":
    # ============================================================ 3D: hopfione di Ranada
    log("3D hopfione (campo annodato di Ranada / Irvine-Bouwmeester), scala a = 1")
    n3, L3 = 128, 44.0
    (X, Y, Z), K = grid3d(n3, L3)
    dV = (L3 / n3) ** 3
    r2 = X**2 + Y**2 + Z**2
    pref = 8.0 / (np.pi * (1 + r2) ** 3)
    Bh = pref * np.array([2 * (Y - X * Z), -2 * (X + Y * Z), X**2 + Y**2 - Z**2 - 1])
    Eh = pref * np.array([X**2 - Y**2 - Z**2 + 1, 2 * (X * Y - Z), 2 * (X * Z + Y)])
    Ek = np.array([np.fft.fftn(c) for c in Eh])
    Bk = np.array([np.fft.fftn(c) for c in Bh])
    divE = np.max(np.abs(sum(np.fft.ifftn(1j * K[i] * Ek[i]) for i in range(3))))
    divB = np.max(np.abs(sum(np.fft.ifftn(1j * K[i] * Bk[i]) for i in range(3))))
    EdotB = np.max(np.abs((Eh * Bh).sum(0))) / np.max((Eh**2).sum(0))
    # coerenza con Maxwell: dE/dt = curl B deve dare un campo trasverso consistente (controllo di segno)
    log(f"  max|div E|={divE:.1e} max|div B|={divB:.1e} max|E.B|/max|E|^2={EdotB:.1e} (|E|=|B| ovunque: "
        f"{np.max(np.abs((Eh**2).sum(0)-(Bh**2).sum(0))):.1e})")
    Fk = Ek + 1j * Bk
    rr = np.sqrt(r2)
    ts = np.linspace(0, 16, 17)
    hop = dict(t=ts.tolist(), f_r1=[], f_r2=[], f_r4=[], helicity_m=[], helicity_e=[], energy=[], rpeak=[])
    rb = np.linspace(0, 24, 97)
    for t in ts:
        Ft = maxwell3d_evolve(Fk, K, t)
        Ekt, Bkt = Ft.real * 0, None
        # separa E, B: E_k = (F_k + conj(F_{-k}))/2 nello spazio reale = Re(F), B = Im(F)
        Fx = np.array([np.fft.ifftn(c) for c in Ft])
        Et, Bt = Fx.real, Fx.imag
        e = 0.5 * ((Et**2).sum(0) + (Bt**2).sum(0))
        tot = e.sum() * dV
        hop["energy"].append(float(tot))
        hop["f_r1"].append(float(e[rr < 1].sum() * dV / tot))
        hop["f_r2"].append(float(e[rr < 2].sum() * dV / tot))
        hop["f_r4"].append(float(e[rr < 4].sum() * dV / tot))
        hm, he = helicity3d(np.array([np.fft.fftn(c) for c in Et]), np.array([np.fft.fftn(c) for c in Bt]),
                            K, dV, X.size)
        hop["helicity_m"].append(float(hm)); hop["helicity_e"].append(float(he))
        hist, _ = np.histogram(rr, bins=rb, weights=e)
        hop["rpeak"].append(float(0.5 * (rb[np.argmax(hist)] + rb[np.argmax(hist) + 1])))
    th = {k: t_half(ts, np.array(hop[k])) for k in ("f_r1", "f_r2", "f_r4")}
    hop["t_half"] = th
    hop["energy_drift"] = hop["energy"][-1] / hop["energy"][0] - 1
    hop["helicity_drift"] = (hop["helicity_m"][-1] + hop["helicity_e"][-1]) / (hop["helicity_m"][0] + hop["helicity_e"][0]) - 1
    res["hopfion"] = hop
    for i in range(0, len(ts), 2):
        log(f"  t={ts[i]:5.1f}  frazione r<1: {hop['f_r1'][i]:.4f}  r<2: {hop['f_r2'][i]:.4f}  r<4: {hop['f_r4'][i]:.4f}"
            f"  picco del guscio a r={hop['rpeak'][i]:5.2f}  elicita' (m,e)=({hop['helicity_m'][i]:.4f},{hop['helicity_e'][i]:.4f})")
    log(f"  t1/2: r<1 {th['f_r1']:.2f}, r<2 {th['f_r2']:.2f}, r<4 {th['f_r4']:.2f};  deriva energia {hop['energy_drift']:.1e}, "
        f"deriva elicita' {hop['helicity_drift']:.1e}")


with open(os.path.join(OUT, f"trap_null_{PART}.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
log(f"scritto lab/out/trap_null_{PART}.json")
