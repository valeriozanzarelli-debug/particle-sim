"""vl_campi.py - VIDEO ATOMO: fields used to move and colour the balls (from the coupled light+atom export).

Electron (48^3, a = 0.5 lambda_C):
  * rho_s = |psi|^2 smoothed with a gaussian sigma = 2 cells (1 lambda_C, FWHM 2.35 lambda_C): removes the
    period-4-cell checkerboard of the 4 lattice valleys (docs/BATTIMENTO_1S.md). Declared on screen.
  * the exported current j (cell-centred bilinear R+sR - L+sL, averaged over 16 ticks) is NOT the lattice face
    flux: its divergence does not match d rho/dt (measured below). So the flow is corrected, per frame interval:
        J = j_s - rho_s grad(chi),   div(rho_s grad chi) = div j_s + d rho_s/dt
    on the faces of the cells (j_s, rho_s: mean of the two frames). This is the smallest change of the
    velocity field (minimum of int |dJ|^2/rho) that makes the flow carry rho_s exactly: the rotation stays the
    rule's, the compression comes from continuity. Balls use the face fluxes inside each cell (RT0, like
    SURFISTI) and rho_s of the cell, linear in time: a uniform cloud in a cell stays uniform.
Light (64^3, h = 8 lambda_C):
  * only E is exported. B is rebuilt with Faraday dB/dt = -curl E (spectral curl, trapezoid in time, dt = 8),
    starting from B0 = x^ x E0 (the pulse moves along +x); the static residual this leaves where the pulse
    passed is removed (mean of B over frames 60-200, when only the atom's own light is left).

usage: python3 lab/vl_campi.py         -> scratch arrays (not in git) + lab/out/vl_campi.json
"""
import json
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, QUI)
from atomo_export_leggi import carica  # noqa: E402

SCR = os.environ.get("VL_SCR", "/tmp/claude-0/-home-claude-particle-sim/a13a264e-04a0-55c0-8c14-0e293e6cf6fd/scratchpad/vl")
SIG = float(os.environ.get("VL_SIG", 2.0))   # smoothing, cells (amendment 3: chosen from the checkerboard of the export)
A = 0.5            # atom cell, lambda_C
H = 8.0            # light cell (downsampled), lambda_C
DT = 8.0           # time between frames, lambda_C/c
R0 = 40.0          # sphere where the atom's light balls are born (lambda_C)
T_ON = 240.0       # from here on the pulse has left the sphere R0 (flux through R0 measured, see json)
RES_FRAMES = (60, 201)
FLUSSO = os.environ.get("VL_FLUSSO", "auto")   # amendment 1 (lab/out/vl_amend1.json): "continuita" | "corrente" | "auto"
# auto: use the exported current (minimally corrected) only if it is already almost a flux, i.e. the correction
# needed to obey continuity is <= 25 % of it on test intervals; otherwise continuity only (no fake rotation)


def kgrid(n, h):
    k = 2 * np.pi * np.fft.fftfreq(n, h)
    return np.meshgrid(k, k, k, indexing="ij")


def grad_spec(f, K):
    fh = np.fft.fftn(f)
    return np.stack([np.fft.ifftn(1j * K[i] * fh).real for i in range(3)])


def div_spec(v, K):
    return np.fft.ifftn(sum(1j * K[i] * np.fft.fftn(v[i]) for i in range(3))).real


def inv_lap(f, K):
    k2 = K[0] ** 2 + K[1] ** 2 + K[2] ** 2
    fh = np.fft.fftn(f)
    k2[0, 0, 0] = 1.0
    out = -fh / k2
    out[0, 0, 0] = 0.0
    return np.fft.ifftn(out).real


def curl_spec(F, K):
    Fh = [np.fft.fftn(F[i]) for i in range(3)]
    c = [1j * (K[1] * Fh[2] - K[2] * Fh[1]), 1j * (K[2] * Fh[0] - K[0] * Fh[2]), 1j * (K[0] * Fh[1] - K[1] * Fh[0])]
    return np.stack([np.fft.ifftn(x).real for x in c])


def _faces(c, axis):
    """cell-centred field -> interior face values (average of the two cells), zero on the box boundary"""
    n = c.shape[axis]
    sh = list(c.shape); sh[axis] = n + 1
    out = np.zeros(sh, c.dtype)
    sl = [slice(None)] * 3
    sl[axis] = slice(1, n)
    a = np.take(c, np.arange(n - 1), axis=axis)
    b = np.take(c, np.arange(1, n), axis=axis)
    out[tuple(sl)] = 0.5 * (a + b)
    return out


def _div(F):
    return sum(np.diff(F[ax], axis=ax) for ax in range(3)) / A


def _grad_faces(chi):
    out = []
    for ax in range(3):
        n = chi.shape[ax]
        sh = list(chi.shape); sh[ax] = n + 1
        g = np.zeros(sh)
        sl = [slice(None)] * 3
        sl[ax] = slice(1, n)
        g[tuple(sl)] = np.diff(chi, axis=ax) / A
        out.append(g)
    return out


def _cg(rf, f, x0, tol=1e-7, itmax=3000):
    """solve div(rf grad chi) = f (Neumann) by conjugate gradients on -L; rf: face densities (3 arrays)"""
    def L(x):
        g = _grad_faces(x)
        return _div([rf[a] * g[a] for a in range(3)])
    diag = sum(np.take(rf[a], np.arange(rf[a].shape[a] - 1), axis=a) + np.take(rf[a], np.arange(1, rf[a].shape[a]), axis=a)
               for a in range(3)) / A ** 2
    Minv = 1.0 / np.maximum(diag, 1e-30)
    b = -(f - f.mean())
    x = x0.copy()
    r = b - (-L(x))
    z = Minv * r
    p = z.copy()
    rz = (r * z).sum()
    nb = np.sqrt((b ** 2).sum()) + 1e-300
    for it in range(itmax):
        Ap = -L(p)
        al = rz / (p * Ap).sum()
        x += al * p
        r -= al * Ap
        if np.sqrt((r ** 2).sum()) < tol * nb:
            break
        z = Minv * r
        rz2 = (r * z).sum()
        p = z + (rz2 / rz) * p
        rz = rz2
    return x, it, float(np.sqrt((r ** 2).sum()) / nb)


def elettrone(d, info):
    n = 48
    nf = d.n
    rho_s = np.zeros((nf, n, n, n), np.float32)
    js = np.zeros((nf, 3, n, n, n), np.float32)
    for k in range(nf):
        sm = (lambda x: gaussian_filter(x, SIG, mode="constant")) if SIG > 0 else (lambda x: x)
        rho_s[k] = sm(d.rho(k))
        js[k] = np.stack([sm(c) for c in d.j(k)])
    Fx = np.zeros((nf - 1, n + 1, n, n), np.float32)
    Fy = np.zeros((nf - 1, n, n + 1, n), np.float32)
    Fz = np.zeros((nf - 1, n, n, n + 1), np.float32)
    chi = np.zeros((n, n, n))
    global FLUSSO
    if FLUSSO == "auto":
        test = []
        for k in (1, nf // 2):
            jm = 0.5 * (js[k].astype(np.float64) + js[k + 1])
            rm = 0.5 * (rho_s[k].astype(np.float64) + rho_s[k + 1])
            drdt = (rho_s[k + 1].astype(np.float64) - rho_s[k]) / DT
            jf = [_faces(jm[a], a) for a in range(3)]
            rf = [np.maximum(_faces(rm, a), 1e-7 * rm.max()) for a in range(3)]
            for a in range(3):
                sl = [slice(None)] * 3; sl[a] = 0; rf[a][tuple(sl)] = 0; sl[a] = -1; rf[a][tuple(sl)] = 0
            c_, _, _ = _cg(rf, _div(jf) + drdt, np.zeros((n, n, n)))
            g = _grad_faces(c_)
            test.append(float(np.sqrt(sum(((rf[a] * g[a]) ** 2).sum() for a in range(3)) / max(sum((jf[a] ** 2).sum() for a in range(3)), 1e-300))))
        FLUSSO = "corrente" if max(test) <= 0.25 else "continuita"
        info["scelta_flusso_auto"] = {"correzione_relativa_sui_test": test, "scelta": FLUSSO}
    with open(os.path.join(SCR, "flusso.txt"), "w") as fh:
        fh.write(FLUSSO)
    with open(os.path.join(SCR, "sigma.txt"), "w") as fh:
        fh.write(str(SIG))
    miss, corr_share, its, resid = [], [], [], []
    for k in range(nf - 1):
        jm = 0.5 * (js[k].astype(np.float64) + js[k + 1])
        rm = 0.5 * (rho_s[k].astype(np.float64) + rho_s[k + 1])
        drdt = (rho_s[k + 1].astype(np.float64) - rho_s[k]) / DT
        jf = [_faces(jm[a], a) for a in range(3)]
        if FLUSSO == "continuita":   # amendment 1: the cell-centred j is not a lattice flux; keep only continuity
            jfx = jf
            jf = [np.zeros_like(x) for x in jf]
        rf = [np.maximum(_faces(rm, a), 1e-7 * rm.max()) for a in range(3)]
        for a in range(3):   # no flux through the box boundary
            sl = [slice(None)] * 3; sl[a] = 0; rf[a][tuple(sl)] = 0; sl[a] = -1; rf[a][tuple(sl)] = 0
        divj = _div(jf)
        divj_exp = _div(jfx) if FLUSSO == "continuita" else divj
        if k in (10, 25, 60, 150):
            miss.append({"intervallo": k, "sum_abs_div_j": float(np.abs(divj_exp).sum()), "sum_abs_drho_dt": float(np.abs(drdt).sum()),
                         "corr(drho/dt, -div j)": float(np.corrcoef(drdt.ravel(), -divj_exp.ravel())[0, 1])})
        f = divj + drdt
        chi, it, res = _cg(rf, f, chi)
        g = _grad_faces(chi)
        F = [jf[a] - rf[a] * g[a] for a in range(3)]
        its.append(it); resid.append(res)
        if FLUSSO != "continuita":
            num = sum(((rf[a] * g[a]) ** 2).sum() for a in range(3))
            den = sum((jf[a] ** 2).sum() for a in range(3))
            corr_share.append(float(np.sqrt(num / den)))
        Fx[k], Fy[k], Fz[k] = F
        if k == 100:
            cont = _div(F) + drdt
            cont_rel = float(np.abs(cont).sum() / np.abs(drdt).sum())
            xa = (np.arange(n) - n / 2 + 0.5) * A
            X, Y, Z = np.meshgrid(xa, xa, xa, indexing="ij")
            Jc = np.stack([0.5 * (F[0][1:] + F[0][:-1]), 0.5 * (F[1][:, 1:] + F[1][:, :-1]), 0.5 * (F[2][:, :, 1:] + F[2][:, :, :-1])])
            rp = np.hypot(X, Y).clip(1e-9)
            jphi = (-Y * Jc[0] + X * Jc[1]) / rp
            jmag = np.sqrt((Jc ** 2).sum(0))
            vm = jmag / np.maximum(rm, 1e-12)
            w = rm / rm.sum()
            rot = {"quota_azimutale": float(np.abs(jphi).sum() / jmag.sum()), "quota_azimutale_con_segno": float(jphi.sum() / jmag.sum()),
                   "velocita_media_pesata_c": float((w * vm).sum()), "continuita_residuo_rel": cont_rel}
    def valley_amp(f):
        Fh = np.abs(np.fft.fftn(f))
        idx = [(12, 12, 0), (12, 0, 12), (0, 12, 12), (12, 36, 0), (36, 12, 0), (36, 36, 0)]
        return float(sum(Fh[i] for i in idx)) / float(Fh[0, 0, 0])
    info["elettrone"] = {
        "sigma_celle": SIG, "sigma_lambdaC": SIG * A,
        "continuita_corrente_esportata (facce, smussata, media dell'intervallo)": miss,
        "flusso": FLUSSO,
        "correzione_rho_grad_chi_su_j_mediana": float(np.median(corr_share)) if corr_share else None,
        "cg_iterazioni_max": int(max(its)), "cg_residuo_max": float(max(resid)),
        "fotogramma100": rot,
        "scacchiera_periodo4_rel_grezza": valley_amp(d.rho(100)), "scacchiera_periodo4_rel_smussata": valley_amp(rho_s[100]),
    }
    info["elettrone"]["attenuazione_scacchiera"] = info["elettrone"]["scacchiera_periodo4_rel_smussata"] / info["elettrone"]["scacchiera_periodo4_rel_grezza"]
    np.save(os.path.join(SCR, "rho_s.npy"), rho_s)
    np.save(os.path.join(SCR, "Fx.npy"), Fx)
    np.save(os.path.join(SCR, "Fy.npy"), Fy)
    np.save(os.path.join(SCR, "Fz.npy"), Fz)


def pexc_campo(d):
    """field-level estimator of p_exc using only the exported (dominant-component) phase: 1 - |sum sqrt(rho0 rho) e^{i dphi}|^2"""
    w1 = d.man["omega_1s"]
    r0, f0 = d.rho(0), d.fase(0)
    out = []
    for k in range(d.n):
        t = d.fotogrammi["t"][k]
        dphi = d.fase(k) - f0 - (1 - w1) * t
        Aov = (np.sqrt(d.rho(k) * r0) * np.exp(1j * dphi)).sum() * A ** 3
        out.append(1 - abs(Aov) ** 2)
    return np.array(out)


def luce(d, info):
    n = 64
    K = kgrid(n, H)
    nf = d.n
    E = np.stack([d.E(k) for k in range(nf)]).astype(np.float32)
    B = np.zeros_like(E)
    B[0, 1] = -E[0, 2]
    B[0, 2] = E[0, 1]
    Cprev = curl_spec(E[0], K)
    for k in range(1, nf):
        C = curl_spec(E[k], K)
        B[k] = B[k - 1] - 0.5 * DT * (Cprev + C)
        Cprev = C
    Bres = B[RES_FRAMES[0]:RES_FRAMES[1]].mean(0)
    B -= Bres[None]
    X = (2 * np.arange(n) + 0.5 - n) * 4.0
    XX, YY, ZZ = np.meshgrid(X, X, X, indexing="ij")
    R = np.sqrt(XX ** 2 + YY ** 2 + ZZ ** 2)
    rh = np.stack([XX, YY, ZZ]) / np.maximum(R, 1e-9)
    inner = (np.abs(XX) < 200) & (np.abs(YY) < 200) & (np.abs(ZZ) < 200)
    m = inner & (R > 60)
    cors = []
    for k in range(60, nf):
        rxE = np.cross(rh, E[k], axis=0)
        cors.append(float(np.corrcoef(B[k][:, m].ravel(), rxE[:, m].ravel())[0, 1]))
    pw = []
    for k in (5, 10, 20, 30):
        mm = inner
        pw.append(float(np.corrcoef(B[k, 1][mm], -E[k, 2][mm])[0, 1]))
    # flux through the sphere R0
    M = 2000
    i = np.arange(M) + 0.5
    th = np.arccos(1 - 2 * i / M)
    ph = np.pi * (1 + 5 ** 0.5) * i
    nh = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    c = (R0 * nh / 4.0 + 64 - 0.5) / 2.0
    P = []
    for k in range(nf):
        S = np.cross(E[k], B[k], axis=0)
        Sn = sum(map_coordinates(S[a], c, order=1) * nh[a] for a in range(3))
        P.append(float(Sn.mean() * 4 * np.pi * R0 ** 2))
    P = np.array(P)
    e2 = d.man["e_eff"] ** 2
    t = d.fotogrammi["t"]
    on = t >= T_ON
    # energy out through R0 after T_ON (trapezoid), in m c^2 (E is e*E: divide by e_eff^2, Heaviside-Lorentz)
    Pon = P[on]
    U_out = float(np.sum(0.5 * (Pon[1:] + Pon[:-1])) * DT / e2)
    pex = d.fotogrammi["p_exc"]
    k_on = int(np.argmax(on))
    U_atom = float((pex[k_on] - pex[-1]) * d.man["omega_21"])
    U_pulse = float(0.5 * ((E[0] ** 2).sum() + (B[0] ** 2).sum()) * H ** 3 / e2)
    info["luce"] = {
        "B_onda_piana_corr_By_-Ez_fotogrammi_5_10_20_30": pw,
        "B_vs_rxE_R>60_corr_media_60_200": float(np.mean(cors)), "B_vs_rxE_min": float(np.min(cors)),
        "flusso_R0": {"R0": R0, "T_ON": T_ON, "P_fotogrammi": {int(k): P[k] for k in (0, 5, 10, 20, 25, 30, 35, 40, 60, 100, 150, 200)}},
        "energia_impulso_mc2": U_pulse,
        "energia_uscita_da_R0_dopo_T_ON_mc2": U_out,
        "energia_persa_dall_atomo_dp_exc_x_omega21_mc2": U_atom,
        "rapporto_uscita_su_persa": U_out / U_atom,
        "energia_riemessa_su_impulso": U_out / U_pulse,
    }
    np.save(os.path.join(SCR, "E.npy"), E)
    np.save(os.path.join(SCR, "B.npy"), B)
    np.save(os.path.join(SCR, "P_R0.npy"), P)


if __name__ == "__main__":
    os.makedirs(SCR, exist_ok=True)
    d = carica()
    from vl_comune import impronta
    imp = impronta()
    info = {"fonte": d.dir, "commit_dati": "addffd7 (al 13/10, 1s esatto)", "impronta_esportazione": imp[0], "manifest_mtime_utc": imp[1]}
    elettrone(d, info)
    pc = pexc_campo(d)
    np.save(os.path.join(SCR, "pexc_campo.npy"), pc)
    info["pexc_stimatore_fase_dominante_sul_reticolo"] = {int(k): [float(pc[k]), float(d.fotogrammi["p_exc"][k])] for k in (25, 50, 100, 150, 200)}
    luce(d, info)
    with open(os.path.join(QUI, "out", "vl_campi.json"), "w") as f:
        json.dump(info, f, indent=1)
    print(json.dumps(info, indent=1))
