"""vl_campi.py - VIDEO ATOMO: fields used to move and colour the balls (from the coupled light+atom export).

Electron (48^3, a = 0.5 lambda_C):
  * rho_s = |psi|^2 smoothed with a gaussian sigma = 2 cells (1 lambda_C, FWHM 2.35 lambda_C): removes the
    period-4-cell checkerboard of the 4 lattice valleys (docs/BATTIMENTO_1S.md). Declared on screen.
  * the exported current j (cell-centred bilinear R+sR - L+sL, averaged over 16 ticks) is NOT the lattice face
    flux: its divergence does not match d rho/dt (measured below). So the flow is rebuilt:
        J = (divergence-free part of smoothed j)  +  grad(phi),   lap(phi) = -(rho_s[k+1] - rho_s[k]) / 8
    i.e. the rotation comes from the rule's current, the compression from continuity (L2-minimal correction,
    the same choice as SURFISTI S5). With rho linear in time inside a frame interval, the balls moved by
    v = J / rho_s carry rho_s exactly (up to the integrator).
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
SIG = 2.0          # smoothing of the electron fields, cells (1 lambda_C)
A = 0.5            # atom cell, lambda_C
H = 8.0            # light cell (downsampled), lambda_C
DT = 8.0           # time between frames, lambda_C/c
R0 = 40.0          # sphere where the atom's light balls are born (lambda_C)
T_ON = 240.0       # from here on the pulse has left the sphere R0 (flux through R0 measured, see json)
RES_FRAMES = (60, 201)


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


def elettrone(d, info):
    n = 48
    K = kgrid(n, A)
    nf = d.n
    rho_s = np.zeros((nf, n, n, n), np.float32)
    jrot = np.zeros((nf, 3, n, n, n), np.float32)
    share_grad, div_vs_dt = [], []
    for k in range(nf):
        rho_s[k] = gaussian_filter(d.rho(k), SIG, mode="constant")
        js = np.stack([gaussian_filter(c, SIG, mode="constant") for c in d.j(k)])
        g = grad_spec(inv_lap(div_spec(js, K), K), K)
        jrot[k] = js - g
        if k > 0:
            share_grad.append(float(np.sqrt((g ** 2).sum() / max((js ** 2).sum(), 1e-30))))
    G = np.zeros((nf - 1, 3, n, n, n), np.float32)
    for k in range(nf - 1):
        drdt = (rho_s[k + 1] - rho_s[k]) / DT
        G[k] = grad_spec(inv_lap(-drdt, K), K)
        if k >= 1:
            js_div = div_spec(jrot[k] + 0, K)  # ~0 by construction
            div_vs_dt.append(float(np.abs(js_div).sum() / max(np.abs(drdt).sum(), 1e-30)))
    # how badly the exported current misses continuity (raw, smoothed, central differences)
    miss = []
    for k in (10, 25, 60, 150):
        r0_, r1_ = gaussian_filter(d.rho(k), SIG, mode="constant"), gaussian_filter(d.rho(k + 1), SIG, mode="constant")
        js = np.stack([gaussian_filter(c, SIG, mode="constant") for c in d.j(k)])
        dv = div_spec(js, K)
        drdt = (r1_ - r0_) / DT
        miss.append({"frame": k, "sum_abs_div_j": float(np.abs(dv).sum()), "sum_abs_drho_dt": float(np.abs(drdt).sum()),
                     "corr(drho/dt, -div j)": float(np.corrcoef(drdt.ravel(), -dv.ravel())[0, 1])})
    # azimuthal share and speed of the rotation (frame 100)
    xa = (np.arange(n) - n / 2 + 0.5) * A
    X, Y, Z = np.meshgrid(xa, xa, xa, indexing="ij")
    k = 100
    jr = jrot[k]
    rp = np.hypot(X, Y).clip(1e-9)
    jphi = (-Y * jr[0] + X * jr[1]) / rp
    jm = np.sqrt((jr ** 2).sum(0))
    w = rho_s[k] / rho_s[k].sum()
    vm = jm / np.maximum(rho_s[k], 1e-12)
    # checkerboard: power at the period-4 valley wavevectors, raw vs smoothed
    def valley_amp(f):
        F = np.abs(np.fft.fftn(f))
        idx = [(12, 12, 0), (12, 0, 12), (0, 12, 12), (12, 36, 0), (36, 12, 0), (36, 36, 0)]
        return float(sum(F[i] for i in idx)), float(F[0, 0, 0])
    va_raw, v0_raw = valley_amp(d.rho(100))
    va_s, v0_s = valley_amp(rho_s[100])
    info["elettrone"] = {
        "sigma_celle": SIG, "sigma_lambdaC": SIG * A,
        "continuita_corrente_esportata": miss,
        "peso_parte_gradiente_della_j_esportata_mediana": float(np.median(share_grad)),
        "jrot_quota_azimutale_fotogramma100": float(np.abs(jphi).sum() / jm.sum()),
        "jrot_quota_azimutale_con_segno": float(jphi.sum() / jm.sum()),
        "velocita_media_pesata_rho_fotogramma100_c": float((w * vm).sum()),
        "scacchiera_periodo4_rel_grezza": va_raw / v0_raw, "scacchiera_periodo4_rel_smussata": va_s / v0_s,
        "attenuazione_scacchiera": (va_s / v0_s) / (va_raw / v0_raw),
    }
    np.save(os.path.join(SCR, "rho_s.npy"), rho_s)
    np.save(os.path.join(SCR, "jrot.npy"), jrot)
    np.save(os.path.join(SCR, "G.npy"), G)


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
    info = {"fonte": d.dir, "commit_dati": "addffd7"}
    elettrone(d, info)
    pc = pexc_campo(d)
    np.save(os.path.join(SCR, "pexc_campo.npy"), pc)
    info["pexc_stimatore_fase_dominante_sul_reticolo"] = {int(k): [float(pc[k]), float(d.fotogrammi["p_exc"][k])] for k in (25, 50, 100, 150, 200)}
    luce(d, info)
    with open(os.path.join(QUI, "out", "vl_campi.json"), "w") as f:
        json.dump(info, f, indent=1)
    print(json.dumps(info, indent=1))
