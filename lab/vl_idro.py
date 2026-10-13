"""vl_idro.py - VIDEO IDROGENO: zoom 'dal grande al piccolo' on real hydrogen (Z = 1) at the Bohr scale.

Data: /mnt/project-files/ripresa/dati/idrogeno_bohr/ (docs/IDROGENO_BOHR.md). Pre-registration: lab/out/vl_idro_prereg.json.

    python3 lab/vl_idro.py campi      # continuity-corrected face fluxes of the far superposition (scratch)
    python3 lab/vl_idro.py palline    # far balls over 3 Lyman-alpha periods + still 1s balls (coarse, fine)
    python3 lab/vl_idro.py esami      # H1-H6 -> lab/out/vl_idro_esami.json
    python3 lab/vl_idro.py render A B # video frames A..B-1 (resumable)
    python3 lab/vl_idro.py video      # mp4s + stills

Parts: A far view (1s + 2p_z swinging, balls moved by the rule's current minimally corrected for continuity),
B cross-fade to the still 1s and zoom from 6 a0 to 0.22 a0 (cells a0/8 -> a0/16 -> a0/256), C the 16 fine ticks
on the proton (colour = |R|^2 - |L|^2, amplified and declared). Every pixel of the scene is a soft ball.
"""
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import gaussian_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vl_comune import SCR as SCR0, OUT  # noqa: E402

BASE = "/mnt/project-files/ripresa/dati/idrogeno_bohr"
DIR = os.path.join(BASE, "fase2")          # amendment 1: phase-2 data (manifest -> fase2, commit 056007a)
SCR = os.path.join(SCR0, "idro2")
MAN = json.load(open(os.path.join(BASE, "manifest.json")))
M2 = MAN["fase2"]
LO = MAN["lontano"]
NA, AA = LO["griglia"]["n"], LO["griglia"]["passo_a0"]        # 96, a0/8
NF = LO["fotogrammi"]                                         # 24
_DIP = np.loadtxt(os.path.join(DIR, "dipolo_lontano.txt"))
DT = float(_DIP[1, 0] - _DIP[0, 0])                           # 0.701804 (phase 2)
T = NF * DT                                                   # 16.8433 a.u.
PER_S = T * LO["periodo_s"] / LO["periodo_au"]                # seconds per period
NG, AG = MAN["zoom"]["griglia_grossa"]["n"], MAN["zoom"]["griglia_grossa"]["passo_a0"]   # 64, a0/16
NFI, AF = M2["zoom"]["rho_fine.f16"][1], M2["zoom"]["passo_a0"]      # 112, a0/256
NT = M2["zoom"]["rho_fine.f16"][0]                                    # 16
HB = NFI / 2 * AF                                                     # half side of the fine block, 0.21875 a0
NB = int(round(2 * HB / AG))                                          # 7 coarse cells
R_CUT = 0.21                                                          # drawn sphere
CHIR_FULL = 4e-3
N_A = N_G = N_F = 100000
SEED = 20261014
PER = 3                 # periods shown
VPP = 96                # video frames per period
NV_A = PER * VPP
NV_X = 48               # cross-fade
NV_Z = 192              # zoom
TPF = 12                # video frames per fine tick
NV_C = NT * TPF
NV = NV_A + NV_X + NV_Z + NV_C
W, H = 1920, 1080


def impronta():
    h = hashlib.sha256()
    h.update(open(os.path.join(BASE, "manifest.json"), "rb").read())
    h.update(open(os.path.join(DIR, "dipolo_lontano.txt"), "rb").read())
    return h.hexdigest()[:10]


def load(name, shape):
    return np.fromfile(os.path.join(DIR, name), "<f2").reshape(shape)


# ---------------------------------------------------------------- face fluxes (spacing a)
def faces(c, axis):
    n = c.shape[axis]
    sh = list(c.shape); sh[axis] = n + 1
    out = np.zeros(sh)
    sl = [slice(None)] * 3; sl[axis] = slice(1, n)
    out[tuple(sl)] = 0.5 * (np.take(c, np.arange(n - 1), axis=axis) + np.take(c, np.arange(1, n), axis=axis))
    return out


def div(F, a):
    return sum(np.diff(F[ax], axis=ax) for ax in range(3)) / a


def grad_faces(chi, a):
    out = []
    for ax in range(3):
        n = chi.shape[ax]
        sh = list(chi.shape); sh[ax] = n + 1
        g = np.zeros(sh)
        sl = [slice(None)] * 3; sl[ax] = slice(1, n)
        g[tuple(sl)] = np.diff(chi, axis=ax) / a
        out.append(g)
    return out


def cg(rf, f, x0, a, tol=1e-7, itmax=3000):
    def L(x):
        g = grad_faces(x, a)
        return div([rf[k] * g[k] for k in range(3)], a)
    diag = sum(np.take(rf[k], np.arange(rf[k].shape[k] - 1), axis=k) + np.take(rf[k], np.arange(1, rf[k].shape[k]), axis=k)
               for k in range(3)) / a ** 2
    Minv = 1.0 / np.maximum(diag, 1e-30)
    b = -(f - f.mean())
    x = x0.copy()
    r = b + L(x)
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


def campi():
    os.makedirs(SCR, exist_ok=True)
    r = load("rho_lontano.f16", (NF, NA, NA, NA)).astype(np.float64)
    j = load("j_lontano.f16", (NF, 3, NA, NA, NA)).astype(np.float64)
    F = [np.zeros((NF, NA + 1, NA, NA), np.float32), np.zeros((NF, NA, NA + 1, NA), np.float32), np.zeros((NF, NA, NA, NA + 1), np.float32)]
    chi = np.zeros((NA, NA, NA))
    share, its, rot = [], [], []
    x = (np.arange(NA) - NA / 2) * AA
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    t_start = time.time()
    for k in range(NF):
        ck = os.path.join(SCR, f"campi_k{k:02d}.npz")
        if os.path.exists(ck):
            z = np.load(ck)
            for a in range(3):
                F[a][k] = z[f"F{a}"]
            share.append(float(z["share"])); its.append(int(z["it"])); rot.append(float(z["rot"])); chi = z["chi"]
            continue
        if time.time() - t_start > 450:
            print("checkpoint: rilanciare", flush=True)
            return
        k1 = (k + 1) % NF
        jm = 0.5 * (j[k] + j[k1])
        rm = 0.5 * (r[k] + r[k1])
        drdt = (r[k1] - r[k]) / DT
        jf = [faces(jm[a], a) for a in range(3)]
        rf = [np.maximum(faces(rm, a), 1e-7 * rm.max()) for a in range(3)]
        for a in range(3):
            sl = [slice(None)] * 3; sl[a] = 0; rf[a][tuple(sl)] = 0; sl[a] = -1; rf[a][tuple(sl)] = 0
        f = div(jf, AA) + drdt
        chi, it, res = cg(rf, f, chi, AA)
        g = grad_faces(chi, AA)
        Fk = [jf[a] - rf[a] * g[a] for a in range(3)]
        for a in range(3):
            F[a][k] = Fk[a]
        num = sum(((rf[a] * g[a]) ** 2).sum() for a in range(3))
        den = sum((jf[a] ** 2).sum() for a in range(3))
        share.append(float(np.sqrt(num / den)))
        its.append(it)
        Jc = np.stack([0.5 * (Fk[0][1:] + Fk[0][:-1]), 0.5 * (Fk[1][:, 1:] + Fk[1][:, :-1]), 0.5 * (Fk[2][:, :, 1:] + Fk[2][:, :, :-1])])
        rp = np.hypot(X, Y).clip(1e-9)
        rot.append(float(((-Y * Jc[0] + X * Jc[1]) / rp).sum() / np.sqrt((Jc ** 2).sum(0)).sum()))
        np.savez(ck, F0=Fk[0], F1=Fk[1], F2=Fk[2], share=share[-1], it=it, rot=rot[-1], chi=chi)
        print(k, it, res, share[-1], flush=True)
    for a, c in enumerate("xyz"):
        np.save(os.path.join(SCR, f"F{c}.npy"), F[a])
    info = {"impronta": impronta(), "correzione_su_j_mediana": float(np.median(share)), "correzione_max": float(max(share)),
            "cg_iterazioni_max": int(max(its)), "rotazione_netta_attorno_a_z_media": float(np.mean(rot))}
    json.dump(info, open(os.path.join(SCR, "campi.json"), "w"), indent=1)
    print(info)


# ---------------------------------------------------------------- balls
class Lontano:
    def __init__(self):
        self.r = load("rho_lontano.f16", (NF, NA, NA, NA)).astype(np.float64).reshape(NF, -1)
        self.F = [np.load(os.path.join(SCR, f"F{c}.npy")).astype(np.float64) for c in "xyz"]
        self.clip = 0

    def v(self, p, t):
        s = (t % T) / DT
        k = int(np.floor(s)) % NF
        u = s - np.floor(s)
        k1 = (k + 1) % NF
        x = np.clip(p / AA + NA / 2 + 0.5, 0, NA - 1e-9)
        c = np.floor(x).astype(np.int64)
        f = x - c
        lin = (c[0] * NA + c[1]) * NA + c[2]
        rho = np.maximum((1 - u) * self.r[k, lin] + u * self.r[k1, lin], 1e-12)
        v = np.empty_like(p)
        for a in range(3):
            hi = [c[0], c[1], c[2]]
            hi[a] = hi[a] + 1
            v[a] = ((1 - f[a]) * self.F[a][k][c[0], c[1], c[2]] + f[a] * self.F[a][k][hi[0], hi[1], hi[2]]) / rho
        sp = np.sqrt((v ** 2).sum(0))
        m = sp > 20.0          # alpha c units: 20 = 0.15 c (code only)
        self.clip += int(m.sum())
        v[:, m] *= 20.0 / sp[m]
        return v


def semina(rho, n, rng, a, off):
    p = rho.ravel() / rho.sum()
    cells = rng.choice(p.size, size=n, p=p)
    i, j, k = np.unravel_index(cells, rho.shape)
    return (np.stack([i, j, k]) + rng.random((3, n)) - off) * a


def palline():
    os.makedirs(SCR, exist_ok=True)
    rng = np.random.default_rng(SEED)
    r0 = load("rho_lontano.f16", (NF, NA, NA, NA))[0].astype(np.float64)
    p = semina(r0, N_A, rng, AA, NA / 2 + 0.5)
    sim = Lontano()
    out = np.zeros((NV_A + 1, 3, N_A), np.float32)
    out[0] = p
    hdt = T / VPP / 4
    lim = (NA / 2) * AA - 1e-6
    t0 = time.time()
    for iv in range(NV_A):
        for s in range(4):
            t = iv * T / VPP + s * hdt
            v1 = sim.v(p, t)
            v2 = sim.v(p + 0.5 * hdt * v1, t + 0.5 * hdt)
            p = np.clip(p + hdt * v2, -lim, lim)
        out[iv + 1] = p
    np.save(os.path.join(SCR, "lontano_pos.npy"), out)
    # the still 1s: coarse balls outside the fine block, fine balls inside
    rg = load("rho_grosso.f16", (NG, NG, NG)).astype(np.float64)
    rg_cut = rg.copy()
    c0 = NG // 2 - NB // 2
    rg_cut[c0:c0 + NB, c0:c0 + NB, c0:c0 + NB] = 0
    pg = semina(rg_cut, N_G, rng, AG, NG / 2 + 0.5)
    rf = load("rho_fine.f16", (NT, NFI, NFI, NFI))[0].astype(np.float64)
    pf = semina(rf, N_F, rng, AF, NFI / 2)
    mg, mf = rg_cut.sum() * AG ** 3, rf.sum() * AF ** 3
    np.savez(os.path.join(SCR, "fermo.npz"), pg=pg.astype(np.float32), pf=pf.astype(np.float32), mg=mg, mf=mf)
    print(json.dumps({"secondi": time.time() - t0, "clip": sim.clip, "massa_grossa_fuori_blocco": mg, "massa_fine": mf}))


# ---------------------------------------------------------------- exams
def cic(p, n, a, off):
    x = p / a + off
    i0 = np.floor(x).astype(int)
    f = x - i0
    g = np.zeros(n ** 3)
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                ii = [np.clip(i0[0] + dx, 0, n - 1), np.clip(i0[1] + dy, 0, n - 1), np.clip(i0[2] + dz, 0, n - 1)]
                w = (f[0] if dx else 1 - f[0]) * (f[1] if dy else 1 - f[1]) * (f[2] if dz else 1 - f[2])
                g += np.bincount((ii[0] * n + ii[1]) * n + ii[2], w, n ** 3)
    return g.reshape(n, n, n)


def rho_t(R, t):
    """exact two-harmonic interpolation of the far density"""
    c = np.fft.rfft(R, axis=0)
    w = 2 * np.pi * np.arange(c.shape[0]) / T
    return (c[0].real + 2 * (c[1] * np.exp(1j * w[1] * t)).real) / R.shape[0]


def esami():
    R = load("rho_lontano.f16", (NF, NA, NA, NA)).astype(np.float64)
    P = np.load(os.path.join(SCR, "lontano_pos.npy"), mmap_mode="r")
    x = (np.arange(NA) - NA / 2) * AA
    mass = R[0].sum() * AA ** 3
    cors, tvs, dzb, dzg = [], [], [], []
    for k in range(PER * NF + 1):
        iv = k * VPP // NF
        p = np.asarray(P[iv], np.float64)
        r = R[k % NF]
        h = gaussian_filter(cic(p, NA, AA, NA / 2), 1.0, mode="constant")
        g = gaussian_filter(r, 1.0, mode="constant")
        h /= h.sum(); g /= g.sum()
        cors.append(float(np.corrcoef(h.ravel(), g.ravel())[0, 1]))
        tvs.append(float(0.5 * np.abs(h - g).sum()))
        dzb.append(float(mass * p[2].mean()))
        dzg.append(float((r.sum((0, 1)) * x).sum() * AA ** 3))
    dzb, dzg = np.array(dzb), np.array(dzg)
    p0 = np.asarray(P[0], np.float64)
    ret = [float(np.sqrt(((np.asarray(P[n * VPP], np.float64) - p0) ** 2).sum(0).mean())) for n in (1, 2, 3)]
    fz = np.load(os.path.join(SCR, "fermo.npz"))
    rf = load("rho_fine.f16", (NT, NFI, NFI, NFI))[0].astype(np.float64)
    h = gaussian_filter(cic(fz["pf"].astype(np.float64), NFI, AF, NFI / 2 - 0.5), 1.0, mode="constant")
    g = gaussian_filter(rf, 1.0, mode="constant")
    c5 = float(np.corrcoef(h.ravel() / h.sum(), g.ravel() / g.sum())[0, 1])
    pf = fz["pf"].astype(np.float64)
    rr = np.sqrt((pf ** 2).sum(0))
    R4 = 4 * AF
    dens_c = float(fz["mf"] * (rr < R4).sum() / N_F / (4 / 3 * np.pi * R4 ** 3))
    xf = (np.arange(NFI) - NFI / 2 + 0.5) * AF
    XX, YY, ZZ = np.meshgrid(xf, xf, xf, indexing="ij")
    rc = np.sqrt(XX ** 2 + YY ** 2 + ZZ ** 2)
    edges = np.linspace(0, R_CUT, 10)
    nb_ = np.histogram(rr, edges)[0] * fz["mf"] / N_F
    ng_ = np.histogram(rc, edges, weights=rf)[0] * AF ** 3
    prof = [[float(a_), float(b_)] for a_, b_ in zip(nb_, ng_)]
    dmax = float(np.max(np.abs(nb_ / ng_ - 1)))
    cf = load("chir_fine.f16", (NT, NFI, NFI, NFI)).astype(np.float32)
    s = np.sign(cf)
    flips = [float((s[t] != s[t + 1]).mean()) for t in range(NT - 1)]
    ci = json.load(open(os.path.join(SCR, "campi.json")))
    res = {
        "impronta": impronta(), "manifest_commit": M2.get("commit"),
        "H1": {"corr_min": min(cors), "TV0": tvs[0], "TV_max": max(tvs), "PASS": bool(min(cors) >= 0.97 and max(tvs) <= tvs[0] + 0.05)},
        "H2": {"pearson": float(np.corrcoef(dzb, dzg)[0, 1]), "ampiezza_palline": float(np.abs(dzb).max()),
               "ampiezza_griglia": float(np.abs(dzg).max()), "ampiezza_file": float(np.abs(_DIP[:, 1]).max()),
               "PASS": bool(np.corrcoef(dzb, dzg)[0, 1] >= 0.99 and abs(np.abs(dzb).max() / np.abs(dzg).max() - 1) <= 0.03)},
        "H3": {"correzione_mediana": ci["correzione_su_j_mediana"], "correzione_max": ci["correzione_max"],
               "rotazione_netta": ci["rotazione_netta_attorno_a_z_media"], "PASS": bool(ci["correzione_su_j_mediana"] <= 0.25)},
        "H4": {"spostamento_rms_dopo_1_2_3_periodi_a0": ret},
        "H5": {"pearson_fine": c5, "densita_nucleo_palline": dens_c, "densita_nucleo_dati": float(rf[NFI // 2 - 2:NFI // 2 + 2, NFI // 2 - 2:NFI // 2 + 2, NFI // 2 - 2:NFI // 2 + 2].mean()),
               "vecchio_H5_PASS": bool(c5 >= 0.97), "profilo": prof, "scarto_max_profilo": dmax, "PASS": bool(dmax <= 0.05),
               "nota": "emendamento 1: esame = profilo radiale in 9 gusci (r < 0.21 a0), scarto massimo <= 5 %"},
        "H6": {"frazione_cambio_segno": [min(flips), max(flips)], "PASS": bool(min(flips) >= 0.5), "previsto": "FAIL"},
        "serie": {"corr": cors, "tv": tvs, "dip_palline": list(dzb), "dip_griglia": list(dzg)},
    }
    json.dump(res, open(os.path.join(OUT, "vl_idro_esami.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "serie"}, indent=1))


# ---------------------------------------------------------------- render
import vl_render as VR  # noqa: E402


def cam(iv):
    """camera: (eye, centre, fov) and half-width of the view at the centre (a0)"""
    if iv < NV_A + NV_X:
        hw = 6.0
        u = iv / (NV_A + NV_X)
        az = np.radians(-60 + 25 * u)
    elif iv < NV_A + NV_X + NV_Z:
        u = (iv - NV_A - NV_X) / (NV_Z - 1)
        e = 0.5 - 0.5 * np.cos(np.pi * u)                       # ease in-out
        hw = 6.0 * (0.22 / 6.0) ** e
        az = np.radians(-35 + 20 * e)
    else:
        hw = 0.22
        u = (iv - NV_A - NV_X - NV_Z) / NV_C
        az = np.radians(-15 + 40 * u)
    fov = 36.0
    D = hw / np.tan(np.radians(fov / 2))
    el = np.radians(14)
    eye = D * np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    return VR.camera(eye, (0, 0, 0), fov, W, H), hw, D


def draw_balls(layers, p, wt, col, rad_w, camv):
    """bin balls by their on-screen radius (world radius / distance) and accumulate"""
    sx, sy, z = VR.project(p, camv)
    f = camv[4]
    rpx = rad_w * f / np.maximum(z, 1e-6)
    b = np.clip(np.round(np.log2(np.maximum(rpx, 1.0)) * 2).astype(int), 0, len(layers) - 1)
    for ib in np.unique(b):
        m = b == ib
        VR.splat(layers[ib], sx[m], sy[m], wt[m], col[m], W, H)


C_BASE = np.array([0.80, 0.82, 1.00])
C_UP = np.array([1.00, 0.62, 0.10])
C_DN = np.array([0.18, 0.52, 1.00])


def colour(s):
    a = np.clip(np.abs(s), 0, 1) ** 0.8
    return (1 - a)[:, None] * C_BASE + a[:, None] * np.where((s > 0)[:, None], C_UP, C_DN)


class Rend:
    def __init__(self):
        self.P = np.load(os.path.join(SCR, "lontano_pos.npy"), mmap_mode="r")
        self.R = load("rho_lontano.f16", (NF, NA, NA, NA)).astype(np.float64)
        self.Rbar = self.R.mean(0)
        self.fz = np.load(os.path.join(SCR, "fermo.npz"))
        self.rf = load("rho_fine.f16", (NT, NFI, NFI, NFI)).astype(np.float32)
        self.cf = load("chir_fine.f16", (NT, NFI, NFI, NFI)).astype(np.float32)
        self.dip = np.loadtxt(os.path.join(DIR, "dipolo_lontano.txt"))
        # brightness per ball = probability carried by the ball, in units of a fine ball (same density -> same light)
        self.wA = (self.R[0].sum() * AA ** 3 / N_A) / (float(self.fz["mf"]) / N_F)

    def drel(self, p, t):
        r = rho_t(self.R, t)
        d = (r - self.Rbar) / (0.5 * self.Rbar + 0.05 * self.Rbar.max())
        return VR.trilinear(d[None], p / AA + NA / 2)[0] / 0.25

    def chirrel(self, p, tick):
        c = self.cf[tick] / np.maximum(self.rf[tick], 1e-6)
        x = np.clip(np.floor(p / AF + NFI / 2).astype(int), 0, NFI - 1)
        return c[x[0], x[1], x[2]] / CHIR_FULL

    def frame(self, iv):
        camv, hw, D = cam(iv)
        layers = np.zeros((14, 4, H, W))
        txt = {}
        if iv <= NV_A + NV_X:                       # far balls (A, then fading in the cross-fade)
            ia = min(iv, NV_A)
            t = ia * T / VPP
            al = 1.0 if iv <= NV_A else 1 - (iv - NV_A) / NV_X
            for kk, wk in enumerate((1.0, 0.5, 0.25)):
                j = ia - kk
                if j < 0:
                    break
                p = np.asarray(self.P[j], np.float64)[:, :60000]
                s = self.drel(p, j * T / VPP)
                draw_balls(layers, p, np.full(p.shape[1], al * wk * self.wA), colour(s), 0.035, camv)
            txt["t"] = t
        if iv > NV_A:                                # the still 1s, coarse + fine
            be = min(1.0, (iv - NV_A) / NV_X)
            pg = self.fz["pg"].astype(np.float64)
            pf = self.fz["pf"].astype(np.float64)
            wg = self.fz["mg"] / N_G
            wf = self.fz["mf"] / N_F
            fade_g = np.clip((np.log(hw) - np.log(0.35)) / (np.log(1.2) - np.log(0.35)), 0, 1)
            if fade_g > 0:
                draw_balls(layers, pg, np.full(pg.shape[1], be * fade_g * wg / wf), np.tile(C_BASE, (pg.shape[1], 1)), 0.6 * AG, camv)
            if iv >= NV_A + NV_X + NV_Z:
                tick = (iv - NV_A - NV_X - NV_Z) // TPF
                rr = np.sqrt((pf ** 2).sum(0))
                m = rr < R_CUT                          # sphere cut inside the fine block (a view)
                s = self.chirrel(pf[:, m], tick)
                draw_balls(layers, pf[:, m], np.ones(m.sum()), colour(s), 0.6 * AF, camv)
                txt["tick"] = tick
            else:
                draw_balls(layers, pf, np.full(pf.shape[1], be), np.tile(C_BASE, (pf.shape[1], 1)), 0.6 * AF, camv)
        img = np.zeros((4, H, W))
        for ib in range(layers.shape[0]):
            if layers[ib][3].any():
                sp = 0.7 * 2 ** (ib / 2)
                img += np.stack([VR.blur(layers[ib][c], sp) for c in range(4)])
        D_ = img[3]
        ref = np.percentile(D_[D_ > 0], 99.6) if (D_ > 0).any() else 1.0   # automatic exposure (declared)
        col = img[:3] / np.maximum(D_, 1e-12)
        lum = 1 - np.exp(-1.0 * D_ / ref)
        glow = np.stack([VR.blur(img[c], 14.0) for c in range(4)])
        gl = (glow[:3] / np.maximum(glow[3], 1e-12)) * (1 - np.exp(-0.8 * glow[3] / max(ref * 0.5, 1e-12)))
        out = 1 - (1 - np.clip(col * lum, 0, 1)) * (1 - 0.5 * np.clip(gl, 0, 1))
        bg = np.array([0.012, 0.013, 0.022])[:, None, None]
        out = np.clip(out + bg * (1 - out), 0, 1) ** (1 / 1.1)
        im = Image.fromarray((out.transpose(1, 2, 0) * 255).astype(np.uint8))
        self.testo(im, iv, hw, txt)
        return im

    def testo(self, im, iv, hw, txt):
        d = ImageDraw.Draw(im)
        F = VR.font
        d.text((28, 22), "Idrogeno vero: dal grande al piccolo", fill=(245, 245, 252), font=F(26, True))
        d.text((28, 58), "Z = 1, α vera · lontano la regola imparata dall'agglomeratore, vicino al protone la zig-zag", fill=(170, 175, 205), font=F(15))
        pm = 2 * hw * 52.9465
        d.text((28, 86), f"campo visivo {2 * hw:.3g} a0 = {pm:.3g} pm", fill=(230, 232, 245), font=F(17))
        y = H - 176
        if iv <= NV_A:
            t = txt["t"]
            zt = float(np.interp(t % T, self.dip[:, 0], self.dip[:, 1], period=T))
            ps = PER_S
            slow = 1.0 / (24.0 / VPP * ps)
            ex = int(np.floor(np.log10(slow)))
            if slow / 10 ** ex >= 9.5:
                ex += 1
            sup = str(ex).translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))
            lines = [(f"1s + 2p_z insieme: la nuvola oscilla su e giù · un periodo di Lyman-α = {ps * 1e15:.2f} fs (riga a 122,21 nm)", (255, 214, 140)),
                     (f"t = {t / T:.2f} periodi   ·   <z> = {zt:+.2f} a0 (dal calcolo)   ·   rallentato {slow / 10 ** ex:.0f}·10{sup} volte", (230, 232, 245)),
                     ("palline mosse dalla corrente della regola (corretta al minimo per la continuità) · oro = più densità della media, blu = meno", (180, 185, 210)),
                     ("regola con parità esatta (mescolamento s–p 2·10⁻⁹) · celle a0/8 (6,6 pm) · nessuno smussamento", (150, 155, 180))]
        elif iv < NV_A + NV_X + NV_Z:
            lines = [("ora lo stato fermo 1s: in un livello puro non c'è corrente, le palline stanno ferme (debito: le copie del reticolo)", (255, 214, 140)),
                     ("la corrente residua 1,7·10⁻³ dei dati è un transitorio (si parte dallo stato della griglia grossa), non fisica: non usata", (255, 160, 140)),
                     ("zoom: celle a0/16 (3,3 pm) fuori, blocco fine a0/256 (0,21 pm) attorno al protone · esposizione automatica", (230, 232, 245)),
                     ("le palline grosse svaniscono quando sono più grandi del campo visivo; dentro il blocco quelle fini", (180, 185, 210))]
        else:
            tk = txt["tick"]
            lines = [(f"sul protone · tick {tk + 1}/{NT} (1 tick = 6,9·10⁻²² s: un passo di una cella a c, poi la rotazione)", (255, 214, 140)),
                     ("colore = |R|² − |L|² sulla densità (chi corre in un verso o nell'altro), pieno a 4·10⁻³: amplificato ~250 volte", (230, 232, 245)),
                     ("nei dati è un dipolo quasi fermo lungo z (±2·10⁻³), NON un tremolio che cambia segno a ogni tick (esame H6 fallito)", (255, 160, 140)),
                     ("energia dell'1s non ancora convergente con la grandezza del blocco: −0,043 ± 0,004 % (fallimento dichiarato)", (255, 160, 140)),
                     ("si disegna solo la sfera r < 0,21 a0 dentro il blocco fine (vista tagliata) · al centro 0,314 a0⁻³ (1/π = 0,318)", (150, 155, 180))]
        for txt_, col in lines:
            d.text((28, y), txt_, fill=col, font=F(16))
            y += 26
        d.text((W - 300, H - 40), f"dati fase 2 · {impronta()}", fill=(110, 115, 140), font=F(13))


def render(a, b):
    d = os.path.join(SCR, "frames")
    os.makedirs(d, exist_ok=True)
    R = Rend()
    t0 = time.time()
    for iv in range(a, min(b, NV)):
        f = os.path.join(d, f"{iv:05d}.png")
        if os.path.exists(f):
            continue
        R.frame(iv).save(f + ".tmp.png")
        os.replace(f + ".tmp.png", f)
        if time.time() - t0 > 470:
            break


def video():
    d = os.path.join(SCR, "frames")
    for nome, vf, crf in (("1080p", "", 18), ("720p", ",scale=1280:720:flags=lanczos", 20)):
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "24", "-i", os.path.join(d, "%05d.png"), "-vf",
                        "tpad=start_duration=1:start_mode=clone:stop_duration=2.5:stop_mode=clone" + vf, "-c:v", "libx264",
                        "-preset", "slow", "-crf", str(crf), "-pix_fmt", "yuv420p", "-movflags", "+faststart",
                        os.path.join(SCR, f"video_idrogeno_bohr_{nome}.mp4")], check=True)
    for nome, iv in (("lontano", NV_A // 4), ("zoom", NV_A + NV_X + NV_Z // 2), ("protone", NV - 30)):
        Image.open(os.path.join(d, f"{iv:05d}.png")).save(os.path.join(SCR, f"video_idrogeno_bohr_{nome}.png"), optimize=True)


if __name__ == "__main__":
    m = sys.argv[1]
    if m == "campi":
        campi()
    elif m == "palline":
        palline()
    elif m == "esami":
        esami()
    elif m == "render":
        render(int(sys.argv[2]), int(sys.argv[3]))
    elif m == "still":
        Rend().frame(int(sys.argv[2])).save(sys.argv[3])
    elif m == "video":
        video()
