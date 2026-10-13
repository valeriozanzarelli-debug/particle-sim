"""locale_f3_valleys.py - phase-3 amendment (lab/out/locale_f3_amend_prereg.json): does the twisted local Wilson remove
the 4-valley splitting of the 3-D 1s (docs/BATTIMENTO_1S.md)?  Reduced model: 48^3, theta = 0.5, Z = 60, nucleus on the
corner of 8 cells or on a cell centre; levels by filter diagonalization (Rayleigh-Ritz).
Writes lab/out/locale_f3_valleys.json (one entry per configuration; exit 3 = run again)."""
import os, sys, json, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locale_core import RealWalk
from locale_cands import fermion
from zz_core import ALPHA

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "locale_f3_valleys.json")
T0 = time.time()
R = json.load(open(OUT)) if os.path.exists(OUT) else {}
NB, TH, Z = 48, 0.5, 60
ZA = Z * ALPHA
AB = 1 / (ZA * TH)                         # Bohr radius in cells (m = theta per cell)
NPER, OM = 700, np.arange(0.80, 0.9501, 0.005)
CUBE_AVG_INV_R = 2.3800                    # average of 1/r over a unit cube centred on the nucleus


def potential(place):
    x = np.arange(NB) - NB // 2
    off = 0.5 if place == "corner" else 0.0
    X, Y, Zc = np.meshgrid(x - off, x - off, x - off, indexing="ij")
    r = np.sqrt(X ** 2 + Y ** 2 + Zc ** 2)
    with np.errstate(divide="ignore"):
        V = -ZA / r
    if place == "centre":
        V[r == 0] = -ZA * CUBE_AVG_INV_R
    return V, r


VW = False
VQ = {"G": (0, 0, 0), "N_xy": (12, 12, 0), "N_xz": (12, 0, 12), "N_yz": (0, 12, 12),
      "X_x": (12, 0, 0), "X_y": (0, 12, 0), "X_z": (0, 0, 12), "R_xyz": (12, 12, 12)}


def valley_weights(v):
    """Fourier weights of a level on the 8 valleys of its copy (class (0,0,0) sites, spacing 2: q = 12 <-> k = pi/2)"""
    sub = v[:, ::2, ::2, ::2]
    F = np.sum(np.abs(np.fft.fftn(sub, axes=(1, 2, 3))) ** 2, axis=0)
    F /= F.sum()
    q = np.arange(24)
    out = {}
    for n, c in VQ.items():
        d = [np.minimum(np.abs(q - ci), 24 - np.abs(q - ci)) for ci in c]
        D2 = d[0][:, None, None] ** 2 + d[1][None, :, None] ** 2 + d[2][None, None, :] ** 2
        out[n] = float(F[D2 <= 25].sum())
    return out


def run(wname, h, place):
    w = fermion("F0", 3, TH) if wname == "F0" else fermion("F2", 3, TH, h=h)
    V, r = potential(place)
    rw = RealWalk(w, NB, 3, V=V * w.T)
    x = np.arange(NB) - NB // 2
    X, Y, Zc = np.meshgrid(x, x, x, indexing="ij")
    cls = (X % 2 == 0) & (Y % 2 == 0) & (Zc % 2 == 0)
    psi0 = np.zeros((4, NB, NB, NB), complex)
    sp = np.array([1, 0, 1, 0]) / np.sqrt(2)
    psi0[:] = sp[:, None, None, None] * (np.exp(-r / AB) * cls)[None]
    psi0 /= np.linalg.norm(psi0)
    acc = np.zeros((len(OM),) + psi0.shape, complex)
    hann = 0.5 - 0.5 * np.cos(2 * np.pi * (np.arange(NPER) + 0.5) / NPER)
    psi = psi0.copy()
    for t in range(NPER):
        ph = hann[t] * np.exp(1j * OM * t)
        acc += ph[:, None, None, None, None] * psi[None]
        psi = rw.period(psi)
    F = acc.reshape(len(OM), -1)
    G = F.conj() @ F.T
    ev, U = np.linalg.eigh(G)
    keep = ev > 1e-11 * ev.max()
    Q = (U[:, keep] / np.sqrt(ev[keep])).T @ F            # orthonormal rows
    UQ = np.array([rw.period(q.reshape(psi0.shape)).ravel() for q in Q])
    H = Q.conj() @ UQ.T
    lam, Vv = np.linalg.eig(H)
    levels = []
    p0 = psi0.ravel()
    for i in range(len(lam)):
        v = Vv[:, i] @ Q
        v /= np.linalg.norm(v)
        Uv = Vv[:, i] @ UQ
        Uv /= np.linalg.norm(Vv[:, i] @ Q)
        res = float(np.linalg.norm(Uv - lam[i] * v))
        eps = float(-np.angle(lam[i]) / (TH * w.T))         # E/m
        wgt = float(abs(np.vdot(v, p0)) ** 2)
        rho = np.sum(np.abs(v.reshape(psi0.shape)) ** 2, axis=0)
        lv = dict(E=eps, weight=wgt, resid=res, r_mean=float(np.sum(rho * r)), w_r_lt_2=float(np.sum(rho[r < 2])))
        if VW and wgt >= 0.02:
            lv["valleys"] = valley_weights(v.reshape(psi0.shape))
        levels.append(lv)
    levels.sort(key=lambda d: d["E"])
    win = [l for l in levels if 0.80 <= l["E"] <= 0.95 and l["resid"] <= 1e-6]
    big = [l for l in win if l["weight"] >= 0.02]
    tot = sum(l["weight"] for l in win)
    return dict(levels=levels, n_window=len(win), n_big=len(big), big=big,
                main_share=(max(l["weight"] for l in win) / tot) if tot > 0 else None, window_weight=tot,
                norm_end=float(np.linalg.norm(psi)), T=w.T)


CONFS = [("F0", 0, "corner"), ("F2", 0.1, "corner"), ("F2", 0.6, "corner"),
         ("F0", 0, "centre"), ("F2", 0.1, "centre"), ("F2", 0.6, "centre")]
if len(sys.argv) > 1 and sys.argv[1] == "vw":     # declared addition: valley weights of the levels (information)
    VW = True
    CONFS = [("F0", 0, "corner"), ("F2", 0.1, "corner")]
for wn, h, pl in CONFS:
    key = f"{wn}|{h}|{pl}" + ("|vw" if VW else "")
    if key in R:
        continue
    if time.time() - T0 > 300:
        print("rerun"); sys.exit(3)
    t1 = time.time()
    R[key] = run(wn, h, pl)
    print(key, f"{time.time() - t1:.0f}s", "n_big", R[key]["n_big"],
          [(round(l["E"], 5), round(l["weight"], 3), f"{l['resid']:.0e}") for l in R[key]["big"]], flush=True)
    with open(OUT, "w") as f:
        json.dump(R, f, indent=1)
print("ALL DONE")
