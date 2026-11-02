"""vl_esami.py - VIDEO ATOMO: the pre-registered exams (lab/out/vl_prereg.json) on the moved balls.

    python3 lab/vl_esami.py          -> lab/out/vl_esami.json (+ series for the figure in lab/out/vl_serie.json)
"""
import json
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vl_comune import SCR, OUT, A, NA, NL, H  # noqa: E402
from atomo_export_leggi import carica  # noqa: E402


def cic(p, n, idx):
    """cloud-in-cell deposit of balls at fractional cell-centre indices idx(p) on an n^3 grid"""
    x = idx(p)
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


def e1_frame(p, rho_s):
    h = gaussian_filter(cic(p, NA, lambda q: q / A + NA / 2 - 0.5), 1.0, mode="constant")
    r = gaussian_filter(rho_s.astype(np.float64), 1.0, mode="constant")
    h /= h.sum(); r /= r.sum()
    return float(np.corrcoef(h.ravel(), r.ravel())[0, 1]), float(0.5 * np.abs(h - r).sum())


def pball(p, d, k, rho_s_k):
    w1 = d.man["omega_1s"]
    t = d.fotogrammi["t"][k]
    c = np.clip(np.floor(p / A + NA / 2).astype(int), 0, NA - 1)
    r0, f0 = d.rho(0), d.fase(0)
    rk, fk = d.rho(k), d.fase(k)
    sel = (c[0], c[1], c[2])
    val = np.sqrt(r0[sel] * rk[sel]) * np.exp(1j * (fk[sel] - f0[sel] - (1 - w1) * t)) / np.maximum(rho_s_k[sel], 1e-14)
    return float(1 - abs(val.mean()) ** 2)


if __name__ == "__main__":
    d = carica()
    fr = d.fotogrammi
    dati = np.load(os.path.join(SCR, "el_dati.npy") if len(sys.argv) < 2 else sys.argv[1], mmap_mode="r")
    rho_s = np.load(os.path.join(SCR, "rho_s.npy"), mmap_mode="r")
    pc = np.load(os.path.join(SCR, "pexc_campo.npy"))
    ks = [k for k in range(201) if np.isfinite(dati[k, 0, 0])]
    cor, tv, dip, pb = {}, {}, {}, {}
    for k in ks:
        p = np.asarray(dati[k], np.float64)
        cor[k], tv[k] = e1_frame(p, rho_s[k])
        dip[k] = float(p[2].mean())
        pb[k] = pball(p, d, k, np.asarray(rho_s[k]))
    K = np.array(ks)
    C = np.array([cor[k] for k in ks]); T = np.array([tv[k] for k in ks])
    from vl_comune import impronta
    res = {"fotogrammi_valutati": len(ks), "impronta_esportazione": impronta()[0], "manifest_mtime_utc": impronta()[1]}
    res["E1"] = {"corr_min": float(C.min()), "corr_min_fotogramma": int(K[C.argmin()]), "TV0": tv[0] if 0 in tv else None,
                 "TV_max": float(T.max()), "TV_max_fotogramma": int(K[T.argmax()]),
                 "PASS": bool(C.min() >= 0.97 and T.max() <= T[0] + 0.05)}
    m = K >= 5
    Dz = np.array([dip[k] for k in K[m]]); Mz = fr["dipolo_z"][K[m]]
    res["E2"] = {"pearson": float(np.corrcoef(Dz, Mz)[0, 1]), "rms_lambdaC": float(np.sqrt(((Dz - Mz) ** 2).mean())),
                 "ampiezza_max_palline": float(np.abs(Dz).max()), "ampiezza_max_manifest": float(np.abs(Mz).max())}
    res["E2"]["PASS"] = bool(res["E2"]["pearson"] >= 0.98 and res["E2"]["rms_lambdaC"] <= 0.10)
    Pb = np.array([pb[k] for k in K[m]])
    res["E3a"] = {"max_diff_palline_vs_reticolo": float(np.abs(Pb - pc[K[m]]).max())}
    res["E3a"]["PASS"] = bool(res["E3a"]["max_diff_palline_vs_reticolo"] <= 0.02)
    m2 = K >= 25
    Pb2 = np.array([pb[k] for k in K[m2]])
    res["E3b"] = {"max_diff_palline_vs_manifest": float(np.abs(Pb2 - fr["p_exc"][K[m2]]).max()) if m2.any() else None,
                  "media_diff": float((Pb2 - fr["p_exc"][K[m2]]).mean()) if m2.any() else None}
    res["E3b"]["PASS"] = bool(m2.any() and res["E3b"]["max_diff_palline_vs_manifest"] <= 0.02)
    # E6 light balls vs u
    e6 = {}
    E = np.load(os.path.join(SCR, "E.npy"), mmap_mode="r"); B = np.load(os.path.join(SCR, "B.npy"), mmap_mode="r")
    for k in (10, 20, 30):
        f = os.path.join(SCR, f"lu_dati_k{k:03d}.npy")
        if not os.path.exists(f):
            continue
        p = np.load(f).astype(np.float64)
        h = cic(p, NL, lambda q: (q / 4.0 + NL - 0.5) / 2.0)
        u = 0.5 * ((np.asarray(E[k], np.float64) ** 2).sum(0) + (np.asarray(B[k], np.float64) ** 2).sum(0))
        X = (2 * np.arange(NL) + 0.5 - NL) * 4.0
        ins = np.abs(X) + H / 2 < 208
        mm = ins[:, None, None] & ins[None, :, None] & ins[None, None, :]
        e6[k] = float(np.corrcoef(h[mm], u[mm])[0, 1])
    res["E6"] = {"pearson": e6, "PASS": bool(len(e6) == 3 and min(e6.values()) >= 0.90)}
    camp = json.load(open(os.path.join(OUT, "vl_campi.json")))
    L = camp["luce"]
    res["E4"] = {"media": L["B_vs_rxE_R>60_corr_media_60_200"], "min": L["B_vs_rxE_min"],
                 "PASS": bool(L["B_vs_rxE_R>60_corr_media_60_200"] >= 0.9 and L["B_vs_rxE_min"] >= 0.8)}
    res["E5"] = {"rapporto": L["rapporto_uscita_su_persa"], "PASS": bool(0.67 <= L["rapporto_uscita_su_persa"] <= 1.5)}
    # E8 (amendment 3): turning speed of the balls around z vs Z alpha c sin(theta), reported only
    za = d.man["Z"] * d.man["alpha"]
    e8 = {}
    for k in (1, 25, 50, 100, 150, 199):
        if not (np.isfinite(dati[k, 0, 0]) and np.isfinite(dati[k + 1, 0, 0])):
            continue
        a_ = np.asarray(dati[k], np.float64); b_ = np.asarray(dati[k + 1], np.float64)
        dphi = np.angle(np.exp(1j * (np.arctan2(b_[1], b_[0]) - np.arctan2(a_[1], a_[0]))))
        rp = np.hypot(a_[0], a_[1]); r = np.sqrt((a_ ** 2).sum(0))
        e8[k] = float((dphi / 8.0 * rp).mean() / (za * rp / np.maximum(r, 1e-9)).mean())
    res["E8_giro_su_Dirac"] = e8
    for f in ("elettrone_DONE.json", "luce_DONE.json"):
        if os.path.exists(os.path.join(SCR, f)):
            res[f[:-10] + "_corsa"] = json.load(open(os.path.join(SCR, f)))
    print(json.dumps(res, indent=1))
    if len(sys.argv) < 2:
        with open(os.path.join(OUT, "vl_esami.json"), "w") as f:
            json.dump(res, f, indent=1)
        with open(os.path.join(OUT, "vl_serie.json"), "w") as f:
            json.dump({"k": ks, "corr": [cor[k] for k in ks], "tv": [tv[k] for k in ks], "dip_palline": [dip[k] for k in ks],
                       "p_palline": [pb[k] for k in ks], "p_reticolo": [float(pc[k]) for k in ks],
                       "p_manifest": [float(fr["p_exc"][k]) for k in ks], "dip_manifest": [float(fr["dipolo_z"][k]) for k in ks]}, f)
