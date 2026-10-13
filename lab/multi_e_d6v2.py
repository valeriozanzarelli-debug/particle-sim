"""multi_e_d6v2.py - esame D6' (pre-registrazione v2): cuspide del MISTO con guida senza cuspide.
Primario: g(r) = I_mix(r) / I_VMC(r), fit pesato g0 (1 + c r + d r^2) su r < 0.25 a0 (25 gusci da 0.01), pesi = volume
del guscio, modello mediato esattamente sul guscio. Secondario (riportato): cubico pesato su h_mix.
tau = 0.0025, 0.00125, 0.000625 (semi nuovi, chiave _v2); c(tau) lineare in tau (quadratico solo se p < 0.05).
Verdetto: |c(0) - m/2| <= 0.05, solo se sigma(c(0)) <= 0.03 (altrimenti INCONCLUSIVO).
Run: python lab/multi_e_d6v2.py"""
import os
import numpy as np
from scipy.stats import chi2 as chi2d
from multi_e_core import OUT, load_json, save_json

hb = 0.01; nb = 25; NG = 20
e = np.arange(nb + 1) * hb; lo, hi = e[:-1], e[1:]
vol = 4 * np.pi * (hi ** 3 - lo ** 3) / 3
mom = lambda n: 3 * (hi ** (n + 3) - lo ** (n + 3)) / ((n + 3) * (hi ** 3 - lo ** 3))
w = np.sqrt(vol)


def blocks(name):
    z = np.load(os.path.join(OUT, f"multi_e_dmc_{name}.npz"), allow_pickle=True)
    HB = np.array(z["HB"]); return HB[:, :nb], HB[:, -1]


def fit(y, deg):
    X = np.stack([mom(k) for k in range(deg + 1)], 1)
    a, *_ = np.linalg.lstsq(X * w[:, None], y * w, rcond=None)
    return a[1] / a[0]


def c_ratio(Hm, Wm, Hv, Wv):
    return fit((Hm.sum(0) / Wm.sum()) / (Hv.sum(0) / Wv.sum()), 2)


def c_cubic(Hm, Wm):
    return fit(Hm.sum(0) / (Wm.sum() * vol), 3)


def jack(f, Hm, Wm, Hv=None, Wv=None):
    gm = np.array_split(np.arange(len(Wm)), NG)
    gv = np.array_split(np.arange(len(Wv)), NG) if Hv is not None else None
    v = []
    for i in range(NG):
        a = (np.delete(Hm, gm[i], 0), np.delete(Wm, gm[i]))
        b = (np.delete(Hv, gv[i], 0), np.delete(Wv, gv[i])) if Hv is not None else ()
        v.append(f(*a, *b))
    v = np.array(v)
    return float(np.sqrt((NG - 1) / NG * np.sum((v - v.mean()) ** 2)))


res = load_json("multi_e_dmc.json")
m = res["He_u0_dmc_tau0.01"]["m_kin"]
Hv, Wv = blocks("He_u0_vmc_tau0.01_v2")
pts = []
rec = {}
for tau in (0.0025, 0.00125, 0.000625):
    nm = f"He_u0_dmc_tau{tau:g}_v2"
    if not os.path.exists(os.path.join(OUT, f"multi_e_dmc_{nm}.npz")):
        continue
    Hm, Wm = blocks(nm)
    cr = c_ratio(Hm, Wm, Hv, Wv); sr = jack(c_ratio, Hm, Wm, Hv, Wv)
    cc = c_cubic(Hm, Wm); sc = jack(c_cubic, Hm, Wm)
    rec[f"tau{tau:g}"] = dict(tau=tau, prod_au=float(len(Wm) * 2.0), c_ratio=float(cr), sigma_ratio=sr, c_cubic=float(cc), sigma_cubic=sc)
    pts.append((tau, cr, sr, cc, sc))
    print(rec[f"tau{tau:g}"], flush=True)
out = dict(points=rec, target=m / 2, tol=0.05, sigma_required=0.03)
if len(pts) == 3:
    t, c, s, c3, s3 = map(np.array, zip(*pts))
    for lab, y, sy in (("ratio", c, s), ("cubic", c3, s3)):
        Wm_ = np.diag(1 / sy ** 2); A = np.stack([np.ones_like(t), t], 1)
        cov = np.linalg.inv(A.T @ Wm_ @ A); cf = cov @ A.T @ Wm_ @ y
        ch = float(np.sum((y - A @ cf) ** 2 / sy ** 2)); p = float(chi2d.sf(ch, 1))
        if p < 0.05:
            A = np.stack([np.ones_like(t), t, t * t], 1); cov = np.linalg.inv(A.T @ Wm_ @ A); cf = cov @ A.T @ Wm_ @ y
        out[lab] = dict(c0=float(cf[0]), sigma_c0=float(np.sqrt(cov[0, 0])), p_linear=p, fit="lineare" if p >= 0.05 else "quadratico")
    c0, s0 = out["ratio"]["c0"], out["ratio"]["sigma_c0"]
    out["verdict"] = ("INCONCLUSIVO (sigma > 0.03)" if s0 > 0.03 else ("PASS" if abs(c0 - m / 2) <= 0.05 else "FAIL"))
    print("D6':", out["verdict"], out["ratio"], "secondario cubico:", out["cubic"])
res = load_json("multi_e_dmc.json"); res["D6_v2"] = out; save_json("multi_e_dmc.json", res)
