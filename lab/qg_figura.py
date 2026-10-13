"""qg_figura.py - the study figure -> lab/out/qg_figura.png and docs/fig/quark_gluoni.png."""
import json
import os
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.ndimage import map_coordinates, spline_filter  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qg_run import SCR  # noqa: E402

LAB = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(LAB)


def piano(F, d, n=161, ext=6.0):
    q = np.array([[d, 0, 0], [0, d, 0], [0, 0, d]], float)
    c = q.mean(0)
    e1 = (q[0] - c) / np.linalg.norm(q[0] - c)
    nrm = np.ones(3) / np.sqrt(3)
    e2 = np.cross(nrm, e1)
    u = np.linspace(-ext, ext, n)
    U, V = np.meshgrid(u, u, indexing="xy")
    pts = c[:, None] + e1[:, None] * U.ravel() + e2[:, None] * V.ravel()
    val = map_coordinates(spline_filter(F, order=3, mode="grid-wrap"), pts, order=3, mode="grid-wrap", prefilter=False)
    qq = np.array([[(qk - c) @ e1, (qk - c) @ e2] for qk in q])
    return u, val.reshape(n, n), qq


def main():
    an = json.load(open(os.path.join(SCR, "qg_prod_analisi.json")))
    z = np.load(os.path.join(SCR, "qg_prod_campi.npz"))
    ver = json.load(open(os.path.join(LAB, "out", "qg_verdetto.json")))
    fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
    # (a) potential
    a = ax[0, 0]
    R = np.arange(1, len(an["V"]) + 1)
    a.errorbar(R, an["V"], an["V_err"], fmt="o", color="k", label="V(R) = ln W(R,3)/W(R,4)")
    c = an["cornell"]
    r = np.linspace(0.8, 6.5, 100)
    a.plot(r, c["V0"][0] - c["alpha"][0] / r + c["sigma"][0] * r, "C1",
           label=f"Cornell: a√σ = {c['a_sqrt_sigma'][0]:.3f}({c['a_sqrt_sigma'][1] * 1e3:.0f})")
    a.set_xlabel("R / a"); a.set_ylabel("a V"); a.set_title("(a) quark-antiquark: confinamento (esame B)")
    a.legend(fontsize=9)
    # (b), (c) maps in the triangle plane
    for k, (d, T) in enumerate(((4, 4), (3, 4))):
        a = ax[0, 1 + k]
        u, M, qq = piano(z[f"F_{d}_{T}"][2], d)
        im = a.imshow(M, origin="lower", extent=[u[0], u[-1], u[0], u[-1]], cmap="inferno")
        a.plot(qq[:, 0], qq[:, 1], "o", ms=9, mfc="none", mec="cyan", mew=2)
        a.plot([0], [0], "+", color="w", ms=10)
        plt.colorbar(im, ax=a, fraction=0.046)
        fl = an["flux"][f"{d}_{T}"]
        a.set_title(f"(b{k + 1}) f nel piano dei quark, d={d}, T={T}\nR_Y = {fl['R_Y'][0]:.2f} ± {fl['R_Y'][1]:.2f} (+ = centroide)")
        a.set_xlabel("passi di reticolo"); a.set_aspect("equal")
    # (d) centroid vs midpoints vs origin
    a = ax[1, 0]
    keys = list(an["flux"].keys())
    x = np.arange(len(keys))
    for j, (nm, col) in enumerate((("centroide", "C3"), ("medi_lati", "C0"), ("origine", "C7"))):
        v = [an["flux"][k]["f"][nm][0] for k in keys]
        e = [an["flux"][k]["f"][nm][1] for k in keys]
        a.errorbar(x + 0.2 * (j - 1), v, e, fmt="o", color=col, label=nm)
    a.set_xticks(x); a.set_xticklabels([f"d={k.split('_')[0]}, T={k.split('_')[1]}" for k in keys])
    a.set_ylabel("f = -Δ(densità d'azione)"); a.set_title("(c) Y contro Δ: centro contro metà dei lati (esame C)")
    a.legend(fontsize=9)
    # (e) part 2: D vs amplitude
    a = ax[1, 1]
    p2 = ver.get("parte2", {})
    amps = {"a002": 0.02, "a005": 0.05, "a01": 0.1, "a03": 0.3, "a1": 1.0}
    Ds = p2.get("D", {})
    xs = [amps[k] for k in Ds]; ys = [Ds[k] for k in Ds]
    a.loglog(xs, ys, "o-", color="C3", label="colori diversi (T¹, T²)")
    if p2.get("G5_colore", {}).get("D_stesso_1.0"):
        a.loglog([1.0], [p2["G5_colore"]["D_stesso_1.0"]], "s", color="C0", label="stesso colore (T¹, T¹)")
    xx = np.array([0.02, 1.0])
    if "a005" in Ds:
        a.loglog(xx, Ds["a005"] * (xx / 0.05) ** 2, ":", color="gray", label="∝ amp² (previsione)")
    a.axhline(0.10, color="C1", lw=0.8); a.axhline(0.01, color="C2", lw=0.8)
    a.set_xlabel("ampiezza del pacchetto (g A)"); a.set_ylabel("D = energia spostata")
    a.set_title("(d) gluoni: si attraversano o si urtano (esami G3-G6)"); a.legend(fontsize=9)
    # (f) conservation
    a = ax[1, 2]
    for name, col in (("a1", "C3"), ("video_grande", "C1"), ("a005", "C2")):
        try:
            dg = json.load(open(os.path.join(SCR, f"qg_ym_{name}_diag.json")))
        except FileNotFoundError:
            continue
        H = np.array(dg["H"]); g = np.array(dg["gauss_abs"])
        a.semilogy(dg["t"], np.abs(H / H[0] - 1) + 1e-16, color=col, label=f"{name}: |ΔH/H|")
        a.semilogy(dg["t"], np.abs(g / g[0] - 1) + 1e-16, "--", color=col, label=f"{name}: deriva Gauss")
    a.axhline(1e-3, color="k", lw=0.8)
    a.set_xlabel("t"); a.set_title("(e) energia e legge di Gauss conservate (G1, G2)"); a.legend(fontsize=8)
    fig.suptitle("Quark e gluoni: parte 1 = SU(3) in tempo immaginario (istantanea), parte 2 = Yang-Mills classico (film)")
    fig.tight_layout()
    out = os.path.join(LAB, "out", "qg_figura.png")
    fig.savefig(out, dpi=110)
    os.makedirs(os.path.join(ROOT, "docs", "fig"), exist_ok=True)
    shutil.copy(out, os.path.join(ROOT, "docs", "fig", "quark_gluoni.png"))
    print(out)


if __name__ == "__main__":
    main()
