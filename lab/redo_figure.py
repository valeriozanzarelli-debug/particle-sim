"""redo_figure.py - one figure for R1-R3 -> lab/out/redo_figure.png"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
J = lambda f: json.load(open(os.path.join(OUT, f)))

fig, ax = plt.subplots(2, 3, figsize=(15, 8.6))
# ---- R1a Schwinger rate
r1 = J("redo_r1.json")
for th, mk in ((0.1, "o"), (0.05, "s")):
    rows = [v for v in r1["scan"].values() if v["th"] == th]
    x = np.array([v["x"] for v in rows]); P0 = np.array([v["P0"] for v in rows])
    Ppi = np.array([v["Ppi"] for v in rows])
    ax[0, 0].loglog(x, P0, mk, label=f"reticolo, gap k=0, th={th}")
    if th == 0.1:
        ax[0, 0].loglog(x, P0 + Ppi, "x", color="gray", label="conto grezzo (k=0 + doppione k=pi)")
xx = np.logspace(np.log10(0.2), np.log10(80), 200)
ax[0, 0].loglog(xx, np.exp(-np.pi / xx), "k-", label="Schwinger exp(-pi m^2/eE)")
ax[0, 0].set_xlabel("E / E_crit"); ax[0, 0].set_ylabel("coppie per modo che attraversa (Gamma / (eE/2pi))")
ax[0, 0].set_title("R1a: coppie dal vuoto in campo uniforme"); ax[0, 0].legend(fontsize=7)
ax[0, 0].set_ylim(1e-6, 3)
# ---- R1c Sauter
for ell, c in ((5, "C0"), (10, "C1"), (20, "C2"), (40, "C3")):
    s = r1["sauter"][f"ell{ell}"]
    ax[0, 1].plot(s["curve_w"], s["curve_formula"], "-", color=c, lw=1)
    ax[0, 1].plot(s["curve_w"], s["curve_lat"], ".", color=c, ms=3, label=f"ell={ell} (E/Ec={s['E_over_Ecrit']:.2f})")
ax[0, 1].set_yscale("log"); ax[0, 1].set_xlabel("energia w (finestra di Klein m < w < V0-m)")
ax[0, 1].set_ylabel("coppie per modo N(w)")
ax[0, 1].set_title("R1c: gradino di Sauter (punti reticolo, linee formula esatta)"); ax[0, 1].legend(fontsize=7)
# ---- R2 flux
for n, c in ((0, "k"), (1, "C0"), (2, "C1"), (4, "C3")):
    f = os.path.join(OUT, f"redo_r2_n{n}.json")
    if not os.path.exists(f):
        continue
    d = J(f"redo_r2_n{n}.json")
    h = d["hist"]
    ax[0, 2].plot([q["k_std"] for q in h], color=c, label=f"n={n} (g n^2={d['gn2']:.1f})")
    ax[1, 0].plot([q["w_cv"] for q in h], color=c, label=f"n={n}")
ax[0, 2].axhline(0.08906 / 10, ls=":", color="gray", label="esame: spread /10")
ax[0, 2].set_yscale("log"); ax[0, 2].set_xlabel("passo del flusso"); ax[0, 2].set_ylabel("spread della curvatura std(kappa)")
ax[0, 2].set_title("R2: flusso di Ricci + pressione del flusso quantizzato"); ax[0, 2].legend(fontsize=7)
ax[1, 0].set_xlabel("passo del flusso"); ax[1, 0].set_ylabel("disuguaglianza delle lunghezze std(w)/media")
ax[1, 0].set_title("R2: le isole (lunghezze che si separano)"); ax[1, 0].legend(fontsize=7)
# ---- R3
for mode, axx in (("main", ax[1, 1]), ("fine", ax[1, 2])):
    f = os.path.join(OUT, f"redo_r3_{mode}.npz")
    if not os.path.exists(f):
        axx.set_visible(False); continue
    z = np.load(f); r = J(f"redo_r3_{mode}.json")
    t = z["t"]
    for key, c, ls in (("M16_psi2", "k", "-"), ("M16_uniform", "C0", "-"), ("M1_uniform", "C0", "--"),
                       ("M16_half", "C1", "-"), ("M1_half", "C1", "--")):
        axx.plot(t, z[key], color=c, ls=ls, label=key.replace("_", ", "))
    axx.axhline(r["M16_psi2"]["F"], color="gray", ls=":", label="pavimento di rumore F")
    axx.set_yscale("log"); axx.set_xlabel("tick"); axx.set_ylabel("H (distanza da |psi|^2)")
    axx.set_title(f"R3 surfista deterministico ({'n=24, th=0.5' if mode == 'main' else 'n=96, th=0.125'})")
    axx.legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "redo_figure.png"), dpi=110)
print("saved")
