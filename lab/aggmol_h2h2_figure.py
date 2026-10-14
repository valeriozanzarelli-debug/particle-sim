"""aggmol_h2h2_figure.py - figure of phase 3 (force between two H2): lab/out/aggmol_h2h2_figure.png
(a) E_int(R) for the four orientations and the isotropic average; (b) the pieces of the isotropic force;
(c) -E2_disp R^6 of the isotropic average against the Casimir-Polder C6."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
r = json.load(open(os.path.join(OUT, "aggmol_h2h2.json")))
G = r["geom"]
Rs = np.array(sorted({v["R"] for v in G.values()}))
x = np.linspace(Rs[0], Rs[-1], 400)
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
a = ax[0]
names = {"H": "H (paralleli)", "X": "X (incrociati)", "T": "T", "L": "L (in fila)"}
for o, col in zip("HXTL", ["tab:blue", "tab:orange", "tab:green", "tab:red"]):
    y = np.array([G[f"{o}_{R:.1f}"]["E_int"] for R in Rs]) * 1e3
    a.plot(x, CubicSpline(Rs, y)(x), color=col, lw=1, label=names[o]); a.plot(Rs, y, "o", color=col, ms=3)
iso = np.array(r["iso"]["E_int"]) * 1e3
a.plot(x, CubicSpline(Rs, iso)(x), "k", lw=2.2, label="media isotropa")
a.axhline(0, color="gray", lw=0.6)
a.axhspan(-0.132, -0.088, xmin=0, xmax=1, color="gray", alpha=0.12, label="esame: buca 0,110 mEh +-20 %")
a.axvline(6.45, color="gray", ls=":", lw=1)
a.set_ylim(-0.35, 0.4); a.set_xlabel("distanza fra i centri R (a0)"); a.set_ylabel("energia di interazione (mEh)")
v = r["verdict"]
a.set_title(f"H2-H2: buca isotropa {v['W2']['depth_mEh']:.3f} mEh a {v['W3']['R_min']:.2f} a0", fontsize=10)
a.legend(fontsize=7)
b = ax[1]
for c, lab, col in [("E_elst", "elettrostatica", "tab:purple"), ("E_exch1", "ripulsione di Pauli (4 elettroni)", "tab:red"),
                    ("E_ind", "induzione", "tab:olive"), ("E_disp", "dispersione (van der Waals)", "tab:blue"),
                    ("E_int", "totale", "k")]:
    y = np.array(r["iso"][c]) * 1e3
    b.plot(Rs, y, "o-", ms=3, color=col, label=lab)
b.axhline(0, color="gray", lw=0.6); b.set_ylim(-0.6, 0.6)
b.set_xlabel("R (a0)"); b.set_ylabel("mEh"); b.set_title("i pezzi della forza (media isotropa)", fontsize=10)
b.legend(fontsize=7)
c = ax[2]
y = -np.array(r["iso"]["E_disp"]) * Rs ** 6
c.plot(Rs, y, "o-", color="tab:blue", label="-E_disp R^6 (curva)")
c.axhline(r["mono"]["C6_CP"], color="tab:blue", ls="--", label=f"C6 Casimir-Polder {r['mono']['C6_CP']:.2f}")
c.axhspan(10.8, 13.2, color="gray", alpha=0.12, label="esame: 12 +-10 %")
c.set_xlabel("R (a0)"); c.set_ylabel("a.u."); c.set_title("coda: la dispersione va come C6/R^6", fontsize=10)
c.legend(fontsize=7)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "aggmol_h2h2_figure.png"), dpi=110)
print("ok")
