"""aggfull_p7_figure.py - phase 7 figure -> lab/out/aggfull_p7_figure.png"""
import os, sys, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_core import dirac_binding, ALPHA
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
G = json.load(open(os.path.join(OUT, "aggfull_p7_glue_err.json")))
C = json.load(open(os.path.join(OUT, "aggfull_p7_chain.json")))
TH = (1 / ALPHA) / 4096
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
J = [24, 32, 48, 64]
for t, lab, mk in (("th1", "cella catena θ", "o"), ("th2", "θ/2", "s"), ("th4", "θ/4", "^")):
    for n, ls in ((1, "-"), (2, "--")):
        ax[0].plot(J, [100 * G[f"{t}_J{j}_{n}s"] for j in J], mk + ls, label=f"{n}s, {lab}")
ax[0].axhspan(-1, 1, color="g", alpha=0.1); ax[0].axhline(0, color="gray", lw=0.5)
ax[0].set_xlabel("J0 (celle grosse nella catena)"); ax[0].set_ylabel("errore α⁴ (differenza Z=1/Z=2), %")
ax[0].set_title("colla: il plateau sparisce affinando la catena"); ax[0].legend(fontsize=7)
for v, tf, mk in (("c1", 1.0, "o"), ("c2", 0.5, "s")):
    m = np.tan(TH * tf) / (TH * tf)
    for n, ls in ((1, "-"), (2, "--")):
        Z = np.arange(1, 5)
        y = [C[f"{v}_Z{z}_n{n}"]["B"] / m / dirac_binding(n, -1, z * ALPHA) - 1 for z in Z]
        ax[1].plot(Z, 1e6 * np.array(y), mk + ls, label=f"{n}s, cella θ{'/2' if tf < 1 else ''}")
ax[1].axhline(0, color="gray", lw=0.5)
ax[1].set_xlabel("Z"); ax[1].set_ylabel("errore relativo dopo la massa tanθ/θ (10⁻⁶)")
ax[1].set_title("catena da sola: errore ∝ Z (un termine Z³, non α⁴)"); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "aggfull_p7_figure.png"), dpi=110)
