"""senzadoppi_figure.py - figure of the study (lab/out/senzadoppi_figura.png):
(a) the doublers: lowest |E| along Gamma -> (pi/2,0,0) [pi valley] -> (pi/2,pi/2,0) [N] -> (pi/2,pi/2,pi/2) -> Gamma
    for F0 (frozen), PA (one coin per period) and FW1 (frozen + one Wilson per period), theta = 0,1;
    solid = distance from 0 (|E|), dashed = distance from the zone edge pi/T;
(b) pairs in a field E = E_c over the whole zone without projector (r6), in units of one species;
(c) 1s of Z = 60 at theta = 0,25 (r7): spectrum of the survival amplitude, F0 vs FW1."""
import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import senzadoppi_core as s
import senzadoppi_r1 as r1
import senzadoppi_r5 as r5
from senzadoppi_r4 import lines

OUTD = os.path.join(os.path.dirname(__file__), "out")
SCR = "/tmp/claude-0/-home-claude-particle-sim/9de435bc-26c2-5d88-b70d-07efb8033e5b/scratchpad"
h = np.pi / 2
P = [np.zeros(3), np.array([h, 0, 0]), np.array([h, h, 0]), np.array([h, h, h]), np.zeros(3)]
lab = ["Γ", "π-valle\n(π/2,0,0)", "N\n(π/2,π/2,0)", "π-valle\n(π/2,π/2,π/2)", "Γ"]
pts, xs, ticks = [], [], [0.0]
x = 0.0
for a, b in zip(P[:-1], P[1:]):
    L = np.linalg.norm(b - a)
    for t in np.linspace(0, 1, 120, endpoint=False):
        pts.append(a + t * (b - a))
        xs.append(x + t * L)
    x += L
    ticks.append(x)
K = np.array(pts)
xs = np.array(xs)

fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
cols = {"F0": "k", "PA": "tab:orange", "FW1": "tab:blue"}
rules = {"F0": r1.RULES["F0"], "PA": r1.RULES["PA"], "FW1": r5.FW1(0.4)}
for nm, R in rules.items():
    E = s.quasi(R.bloch(K, 0.1), R.T)
    ax[0].plot(xs, np.min(np.abs(E), 1), color=cols[nm], lw=1.8, label=f"{nm}: |E|")
    ax[0].plot(xs, np.pi / R.T - np.max(np.abs(E), 1), color=cols[nm], lw=1.2, ls="--", label=f"{nm}: π/T − |E|")
ax[0].set_xticks(ticks)
ax[0].set_xticklabels(lab, fontsize=8)
ax[0].set_ylabel("energia (unità 1/a)")
ax[0].set_ylim(0, 1.0)
ax[0].axhline(0.1, color="gray", lw=0.6)
ax[0].set_title("(a) doppioni, θ = 0,1 (m = 0,1/a)")
ax[0].legend(fontsize=7, ncol=2)

r6 = json.load(open(os.path.join(OUTD, "senzadoppi_r6.json")))["results"]
names = ["F0", "PA", "PAW", "FW1"]
v1 = [r6[f"{n}_E1.0"]["R"] for n in names]
v2 = [r6[f"{n}_E2.0"]["R"] for n in names]
xb = np.arange(4)
ax[1].bar(xb - 0.2, v1, 0.4, label="E = E_c", color="tab:gray")
ax[1].bar(xb + 0.2, v2, 0.4, label="E = 2 E_c", color="tab:green")
ax[1].set_yscale("log")
ax[1].axhline(1, color="k", lw=0.8)
ax[1].set_xticks(xb)
ax[1].set_xticklabels(names)
ax[1].set_ylabel("coppie / coppie di UNA specie")
ax[1].set_title("(b) coppie dal vuoto, tutta la zona, senza proiettore")
for i in range(4):
    ax[1].text(i - 0.2, v1[i] * 1.1, f"{v1[i]:.2f}", ha="center", fontsize=7)
    ax[1].text(i + 0.2, v2[i] * 1.1, f"{v2[i]:.2f}", ha="center", fontsize=7)
ax[1].legend(fontsize=8)

for nm, col in (("F0", "k"), ("FW1", "tab:blue")):
    fn = os.path.join(SCR, f"senzadoppi_r7_{nm}_c.npy")
    if os.path.exists(fn):
        c = np.load(fn)
        Es, S, _ = lines(c, 0.5, 0.75, 1.0)
        ax[2].plot(Es, S, color=col, label=nm)
ax[2].axvline(0.89905, color="tab:red", lw=0.8, ls=":", label="Dirac 1s (controllo)")
ax[2].set_xlabel("energia (m c²)")
ax[2].set_ylabel("peso della riga")
ax[2].set_title("(c) 1s di Z = 60, θ = 0,25: righe dello stato")
ax[2].legend(fontsize=8)
plt.tight_layout()
fo = os.path.join(OUTD, "senzadoppi_figura.png")
plt.savefig(fo, dpi=130)
print("saved", fo)
