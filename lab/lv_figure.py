"""lv_figure.py - figura dello studio LIMITI_LORENTZ -> lab/out/lv_figure.png"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lv_core as lv

HERE = os.path.dirname(os.path.abspath(__file__))
r2 = json.load(open(os.path.join(HERE, "out", "lv_r2.json")))
fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))

# (a) mappa di A(n) sulla sfera (proiezione lat-lon): ritardo della luce ~ A/6 (ka)^2
lon = np.linspace(-np.pi, np.pi, 361)
lat = np.linspace(-np.pi / 2, np.pi / 2, 181)
LO, LA = np.meshgrid(lon, lat)
n = np.stack([np.cos(LA) * np.cos(LO), np.cos(LA) * np.sin(LO), np.sin(LA)], -1)
A = lv.A_of(n)
im = ax[0].pcolormesh(np.degrees(LO), np.degrees(LA), A / 6, cmap="viridis", shading="auto")
plt.colorbar(im, ax=ax[0], label="ritardo della luce  1 - omega/(ck)  in unita' di (ka)^2")
ax[0].set_title("(a) luce: dove la regola rallenta (0 sugli assi,\n1/18 sulle diagonali; elettrone: 4 volte tanto)")
ax[0].set_xlabel("longitudine [gradi]")
ax[0].set_ylabel("latitudine [gradi]")

# (b) errore relativo vs ka: regola (palindromo) contro un tick solo
ka = np.logspace(-3, -0.3, 60)
for d, col in (([1, 1, 1], "C0"), ([1, 1, 0], "C1")):
    nn = np.array(d, float) / np.linalg.norm(d)
    K = ka[:, None] * nn[None]
    errL = 1 - lv.omega_numeric(K, 0, 1)[:, 2] / ka
    errE = 1 - lv.omega_numeric(K, 0, 0.5)[:, 3] / ka
    ax[1].loglog(ka, errL, "-", color=col, label="luce, regola, %s" % d)
    ax[1].loglog(ka, errE, "--", color=col, label="elettrone m=0, regola, %s" % d)
st = np.array([abs(1 - lv.single_tick_omega(k * np.ones(3) / np.sqrt(3))[0, 2] / k) for k in ka])
ax[1].loglog(ka, st, ":", color="C3", lw=2, label="luce, UN tick fisso x,y,z (non la regola), [1,1,1]")
ax[1].set_xlabel("k a  (= E / E_a)")
ax[1].set_ylabel("|1 - omega / (c k)|")
ax[1].set_title("(b) ordine della violazione: (ka)^2 con il palindromo,\n(ka)^1 senza (esclusa)")
ax[1].legend(fontsize=7)
ax[1].grid(alpha=0.3, which="both")

# (c) a_max per esperimento (in lunghezze di Planck)
b = r2["limiti"]["peggiore (diagonale)"]
labels, vals = [], []
short = {"tempo di volo: LHAASO": "GRB 221009A\ntempo di volo", "tempo di volo: H.E.S.S.": "Mrk 501 (H.E.S.S.)\ntempo di volo",
         "tempo di volo: MAGIC": "GRB 190114C (MAGIC)\ntempo di volo", "decadimento del fotone: LHAASO": "fotone 1,2 PeV\nnon decade",
         "decadimento del fotone: Granchio": "Granchio 1,1 PeV\nnon decade", "cavita': risonatori ottici": "cavita' ottiche",
         "cavita': risonatori a microonde": "cavita' a microonde", "INDICATIVO": "elettroni Granchio\n(indicativo)"}
for k, v in b.items():
    for s, lab in short.items():
        if k.startswith(s):
            labels.append(lab)
            vals.append(v["a_max_lP"])
y = np.arange(len(vals))
ax[2].barh(y, np.log10(vals), color=["C0"] * 3 + ["C3"] * 2 + ["C2"] * 2 + ["C7"])
ax[2].set_yticks(y)
ax[2].set_yticklabels(labels, fontsize=8)
ax[2].axvline(0, color="k", lw=1)
ax[2].text(0.3, len(vals) - 0.75, "<- lunghezza di Planck (0)", fontsize=8)
lab_a = r2["laboratorio"]["a = 0.10 lambda_C"]["a/l_P"]
ax[2].axvline(np.log10(lab_a), color="C1", lw=2)
ax[2].text(np.log10(lab_a) - 0.3, 0.2, "Laboratorio\n(a = 0,1 lambda_C)", fontsize=8, ha="right", color="C1")
ax[2].set_xlabel("log10(a massimo / l_Planck), diagonale")
ax[2].set_title("(c) quanto deve essere fine il reticolo")
ax[2].invert_yaxis()
plt.tight_layout()
p = os.path.join(HERE, "out", "lv_figure.png")
plt.savefig(p, dpi=110)
print(p)
