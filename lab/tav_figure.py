"""tav_figure.py - figura dello studio TAVOLA (lab/out/tav_figure.png):
(a) densita' radiale 4 pi r^2 rho(r) del fondamentale e dello stato P (gusci che si vedono), Li e Be
(b) densita' di transizione r^2 t(r): la 'forma' dell'oscillazione quando l'atomo assorbe il fotone
(c) esami: scarti contro i riferimenti (in unita' della soglia pre-registrata)
(d) una fotografia di tanti elettroni del litio (piano x-z), colorati per spin: il terzo elettrone sta fuori
Legge i dati esportati (cartelle /mnt/project-files/debiti/dati/<elemento>) e lab/out/tav_results.json."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tav_ecg as T

DEST = {"Li": os.path.join(T.OUT, "tav_dati_litio"), "Be": os.path.join(T.OUT, "tav_dati_berillio")}   # copie locali dell'export
NM = {"Li": "litio", "Be": "berillio"}
res = json.load(open(os.path.join(T.OUT, "tav_results.json")))
pre = json.load(open(os.path.join(T.OUT, "tav_prereg.json")))

fig, ax = plt.subplots(2, 2, figsize=(12, 9))
cols = {"Li": "#c0392b", "Be": "#2471a3"}
for el in ("Li", "Be"):
    f = os.path.join(DEST[el], f"{NM[el]}_densita.npz")
    if not os.path.exists(f):
        continue
    d = np.load(f); r = d["r"]
    ax[0, 0].plot(r, 4 * np.pi * r ** 2 * d["rho_fondamentale"], color=cols[el], lw=2, label=f"{el} fondamentale")
    ax[0, 0].plot(r, 4 * np.pi * r ** 2 * d["rho_eccitato_0"], color=cols[el], lw=1.2, ls="--", label=f"{el} stato P")
    ax[0, 1].plot(r, r ** 2 * d["transizione_t"] * np.sign(np.sum(d["transizione_t"])), color=cols[el], lw=2, label=el)
ax[0, 0].set_xlim(0, 12); ax[0, 0].set_xlabel("r (bohr)"); ax[0, 0].set_ylabel("4 pi r^2 rho(r)  (elettroni/bohr)")
ax[0, 0].set_title("(a) gusci: densita' radiale esatta (ECG)"); ax[0, 0].legend(fontsize=8)
ax[0, 1].set_xlim(0, 15); ax[0, 1].axhline(0, color="k", lw=0.5); ax[0, 1].set_xlabel("r (bohr)")
ax[0, 1].set_ylabel("r^2 t(r)"); ax[0, 1].set_title("(b) densita' di transizione S -> P (x cos theta)"); ax[0, 1].legend(fontsize=8)

# (c) esami
labels, vals, colors = [], [], []
thr = {"E_mEh": None}
for el in ("Li", "Be"):
    if el not in res:
        continue
    R = res[el]; D = R["diff"]
    if el == "Li":
        items = [("E", D["E_mEh"], 0.10), ("E ione", D["Eion_mEh"], 0.02), ("E 2P", D["EP_mEh"], 0.10),
                 ("ioniz.", D["IP_meV"], 3.0), ("riga", D["line_meV"], 3.0), ("f_L", D["f_L_rel"] * 100, 1.0),
                 ("f_V", D["f_V_rel"] * 100, 1.0), ("vita", D["tau_rel"] * 100, 2.0), ("<r2>", D["r2_rel"] * 100, 1.0),
                 ("rho(0)", D["rho0_rel"] * 100, 1.0)]
    else:
        items = [("E", D["E_mEh"], 1.5), ("E ione", D["Eion_mEh"], 0.10), ("ioniz.", D["IP_meV"], 30.0),
                 ("riga", D["line_meV"], 30.0), ("f_L", D["f_L_rel"] * 100, 3.0), ("f_V", D["f_V_rel"] * 100, 3.0)]
    for n, v, t in items:
        labels.append(f"{el} {n}"); vals.append(v / t); colors.append(cols[el])
y = np.arange(len(labels))
ax[1, 0].barh(y, vals, color=colors)
ax[1, 0].axvline(1, color="k", ls=":"); ax[1, 0].axvline(-1, color="k", ls=":")
ax[1, 0].set_yticks(y); ax[1, 0].set_yticklabels(labels, fontsize=7); ax[1, 0].invert_yaxis()
ax[1, 0].set_xlabel("scarto / soglia pre-registrata (|x| < 1 = PASS)")
ax[1, 0].set_title("(c) esami contro i numeri noti")

# (d) fotografia del litio
f = os.path.join(DEST["Li"], "litio_fondamentale.npy")
if os.path.exists(f):
    X = np.load(f)[:3000]
    for i, (c, lab) in enumerate((("#e67e22", "su"), ("#e67e22", None), ("#1abc9c", "giu"))):
        ax[1, 1].scatter(X[:, i, 0], X[:, i, 2], s=1.5, alpha=0.35, color=c, label=lab)
    ax[1, 1].set_xlim(-8, 8); ax[1, 1].set_ylim(-8, 8); ax[1, 1].set_aspect("equal")
    ax[1, 1].set_title("(d) litio: 3000 fotografie dei 3 elettroni (x-z)"); ax[1, 1].legend(markerscale=6, fontsize=8)
    ax[1, 1].set_xlabel("x (bohr)"); ax[1, 1].set_ylabel("z (bohr)")
fig.tight_layout()
fig.savefig(os.path.join(T.OUT, "tav_figure.png"), dpi=130)
print("scritto lab/out/tav_figure.png")
