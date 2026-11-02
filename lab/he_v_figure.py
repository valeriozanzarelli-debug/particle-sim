"""he_v_figure.py - figure of VIDEO ELIO -> lab/out/he_v_figure.png and docs/fig/he_v_figure.png"""
import os, sys, json, shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import he_v_render as RV

LAB = os.path.dirname(os.path.abspath(__file__))
es = json.load(open(os.path.join(LAB, "out", "he_v_esami.json")))
fig = plt.figure(figsize=(16, 9.6))
for k, (fn, tit) in enumerate([("he_v_foto_A.png", "A · elio: 24 000 palline per fotogramma"),
                               ("he_v_foto_B.png", "B · buca di Coulomb: luce = g(r2 | r_azzurro)")]):
    ax = fig.add_subplot(2, 3, k + 1)
    ax.imshow(plt.imread(os.path.join(RV.SCR, fn))); ax.set_title(tit, fontsize=10); ax.axis("off")
ax = fig.add_subplot(2, 3, 3)
d = RV.dati(); pr = d["pairs"]; dd = d["dens"]
r = np.linalg.norm(pr.reshape(-1, 3), axis=1)
h, e = np.histogram(r, np.arange(0, 4.0001, 0.05), density=True)
c = 0.5 * (e[1:] + e[:-1])
ax.plot(dd["rho1_esatta_r"], 4 * np.pi * dd["rho1_esatta_r"] ** 2 * dd["rho1_esatta"], "k-", lw=1.5, label="densità esatta (thread principale)")
ax.plot(c, h, "o", ms=3, color="tab:blue", label="palline (400 000)")
ax.set_xlim(0, 4); ax.set_xlabel("r (a0)"); ax.set_ylabel("4πr²ρ (per elettrone)")
ax.set_title(f"E1: L1 = {es['E1']['L1']:.4f}", fontsize=10); ax.legend(fontsize=8)
E4 = np.load(os.path.join(RV.SCR, "he_v_E4.npy"))
for k, (arr, tit) in enumerate([(E4[0], "g dalle coppie (ruotate)"), (E4[1], "g della mappa (usata per la luce)")]):
    ax = fig.add_subplot(2, 3, 4 + k)
    a = np.where(E4[2] > 0, arr, np.nan)
    im = ax.imshow(a.T, origin="lower", extent=[-3, 3, 0, 3], vmin=0.4, vmax=1.3, cmap="coolwarm_r")
    ax.plot([1.0], [0.0], "c*", ms=10)
    ax.set_xlabel("z lungo l'elettrone 1 (a0)"); ax.set_ylabel("distanza dall'asse (a0)")
    ax.set_title(tit + (f"  (E4: Pearson {es['E4']['pearson']:.3f})" if k == 0 else ""), fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.03)
ax = fig.add_subplot(2, 3, 6)
rows = [(k, es[k]["PASS"]) for k in ["E1", "E2", "E3", "E4", "E5"]]
ax.barh(range(len(rows)), [1] * len(rows), color=["tab:green" if p else "tab:red" for _, p in rows])
ax.set_yticks(range(len(rows))); ax.set_yticklabels([f"{k}  {'PASS' if p else 'FAIL'}" for k, p in rows])
ax.invert_yaxis(); ax.set_xticks([]); ax.set_title("esami pre-registrati", fontsize=10)
fig.suptitle("Video elio — due elettroni attorno al nucleo fermo, stato ECG del thread principale", fontsize=13)
fig.tight_layout()
out = os.path.join(LAB, "out", "he_v_figure.png")
fig.savefig(out, dpi=110)
os.makedirs(os.path.join(LAB, "..", "docs", "fig"), exist_ok=True)
shutil.copy(out, os.path.join(LAB, "..", "docs", "fig", "he_v_figure.png"))
print("ok")
