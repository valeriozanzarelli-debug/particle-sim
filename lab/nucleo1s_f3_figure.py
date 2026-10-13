"""nucleo1s_f3_figure.py -> lab/out/nucleo1s_f3_figure.png (phase 3 of NUCLEO_1S)"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA, dirac_formula  # noqa: E402
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
R = json.load(open(os.path.join(OUT, "nucleo1s_3d_results.json")))
D1 = dirac_formula(1, -1, 60 * ALPHA)
rel = lambda om: 100 * ((1 - om) / (1 - D1) - 1)


def U_of(L):
    s = sorted(L, key=lambda x: -x["omega"]); w = o = 0.0
    for x in s:
        w += x["weight"]; o += x["weight"] * x["omega"]
        if w >= 0.4:
            break
    return o / w


def spread(L):
    s = [x["omega"] for x in L if x["weight"] >= 0.05 and 0.80 <= x["omega"] <= 0.95]
    return max(s) - min(s)


fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
e1, e2 = R["learn"]["evals"], R["learn2"]["evals"]
allx = e1 + e2
c = np.array([rel(e["centroid"]) for e in allx])
sc = ax[0].scatter([e["s"][0] for e in allx], [e["s"][1] for e in allx], c=np.clip(c, -8, 8), cmap="coolwarm",
                   s=[12 + 4e4 * min(e["std"], 2e-3) for e in allx], edgecolor="k", lw=0.3)
plt.colorbar(sc, ax=ax[0], label="errore dell'energia di legame del livello medio (%)")
for key, m, lab in (("learn", "s", "criterio 1 (fase 2): solo 'una riga'"), ("learn2", "*", "criterio 2 (fase 3): una riga + livello U")):
    s = R[key]["s"]; ax[0].plot(s[0], s[1], m, ms=16, mfc="none", mec="k", mew=2, label=lab)
ax[0].set_xlabel("s0 (2 siti vicini)"); ax[0].set_ylabel("s1 (6 siti successivi)")
ax[0].set_xlim(0.44, 0.82); ax[0].set_ylim(0.65, 1.32)
ax[0].set_title("(a) prove dell'agglomeratore (pallino grande = valli aperte)"); ax[0].legend(fontsize=8, loc="upper left")
A = [0.5, 0.4, 0.3, 0.24, 0.2]
ax[1].axhspan(-1, 1, color="g", alpha=0.12, label="esame: entro 1 %")
ax[1].plot(A, [rel(R[f"scan_bare_{a}"]["quartet"]["centroid"]) for a in A], "o-", c="0.5", label="nudo: media del multipletto")
ax[1].plot(A, [rel(U_of(R[f"scan_bare_{a}"]["lines"])) for a in A], "o--", c="0.2", label="nudo: coppia cieca U")
ax[1].plot(A, [100 * R[f"scan_learned_{a}"]["rel"] for a in A], "s-", c="tab:blue", label="imparato criterio 1 (fase 2)")
ax[1].plot(A, [100 * R[f"scan_learned2_{a}"]["rel"] for a in A], "*-", c="tab:red", ms=11, label="imparato criterio 2 (fase 3)")
ax[1].set_xlabel("lato della cella a (λ_C)"); ax[1].set_ylabel("errore legame 1s vs Dirac (%)"); ax[1].set_ylim(-7, 21)
ax[1].set_title("(b) 1s 3-D, Z = 60 (a = 0,5 = reticolo d'allenamento)"); ax[1].legend(fontsize=8)
ax[2].semilogy(A, [spread(R[f"scan_bare_{a}"]["lines"]) for a in A], "o-", c="0.5", label="nudo")
ax[2].semilogy(A, [max(spread(R[f"scan_learned_{a}"]["lines"]), 1e-5) for a in A], "s-", c="tab:blue", label="criterio 1")
ax[2].semilogy(A, [spread(R[f"scan_learned2_{a}"]["lines"]) for a in A], "*-", c="tab:red", ms=11, label="criterio 2")
ax[2].axhline(0.002, c="k", ls=":", label="soglia split 0,002 mc²")
ax[2].set_xlabel("a (λ_C)"); ax[2].set_ylabel("split delle valli (mc²)"); ax[2].legend(fontsize=8)
ax[2].set_title("(c) le 4 valli: split del multipletto 1s")
plt.tight_layout(); plt.savefig(os.path.join(OUT, "nucleo1s_f3_figure.png"), dpi=110)
print("ok")
