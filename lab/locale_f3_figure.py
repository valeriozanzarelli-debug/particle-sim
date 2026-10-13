"""locale_f3_figure.py - phase 3 figure: vacuum pairs per valley (F0 vs F2) and valley content of the 3-D 1s levels."""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
O = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
R = json.load(open(os.path.join(O, "locale_f3_run.json")))
V = json.load(open(os.path.join(O, "locale_f3_valleys.json")))
VN = ["G", "N_xy", "N_xz", "N_yz", "X_x", "X_y", "X_z", "R_xyz"]
lab = ["Γ\n(fisico)", "N xy", "N xz", "N yz", "X x\n(π)", "X y\n(π)", "X z\n(π)", "R\n(π)"]
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
a = ax[0]
x = np.arange(8)
for j, (w, c, nm) in enumerate((("F0", "C0", "congelata"), ("F2_h0.1", "C3", "Wilson ritorto h = 0,1"))):
    a.bar(x + (j - 0.5) * 0.38, [max(R[f"S|{w}|1|{v}"]["R"], 1e-5) for v in VN], 0.38, color=c, label=nm)
a.set_yscale("log"); a.axhline(1, color="k", lw=0.8)
a.set_xticks(x); a.set_xticklabels(lab, fontsize=8)
a.set_ylabel("coppie / coppie di Schwinger (una specie)")
a.set_title("vuoto in un campo E = E_c (una copia, θ = 0,1)\nvalli a energia 0 | valli a quasi-energia π (senza massa)")
a.legend(fontsize=8)
a = ax[1]
cols = {"G": "C2", "N": "C1", "π": "C7"}
bars = []
for key, nm in (("F0|0|corner|vw", "congelata"), ("F2|0.1|corner|vw", "Wilson h=0,1")):
    for l in V[key]["levels"]:
        if "valleys" in l and 0.75 < l["E"] < 0.96:
            vv = l["valleys"]
            bars.append((f"{nm}\nE={l['E']:.4f}\npeso {l['weight']:.2f}", vv["G"], vv["N_xy"] + vv["N_xz"] + vv["N_yz"],
                         vv["X_x"] + vv["X_y"] + vv["X_z"] + vv["R_xyz"]))
xb = np.arange(len(bars))
g = np.array([b[1] for b in bars]); n = np.array([b[2] for b in bars]); p = np.array([b[3] for b in bars])
a.bar(xb, g, color=cols["G"], label="valle Γ (k = 0)")
a.bar(xb, n, bottom=g, color=cols["N"], label="3 valli N (k = π/2, π/2, 0)")
a.bar(xb, p, bottom=g + n, color=cols["π"], label="valli π")
a.set_xticks(xb); a.set_xticklabels([b[0] for b in bars], fontsize=7)
a.set_ylabel("peso di Fourier del livello"); a.set_ylim(0, 1.05)
a.set_title("1s 3-D, Z = 60, nucleo sull'angolo di 8 celle (48³, θ = 0,5)\ncongelata: ogni livello 25/25/25/25 · Wilson: tutto in Γ")
a.legend(fontsize=7, loc="lower right")
plt.tight_layout(); plt.savefig(os.path.join(O, "locale_f3_figure.png"), dpi=110)
print("ok")
