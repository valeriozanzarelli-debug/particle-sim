"""nucleo1s_f5_figure.py -> lab/out/nucleo1s_f5_figure.png (phase 5 of NUCLEO_1S)"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "nucleo1s_f5_prereg.json")))
R = json.load(open(os.path.join(OUT, "nucleo1s_f5_results.json")))
R3 = json.load(open(os.path.join(OUT, "nucleo1s_3d_results.json")))
JN = json.load(open(os.path.join(OUT, "nucleo1s_junction3d.json")))
Zs = [20, 40, 60, 80]
ls = {20: R["learn_Z20"]["s"], 40: R["learn_Z40"]["s"], 60: R3["learn2"]["s"], 80: R["learn_Z80"]["s"]}
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
for i, (lab, yl) in enumerate((("s0 (2 siti vicini)", (0, 1.8)), ("s1 (6 siti successivi)", (0.9, 1.25)))):
    ax[i].plot(Zs, [ls[z][i] for z in Zs], "o-", c="k", ms=9, label="imparato (criterio interno, a = 0,5)")
    ax[i].plot(Zs, [P["predictions"]["PA"][f"Z{z}_a0.5"][i] for z in Zs], "s--", c="tab:blue", label="P_A media sulla cella del colore")
    if i == 0:
        for v, c in (("PB1", "tab:red"), ("PB2", "tab:orange")):
            y = [P["predictions"][v][f"Z{z}_a0.5"][0] for z in Zs]
            ax[i].plot([z for z, q in zip(Zs, y) if q is not None], [q for q in y if q is not None], "D:", c=c,
                       label=f"{v} diffusione della prima cella radiale")
        PC = json.load(open(os.path.join(OUT, "nucleo1s_f5_pc.json")))
        y = [PC[f"PC1_Z{z}_a0.01"]["s"] for z in Zs]
        ax[i].plot([z for z, q in zip(Zs, y) if q is not None], [q for q in y if q is not None], "*-", c="tab:green", ms=13,
                   label="P_C1 prima cella senza scala (emendamento 2)")
        ax[i].plot([60], [0.6052], "kx", ms=10, label="Z = 60 imparato, retta in a portata ad a → 0")
    ax[i].set_xlabel("Z"); ax[i].set_ylabel(lab); ax[i].set_ylim(*yl); ax[i].legend(fontsize=7.5)
ax[0].set_title("(a) s0: previsioni geometriche vs imparato"); ax[1].set_title("(b) s1: previsione vs imparato")
ax[0].text(20, 0.05, "Z = 20: s0 poco determinato (valli già quasi chiuse)\nZ = 80: valli non chiuse in 40 prove", fontsize=7)
for r in JN:
    c = r["case"]
    mk = "o" if c["k"][1] == 0 else "^"
    col = "tab:green" if c.get("thj", None) is None else "tab:gray"
    col = "tab:purple" if c["sig"] == 25.0 else col
    col = "tab:brown" if c["a"] == 0.5 else col
    ax[2].semilogy(c["k"][0], r["R_conservative"], mk, c=col, ms=9)
ax[2].axhline(1e-3, c="k", ls=":", lw=0.8); ax[2].axhline(3e-3, c="k", ls="--", lw=0.8)
from matplotlib.lines import Line2D
h = [Line2D([], [], marker="o", ls="", c="tab:green", label="normale, a = 0,25, σ = 10"),
     Line2D([], [], marker="o", ls="", c="tab:purple", label="normale, σ = 25 (emendamento)"),
     Line2D([], [], marker="^", ls="", c="tab:green", label="obliquo k_y = 0,785"),
     Line2D([], [], marker="o", ls="", c="tab:gray", label="senza moneta di massa"),
     Line2D([], [], marker="o", ls="", c="tab:brown", label="a = 0,5")]
ax[2].legend(handles=h, fontsize=7.5); ax[2].set_xlabel("k_x (mc)"); ax[2].set_ylabel("probabilità riflessa")
ax[2].set_title("(c) raccordo 9×9 in 3-D (soglie: 1e-3 normale, 3e-3 obliquo)")
plt.tight_layout(); plt.savefig(os.path.join(OUT, "nucleo1s_f5_figure.png"), dpi=110)
print("ok")
