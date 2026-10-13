"""nucleo1s_f4_figure.py -> lab/out/nucleo1s_f4_figure.png (phase 4 of NUCLEO_1S)"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA, dirac_formula  # noqa: E402
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
R = json.load(open(os.path.join(OUT, "nucleo1s_3d_results.json")))
JN = json.load(open(os.path.join(OUT, "nucleo1s_junction1d.json")))
D1 = dirac_formula(1, -1, 60 * ALPHA)
rel = lambda om: 100 * ((1 - om) / (1 - D1) - 1)


def U_of(L):
    s = sorted(L, key=lambda x: -x["omega"]); w = o = 0.0
    for x in s:
        w += x["weight"]; o += x["weight"] * x["omega"]
        if w >= 0.4:
            break
    return o / w


fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
A = [0.5, 0.4, 0.3, 0.24, 0.2, 0.16]
ax[0].axhspan(-1, 1, color="g", alpha=0.12, label="esame: entro 1 %")
ax[0].plot(A, [rel(U_of(R[f"scan_bare_{a}"]["lines"])) for a in A], "o--", c="0.3", label="nudo: coppia cieca U")
A2 = [0.5, 0.4, 0.3, 0.24, 0.2]
ax[0].plot(A2, [100 * R[f"scan_learned2_{a}"]["rel"] for a in A2], "*-", c="tab:red", ms=9, label="fase 3: numeri di a = 0,5 ovunque")
At = [a for a in (0.4, 0.24, 0.2, 0.16) if f"scan_transfer_{a}" in R]
ax[0].plot(At, [100 * R[f"scan_transfer_{a}"]["rel"] for a in At], "D", c="tab:green", ms=9, label="fase 4: trasferito (legge lineare in a)")
AP = [a for a in (0.4, 0.24, 0.2) if f"scan_transferP_{a}" in R]
ax[0].plot(AP, [100 * R[f"scan_transferP_{a}"]["rel"] for a in AP], "x", c="tab:purple", ms=9, label="fase 4: legge a^(2-2γ) (non giudicata)")
for a, k in ((0.5, "learn2"), (0.3, "learn3")):
    ax[0].plot(a, rel(min(R[k]["evals"], key=lambda e: e["obj"])["centroid"]), "s", mfc="none", mec="k", ms=14, mew=2)
ax[0].text(0.43, -1.5, "quadrati: reticoli\ndi allenamento", ha="center", fontsize=8)
ax[0].set_xlabel("lato della cella a (λ_C)"); ax[0].set_ylabel("errore legame 1s vs Dirac (%)"); ax[0].set_ylim(-2.0, 0.6)
ax[0].set_title("(a) 1s 3-D, Z = 60: il nucleo imparato su due reticoli"); ax[0].legend(fontsize=7.5, loc="lower left")
s5, s3 = np.array(R["learn2"]["s"]), np.array(R["learn3"]["s"])
aa = np.linspace(0.12, 0.55, 50)
for i, (c, lab) in enumerate((("tab:blue", "s0 (2 siti vicini)"), ("tab:orange", "s1 (6 siti successivi)"))):
    ax[1].plot(aa, s5[i] + (s3[i] - s5[i]) * (aa - 0.5) / (0.3 - 0.5), c=c, label=lab)
    ax[1].plot([0.5, 0.3], [s5[i], s3[i]], "s", c=c, mfc="none", ms=12, mew=2)
ax[1].set_xlabel("a (λ_C)"); ax[1].set_ylabel("fattore sul potenziale"); ax[1].legend(fontsize=8)
ax[1].set_title("(b) la legge: retta fra i due reticoli imparati")
for a, mk in ((0.5, "o"), (0.25, "s")):
    for th, ls in ((0.0, "--"), (a / 2, "-")):
        rr = [r for r in JN if r["a"] == a and abs(r["thj"] - th) < 1e-9]
        ax[2].semilogy([r["k0"] for r in rr], [r["R"] for r in rr], mk + ls,
                       label=f"a = {a}, moneta del raccordo {'m a/2' if th else '0'}")
ax[2].set_xlabel("impulso del pacchetto k (mc)"); ax[2].set_ylabel("probabilità riflessa dal raccordo")
ax[2].set_title("(c) prototipo 1-D del raccordo grosso/fine (unitario a 1e-13)"); ax[2].legend(fontsize=7.5)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "nucleo1s_f4_figure.png"), dpi=110)
print("ok")
