"""locale_figure.py - figure of the study REGOLA LOCALE -> lab/out/locale_figure.png (reads locale_r2/r3 json)."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from locale_core import PI
from locale_cands import fermion

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
r2 = json.load(open(os.path.join(OUT, "locale_r2.json")))
r3 = json.load(open(os.path.join(OUT, "locale_r3.json")))
fig, ax = plt.subplots(2, 2, figsize=(12, 9))

# (a) fermion spectrum along a path of the reduced zone
pts = {"Γ": (0, 0, 0), "N": (PI / 2, PI / 2, 0), "H": (PI / 2, 0, 0), "P": (PI / 2, PI / 2, PI / 2)}
path = ["Γ", "N", "H", "Γ", "P"]
K, ticks, s = [], [0], 0
for a, b in zip(path[:-1], path[1:]):
    A, B = np.array(pts[a]), np.array(pts[b])
    for t in np.linspace(0, 1, 120, endpoint=False):
        K.append(A + t * (B - A))
    s += 120; ticks.append(s)
K.append(np.array(pts[path[-1]]))
K = np.array(K)
a = ax[0, 0]
for name, col, lab in (("F0", "C3", "regola congelata"), ("F2", "C0", "con Wilson locale (Y, h=0.1)")):
    w = fermion(name, 3, 0.1)
    eps, _ = w.phases(K)
    E = np.sort(np.abs(eps), axis=1) / w.T
    a.plot(E[:, 0], col, lw=2, label=lab)
    a.plot(E[:, 2], col, lw=0.8, alpha=0.6)
a.axhline(0.1, color="k", ls=":", lw=0.8)
a.text(5, 0.11, "m = 0,1", fontsize=8)
a.set_xticks(ticks); a.set_xticklabels(path)
a.set_ylabel("energia più bassa |E| (unità 1/a)")
a.set_title("(a) doppioni: valli finte a E = m nei punti N (rosso) — tolte dal Wilson locale (blu)")
a.legend(fontsize=8); a.set_ylim(0, 1.2)

# (b) light: worst error vs cells per lambda
a = ax[0, 1]
Ns = [8, 12, 16, 24, 32, 48, 64]
for name, col, lab in (("L0_3D", "C3", "congelata 3-D"), ("L2_3D", "C0", "S4 (passo indietro) 3-D"),
                       ("L0_2D", "C1", "congelata 2-D"), ("L1_2D", "C2", "esagonale 2-D"), ("L4b_3D", "C7", "bcc 3-D")):
    g = r3["M"][name]["grid"]
    a.loglog(Ns, [g[str(n)]["diffraction"] for n in Ns], "o-", color=col, label=lab)
a.axhline(0.01, color="k", ls=":", lw=0.8)
a.set_xlabel("celle per lunghezza d'onda"); a.set_ylabel("errore peggiore sulla diffrazione")
a.set_title("(b) luce: quanto sbaglia l'allargamento del fascio"); a.legend(fontsize=8)

# (c) kinetic mass error vs theta
a = ax[1, 0]
for name, col, lab in (("F0", "C3", "congelata (ordine 2)"), ("F3", "C0", "S4 (ordine 4)"), ("F2", "C2", "Wilson locale")):
    e = r2["D8"][name]["err"]
    th = sorted(float(t) for t in e)
    a.loglog(th, [e[str(t)] if str(t) in e else e[t] for t in th], "o-", color=col, label=lab)
a.set_xlabel("θ = m c a / ħ"); a.set_ylabel("|m_cinetica/m − 1|")
a.set_title("(c) elettrone: errore della massa cinetica"); a.legend(fontsize=8)

# (d) cells per lambda at 1% and cost
a = ax[1, 1]
names = ["L0_3D", "L3_3D", "L4a_3D", "L4b_3D", "L5_3D", "L2_3D", "L0_2D", "L1_2D", "L2_2D"]
Nw = [r3["M"][n]["N1_worst"] for n in names]
cost = [r3["M"][n]["cost"] / r3["M"]["L0_3D" if n.endswith("3D") else "L0_2D"]["cost"] for n in names]
x = np.arange(len(names))
a.bar(x, Nw, color=["C3" if n.startswith("L0") else "C0" if n.startswith("L2") else "C7" for n in names])
for i, c in enumerate(cost):
    a.text(i, Nw[i] + 1, f"costo\n×{c:.2f}", ha="center", fontsize=7)
a.axhline(20, color="k", ls=":", lw=0.8)
a.set_xticks(x); a.set_xticklabels(names, rotation=40, fontsize=8)
a.set_ylabel("celle/λ per l'1% (fase, gruppo, diffrazione)")
a.set_title("(d) luce: celle per λ necessarie e costo relativo")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "locale_figure.png"), dpi=95)
print("saved")
