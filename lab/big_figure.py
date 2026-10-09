"""big_figure.py - lab/out/big_figure.png: R1 reader, R2 control, highways, reconstructed shape, H(t), spectrum."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED = "#222222", "#888888"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
n = 24
fig, ax = plt.subplots(2, 3, figsize=(15, 9.5))

# (a) R1
r1 = json.load(open(os.path.join(OUT, "big_r1.json")))
cls = ["shuffle", "noise", "random_graph", "random_state"]
a = ax[0, 0]
for k, c in enumerate(cls):
    v = [q["rel"][0] if q["rel"] else 0 for q in r1["nulls"][c]]
    a.scatter(np.full(len(v), k) + np.random.default_rng(k).uniform(-0.15, 0.15, len(v)), v, s=12, color=MUTED)
for k, (s, lab) in enumerate((("ring", "anello"), ("torus_wilson", "toro Wilson"), ("line_x_hidden_circle", "linea×cerchio"))):
    v = r1["known_shapes"][s]["rel"][0]
    a.scatter([4 + k], [v], s=40, color=C1, zorder=3)
a.scatter([6], [r1["cylinder_hidden_bar"]["rel"]], s=40, facecolor="white", edgecolor=C1, zorder=3)
a.axhline(0.4, color=C2, lw=1.5, ls="--"); a.text(-0.4, 0.42, "soglia 0,4", color=INK)
a.set_xticks(range(7)); a.set_xticklabels(["mescol.", "rumore", "grafo", "stato", "anello", "toro W.", "lin×cerc."],
                                          rotation=30)
a.set_ylabel("buco più grande: persistenza / diametro")
a.set_title("R1: lettore riparato\n100 nulli nuovi: 0 falsi positivi; forme note 3/3", fontsize=10)
a.text(6.1, r1["cylinder_hidden_bar"]["rel"] + 0.03, "cerchio\nnascosto", fontsize=7, color=MUTED)

# (b) R2
a = ax[0, 1]
z = np.load(os.path.join(OUT, "big_r2_theta0.125.npz"))
for M, c in (("M1", C2), ("M16", C1)):
    a.semilogy(z["t"], z[M], color=c, lw=2, label={"M1": "M = 1 (controllo, niente vortici)", "M16": "M = 16"}[M])
a.axhline(0.000875, color=MUTED, lw=1, ls=":"); a.text(900, 0.00062, "pavimento F", color=MUTED)
a.set_xlabel("tick"); a.set_ylabel("H (distanza da |ψ|²)"); a.legend(frameon=False)
a.set_title("R2 (θ = 0,125, λ = 0): FAIL\nil controllo M = 1 si mescola come M = 16", fontsize=10)

# (c) H(t) in B
a = ax[0, 2]
Hs = {k: [] for k in ("reinf_uniform", "plain_uniform", "reinf_psi2")}
for s in (1, 2, 3):
    zb = np.load(os.path.join(OUT, f"big_b_seed{s}.npz"))
    t = zb["t"]
    for k in Hs:
        Hs[k].append(zb[k + "_H"])
for k, c, lab in (("reinf_uniform", C2, "con rinforzo, partenza uniforme"), ("plain_uniform", C1, "senza rinforzo"),
                  ("reinf_psi2", C3, "con rinforzo, partenza |ψ|²")):
    a.semilogy(t, np.mean(Hs[k], 0), color=c, lw=2, label=lab)
a.axhline(3 * 0.000875, color=MUTED, lw=1, ls="--"); a.text(20, 0.0029, "3F", color=MUTED)
a.set_xlim(0, 1500); a.set_xlabel("tick (media su 3 semi; dopo 1500 tutto resta al pavimento)"); a.set_ylabel("H")
a.legend(frameon=False); a.set_title("B: H(t), media 3 semi\nnessun blocco di Born (non valutato)", fontsize=10)

# (d) highways
a = ax[1, 0]
zb = np.load(os.path.join(OUT, "big_b_seed1.npz"))
W = zb["reinf_uniform_W"]
segs, ws = [], []
for x in range(n):
    for y in range(n):
        if x < n - 1:
            segs.append([(x, y), (x + 1, y)]); ws.append(W[0][x, y])
        if y < n - 1:
            segs.append([(x, y), (x, y + 1)]); ws.append(W[1][x, y])
ws = np.array(ws)
lc = LineCollection(segs, linewidths=0.3 + 3.5 * ws, cmap="Blues", norm=plt.Normalize(0, 1))
lc.set_array(ws); a.add_collection(lc)
a.set_xlim(-1, n); a.set_ylim(-1, n); a.set_aspect("equal"); a.axis("off")
fig.colorbar(lc, ax=a, fraction=0.04, label="peso w del legame")
j = json.load(open(os.path.join(OUT, "big_b_seed1.json")))["reinf_uniform"]
a.set_title(f"B: pesi finali (seme 1, t = 5000)\nGini {j['gini_end']:.2f}, top 10% = {j['top10_end']:.0%}: isole, non strade", fontsize=10)

# (e) reconstructed shape
a = ax[1, 1]
X = zb["shape_X"]; st = zb["shape_sites"]
xs, ys = st // n, st % n
sc = a.scatter(X[:, 0], X[:, 1], c=xs, cmap="Blues", s=8, vmin=-6, vmax=n)
a.set_aspect("equal"); a.set_xlabel("MDS 1"); a.set_ylabel("MDS 2")
a.set_title("B: forma letta dalla MI (seme 1)\n1 pezzo, 0 buchi: un quadrato pieno", fontsize=10)
fig.colorbar(sc, ax=a, fraction=0.04, label="x vero del sito (solo per colorare)")

# (f) spectrum
a = ax[1, 2]
ref = json.load(open(os.path.join(OUT, "big_b_ref.json")))
a.semilogy(np.arange(1, 13), ref["omega_zero"], "o", color=C1, ms=6, label="senza rinforzo (w = 0,3 ovunque)")
for s, mk in ((1, "s"), (2, "^"), (3, "D")):
    g = json.load(open(os.path.join(OUT, f"big_b_seed{s}.json")))["shape_reinf_uniform"]
    a.semilogy(np.arange(1, 13) + 0.1 * s, g["omega_zero"], mk, color=C2, ms=5, mfc="none",
               label="con rinforzo" if s == 1 else None)
a.set_xlabel("livello k"); a.set_ylabel("|ω| più bassi (masse del cammino)")
a.legend(frameon=False, loc="lower right"); a.set_title("B: spettro più basso\nrapporti 1 : 1,07 : 3,1 : 5,3 (niente salti grandi)", fontsize=10)

fig.tight_layout()
fig.savefig(os.path.join(OUT, "big_figure.png"), dpi=130)
print("saved")
