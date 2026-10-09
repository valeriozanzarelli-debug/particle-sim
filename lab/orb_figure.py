"""orb_figure.py - lab/out/orb_figure.png from orb_o1/o2/o3 json."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from orb_core import OUT

o1 = json.load(open(f"{OUT}/orb_o1.json")); o2 = json.load(open(f"{OUT}/orb_o2.json")); o3 = json.load(open(f"{OUT}/orb_o3.json"))
fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
a1, a2 = np.array([1, 0]), np.array([0.5, np.sqrt(3) / 2])

# (A) T2/Z3 fundamental domain
A = ax[0, 0]
L = 9
pts = np.array([i / L * a1 + j / L * a2 for i in range(L + 1) for j in range(L + 1)])
A.scatter(pts[:, 0], pts[:, 1], s=6, c="#bbb")
A.add_patch(Polygon([0 * a1, a1, a1 + a2, a2], fill=False, ec="k", lw=1))
fp = [np.zeros(2), (a1 + a2) / 3, 2 * (a1 + a2) / 3]
A.add_patch(Polygon([fp[0], a1, fp[2], fp[1], fp[0]], fc="#4C72B0", alpha=0.25, ec="#4C72B0"))
for p in fp:
    A.plot(*p, "o", c="#C44E52", ms=10)
A.set_title("T²/Z3: toro esagonale, 3 punti fissi (rosso)\nzona blu = dominio fondamentale (= sfera con 3 punte)", fontsize=10)
A.set_aspect("equal"); A.axis("off")

# (B) T4/Z2 slice: T2/Z2 pillowcase, 4 fixed points per 2-D slice (16 = 4x4 in 4-D)
B = ax[0, 1]
B.add_patch(Polygon([[0, 0], [1, 0], [1, 1], [0, 1]], fill=False, ec="k"))
B.add_patch(Polygon([[0, 0], [0.5, 0], [0.5, 1], [0, 1]], fc="#55A868", alpha=0.25, ec="#55A868"))
for p in [(0, 0), (0.5, 0), (0, 0.5), (0.5, 0.5), (1, 0), (0, 1), (1, 1), (1, 0.5), (0.5, 1)]:
    B.plot(*p, "o", c="#C44E52", ms=8, alpha=1 if max(p) <= 0.5 else 0.3)
B.text(0.52, 0.52, "x → −x", fontsize=9)
B.set_title("T⁴/Z2 (fetta 2-D): 4 punti fissi per fetta,\n4×4 = 16 in 4-D (rete n⁴, n = 4 e 6)", fontsize=10)
B.set_aspect("equal"); B.axis("off"); B.set_xlim(-0.1, 1.1); B.set_ylim(-0.1, 1.1)

# (C) holonomy: invariant spinors
C = ax[0, 2]
names = [("K3 = T⁴/Z2", "K3: T4/Z2 x->-x (angles pi,pi)", True), ("T⁶/Z3", "T6/Z3 (w,w,w)", False),
         ("T²/Z3", "T2/Z3", False), ("controllo Z2\n1 coord.", "control: Z2 on one complex coord (pi,0)", True),
         ("controllo Z3\n1 coord.", "control: T6, Z3 on one T2 (w,0,0)", False)]
vals, loop = [], []
for lab, k, lst in names:
    d = o1[k][-1] if lst else o1[k]
    vals.append(max(l["invariant_spinors"] for l in d["lifts"].values()))
    loop.append(d["with_extra_U1_phase_max_invariant"])
x = np.arange(len(names))
C.bar(x - 0.2, vals, 0.4, color=["#4C72B0", "#4C72B0", "#999", "#999", "#999"], label="spinori invarianti (rotazione vera)")
C.bar(x + 0.2, loop, 0.4, color="none", ec="#DD8452", hatch="//", label="con una fase U(1) in più (trucco)")
C.set_xticks(x); C.set_xticklabels([n[0] for n in names], fontsize=8)
C.set_title("O1: spinore che torna uguale dopo ogni giro\n(riflessione x1→−x1 con massa: nessun operatore)", fontsize=10)
C.legend(fontsize=8); C.set_ylabel("numero")

# (D) Betti / Euler
D = ax[1, 0]
rows = [("T²/Z3", o2["T2/Z3 L=9"]["betti_quotient"], 2, 2), ("T⁴/Z2", o2["T4/Z2 n=4"]["betti_quotient"], 8, 24),
        ("T⁶/Z3", o2["T6/Z3"]["betti_quotient"], 18, 72)]
D.axis("off")
txt = "O2: contati con gli operatori delle palline\n\n"
txt += f"{'':8s}{'punti fissi':>11s}  {'Betti del quoziente':<22s}{'χ':>4s}  {'χ dopo lisciatura':>17s}\n"
fpts = [o2["T2/Z3 L=9"]["fixed_points"], o2["T4/Z2 n=4"]["fixed_points"], o2["T6/Z3"]["fixed_points_brute"]]
for (n, b, chi, chir), f in zip(rows, fpts):
    txt += f"{n:8s}{f:>11d}  {str(tuple(b)):<22s}{chi:>4d}  {chir:>17d}\n"
txt += ("\nK3: 8 + 16 (una sferetta per punto) = 24, b2 = 6+16 = 22\n"
        f"T⁶/Z3: h11 = 9 + 27 = 36, h21 = 0 → 36 famiglie\n"
        f"zeri di Dirac sull'orbifold: K3 {o2['Dirac T4/Z2 n=6']['zero_modes']} (stessa chiralità),"
        f" T²/Z3 {o2['Dirac T2/Z3 L=9']['zero_modes']}\n"
        "la lisciatura (+1, +2 per punto) è messa a mano")
D.text(0, 0.95, txt, va="top", family="monospace", fontsize=9)

# (E) magnetised counts
E = ax[1, 1]
Ms = list(range(1, 13))
for wl, mk, col in [("(0.000,0.000)", "o", "#4C72B0"), ("(0.333,0.667)", "s", "#DD8452"), ("(0.500,0.500)", "^", "#55A868")]:
    for q in range(3):
        xs, ys = [], []
        for M in Ms:
            v = o3["b_magnetised"].get(f"M={M} W={wl}")
            if isinstance(v, dict):
                xs.append(M + (q - 1) * 0.15); ys.append(v["counts"][q])
        E.scatter(xs, ys, marker=mk, color=col, s=25, label=f"W = {wl}" if q == 0 else None)
E.axhline(3, color="#C44E52", lw=0.8, ls="--")
E.set_xlabel("flusso M (quanti) attraverso il toro nascosto"); E.set_ylabel("modi leggeri per settore Z3")
E.set_title("O3: T²/Z3 con flusso + linea di Wilson:\nmodi leggeri nei 3 settori (linea rossa = 3)", fontsize=10)
E.legend(fontsize=8)

# (F) where the modes sit: M = 6
F = ax[1, 2]
for r, wl in enumerate(["(0.000,0.000)", "(0.333,0.667)"]):
    v = o3["b_magnetised"][f"M=6 W={wl}"]
    Wm = np.array([v["weight_at_fixed_points"][str(q)] for q in range(3)])
    F.imshow(Wm, extent=(r * 3.5, r * 3.5 + 3, 3, 0), cmap="Blues", vmin=0, vmax=3.2)
    for q in range(3):
        for f in range(3):
            F.text(r * 3.5 + f + 0.5, q + 0.5, f"{Wm[q, f]:.1f}", ha="center", va="center", fontsize=8)
    F.text(r * 3.5 + 1.5, -0.25, f"W = {wl}\nconteggi {v['counts']}", ha="center", fontsize=8)
F.set_xlim(-0.2, 6.7); F.set_ylim(3.2, -0.9)
F.set_yticks([0.5, 1.5, 2.5]); F.set_yticklabels(["settore 0", "settore 1", "settore 2"])
F.set_xticks([0.5, 1.5, 2.5, 4, 5, 6]); F.set_xticklabels(["P0", "P1", "P2"] * 2)
F.set_title("M = 6: peso dei modi sui 3 punti fissi\nla linea di Wilson sposta i modi da un punto all'altro", fontsize=10)
fig.tight_layout()
fig.savefig(f"{OUT}/orb_figure.png", dpi=110)
print("ok")
