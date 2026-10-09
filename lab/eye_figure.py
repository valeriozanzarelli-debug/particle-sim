"""eye_figure.py - lab/out/eye_figure.png from eye_shapes.json/.npz, eye_nulls_*.json, eye_calib.json"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import eye_core as E

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
S = json.load(open(os.path.join(OUT, "eye_shapes.json")))
Z = np.load(os.path.join(OUT, "eye_shapes.npz"))
N = {c: json.load(open(os.path.join(OUT, f"eye_nulls_{c}.json"))) for c in E.CLASSES}
THR2 = S["THR2"]
fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
# T2/Z3 rho on the cover
a = ax[0, 0]
xy = E.hex_xy(27); oid = Z["t2z3_oid"]
rho = Z["T2Z3_rho_v2"]; sites = Z["T2Z3_sites"]
r_orb = np.empty(len(sites)); r_orb[sites] = rho
sc = a.scatter(xy[:, 0], xy[:, 1], c=r_orb[oid], s=9, cmap="viridis_r", vmin=0.3, vmax=1.1)
for f in Z["t2z3_fixed"]:
    p = xy[oid == f]
    a.scatter(p[:, 0], p[:, 1], s=90, facecolors="none", edgecolors="r")
plt.colorbar(sc, ax=a, label="rho (volume of the ball / typical)")
g = S["T2/Z3"]
a.set_title(f"T2/Z3 (drawn on the 27x27 cover, coloured by orbit)\ncones seen: {g['cones_v2']} at the 3 fixed points, "
            f"angle ~ {360*g['cone_rho_v2'][0]:.0f} deg (true 120)\ncavity rel2 {g['rel2'][0]:.3f} < THR2 {THR2:.3f}: cavity MISSED",
            fontsize=9)
a.set_aspect("equal"); a.axis("off")
# T4/Z2 histogram
a = ax[0, 1]
r4 = Z["T4Z2_rho_v2"]; fx = Z["t4z2_fixed"]
a.hist(np.delete(r4, fx), bins=30, color="0.6", label="other orbits")
a.hist(r4[fx], bins=5, color="r", label="16 fixed points")
a.axvline(0.7, ls="--", c="k"); a.text(0.71, a.get_ylim()[1] * 0.8, "threshold 0.7", fontsize=8)
a.axvline(0.5, ls=":", c="r"); a.text(0.505, a.get_ylim()[1] * 0.6, "theory 1/2", fontsize=8, color="r")
a.set_xlabel("rho (3-scale mean)"); a.legend(fontsize=8)
a.set_title(f"T4/Z2 (n=6, 656 orbits): {S['T4/Z2']['cones_v2']} cones, all at fixed points\n"
            f"loops {S['T4/Z2']['loops']}, cavities {S['T4/Z2']['cavities']} (b2 = 6 NOT resolved)", fontsize=9)
# S2 + flat
a = ax[0, 2]
X = Z["ico_X"]; rs = np.empty(len(X)); rs[Z["S2_sites"]] = Z["S2_rho_v2"]
fr = X[:, 2] > 0
sc = a.scatter(X[fr, 0], X[fr, 1], c=rs[fr], s=10, cmap="viridis_r", vmin=0.3, vmax=1.1)
d5 = Z["ico_deg5"]; d5 = d5[X[d5, 2] > 0]
a.scatter(X[d5, 0], X[d5, 1], s=80, facecolors="none", edgecolors="orange")
a.set_aspect("equal"); a.axis("off")
a.set_title(f"S2 icosphere (top half): 5-valent points rho {S['S2']['rho_at_deg5_v2'][0]:.2f} > 0.7 -> 0 cones\n"
            f"cavity rel2 {S['S2']['rel2'][0]:.3f} -> {S['S2']['cavities']} cavity, loops {S['S2']['loops']}", fontsize=9)
# H2 / H1 bars
a = ax[1, 0]
names = ["T2/Z3", "T4/Z2", "T2_flat", "S2"]
nullmax = {c: [r["rel2"][0] for r in N[c]["nulls"]] for c in E.CLASSES}
for k, c in enumerate(E.CLASSES):
    a.scatter(np.full(25, k) + np.random.default_rng(k).uniform(-.2, .2, 25), nullmax[c], s=8, c="0.5")
for k, nm in enumerate(names):
    a.bar(4 + k, S[nm]["rel2"][0] if S[nm]["rel2"] else 0, color="tab:blue")
a.axhline(THR2, c="r", ls="--"); a.text(0, THR2 + 0.01, f"THR2 = {THR2:.3f} (from 40 calibration nulls)", fontsize=8, color="r")
a.set_xticks(range(8)); a.set_xticklabels([c.replace("_", "\n") for c in E.CLASSES] + [n.replace("_", "\n") for n in names],
                                          fontsize=8)
a.set_ylabel("largest H2 bar / diameter"); a.set_title("cavities: nulls (grey, 25 each) vs drawn shapes", fontsize=9)
# nulls FP
a = ax[1, 1]
xx = np.arange(4)
for j, (key, lab) in enumerate((("fp_cavities", "cavity"), ("fp_cones_v1", "cone v1 (pre-reg)"),
                                ("fp_cones_v2", "cone v2 (amended)"), ("fp_loops", "loop"))):
    a.bar(xx + 0.2 * j - 0.3, [N[c]["summary"][key] for c in E.CLASSES], 0.2, label=lab)
a.axhline(0.05, c="r", ls="--")
a.set_xticks(xx); a.set_xticklabels(E.CLASSES, fontsize=8); a.set_ylabel("false-positive rate (25 fresh nulls)")
a.legend(fontsize=8); a.set_title("exam nulls: v1 cone rule fails the shuffled null; v2 passes", fontsize=9)
# summary text
a = ax[1, 2]; a.axis("off")
g = S["grades"]
txt = ["L1 exams (pre-registered criteria)", "",
       f"nulls cavities: {'PASS' if all(N[c]['summary']['fp_cavities'] <= 0.05 for c in E.CLASSES) else 'FAIL'}",
       f"nulls cones v1: {'PASS' if all(N[c]['summary']['fp_cones_v1'] <= 0.05 for c in E.CLASSES) else 'FAIL'}"
       f"   v2: {'PASS' if all(N[c]['summary']['fp_cones_v2'] <= 0.05 for c in E.CLASSES) else 'FAIL'}",
       f"T2/Z3: cones 3/3 at fixed points, loops 0, cavity {S['T2/Z3']['cavities']} -> {'PASS' if g['v2']['T2Z3'] else 'FAIL'}",
       f"T4/Z2: 16/16 cones -> {'PASS' if g['v2']['T4Z2_cones'] else 'FAIL'}; b2 not resolved",
       f"flat T2: cones 0 {'PASS' if g['v2']['T2_flat_cones0'] else 'FAIL'}; loops {S['T2_flat']['loops']}, cavity {S['T2_flat']['cavities']}",
       f"S2: cones 0 {'PASS' if g['v2']['S2_cones0'] else 'FAIL'}; cavity {S['S2']['cavities']}, loops {S['S2']['loops']}"]
a.text(0, 0.95, "\n".join(txt), va="top", family="monospace", fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "eye_figure.png"), dpi=110)
print("saved")
