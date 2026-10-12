"""knot_figure.py - one figure for K1-K5: lab/out/knot_figure.png"""
import os, sys, json, glob, math, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot_core as kc

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
SCR = os.environ.get("KNOT_SCRATCH", OUT)
J = lambda n: json.load(open(os.path.join(OUT, n)))
fig = plt.figure(figsize=(17, 10))
C = 12 * math.pi**2

# (a) K1
ax = fig.add_subplot(2, 3, 1)
for hh, mk in ((0.2, "o"), (0.15, "s")):
    pts = sorted([(J(f"knot_k1_h{hh:g}_L{L}.json")["half_box"], J(f"knot_k1_h{hh:g}_L{L}.json")["E_over_12pi2"])
                  for L in (32, 48, 64) if os.path.exists(os.path.join(OUT, f"knot_k1_h{hh:g}_L{L}.json"))])
    if pts:
        ax.plot(*zip(*pts), mk + "-", label=f"reticolo h = {hh}")
Rb = np.array([2.4, 3.2, 4, 4.8, 6, 6.4, 8, 10])
if os.path.exists(os.path.join(OUT, "knot_k1_summary.json")):
    S = J("knot_k1_summary.json")
    ax.plot(S["ball_R"], S["ball_E"], "k:", label="ODE in una palla di raggio R")
    ax.axhline(S["continuum"], color="C2", lw=1, label=f"continuo (Richardson) {S['continuum']:.4f}")
ax.axhline(1.2315, color="k", lw=0.8, ls="--", label="letteratura 1.2315")
ax.set_xlabel("mezzo lato della scatola (unità Skyrme)"); ax.set_ylabel("E / 12π²")
ax.set_title("K1  una sola massa: converge come 1/R³ (non ln R)"); ax.legend(fontsize=7)

# (b) K2
ax = fig.add_subplot(2, 3, 2)
vv = np.linspace(0, 0.85, 50)
ax.plot(vv, np.sqrt(1 - vv**2), "k-", lw=0.8, label="1/γ")
for hh, mk in ((0.2, "o"), (0.15, "D")):
    vs, rb, eg, vf = [], [], [], []
    for f in sorted(glob.glob(os.path.join(OUT, f"knot_k2_h{hh:g}_v*.json"))):
        d = json.load(open(f)); vs.append(d["meta"]["v"]); rb.append(d["ratio_b_late"]); eg.append(d["E_over_gE0"])
        vf.append(d["v_fit"] / d["meta"]["v"] if d["meta"]["v"] > 0 else 1)
    if vs:
        ax.plot(vs, rb, mk, color="C0", label=f"h={hh}: σx/σ⊥ barioni (fine corsa)")
        ax.plot(vs, eg, mk, mfc="none", color="C1", label=f"h={hh}: E/(γE0)")
        ax.plot(vs, vf, mk, color="C3", mfc="none", label=f"h={hh}: velocità misurata / v")
ax.set_xlabel("v / c"); ax.set_title("K2  Lorentz in 3D: ok a 0.3; a 0.6-0.8 il reticolo frena", fontsize=9); ax.legend(fontsize=6)

# (c) K3 time series
ax = fig.add_subplot(2, 3, 3)
for f in sorted(glob.glob(os.path.join(OUT, "knot_k3_*.json"))):
    d = json.load(open(f)); r = d["records"]; t = [x["t"] for x in r]
    lab = {"SA": "nodo+antinodo", "SS": "nodo+nodo", "S1": "nodo da solo"}[d["kind"]] + f" v={d['v']}"
    l, = ax.plot(t, [x["E_in"] / r[0]["E_in"] for x in r], label="E dentro r<4.5: " + lab)
    ax.plot(t, [x["absB"] / max(r[0]["absB"], 1e-9) for x in r], "--", color=l.get_color())
ax.set_xlabel("t"); ax.set_ylabel("frazione (— energia, -- ∫|b|)"); ax.set_title("K3  nodo + antinodo → onde"); ax.legend(fontsize=6)

# (d) K3 snapshots, baryon density slice z = 0: top S + antiS (v = 0.6), bottom S + S (v = 0.5, 90-degree scattering)
ax = fig.add_subplot(2, 3, 4)
rows, labs = [], []
for fn, ttl in (("knot_k3_SA_v0.6.json", "nodo+antinodo v=0.6"), ("knot_k3_SS_v0.5.json", "nodo+nodo v=0.5")):
    if os.path.exists(os.path.join(OUT, fn)):
        d = J(fn); ks = [k for k in ("0.0", "4.0", "8.0", "12.0") if k in d["snaps"]]
        ims = [np.array(d["snaps"][k]["b_xy"]).T for k in ks]
        rows.append(np.concatenate([i / max(abs(i).max(), 1e-12) for i in ims], 1)); labs.append(ttl)
if rows:
    ax.imshow(np.concatenate(rows[::-1], 0), cmap="RdBu_r", vmin=-1, vmax=1, origin="lower")
    ax.set_title("K3 densità barionica, piano z=0, t = 0, 4, 8, 12 (ogni quadro normalizzato)" + chr(10) +
                 "sopra: " + labs[0] + " (si annichila)   sotto: " + (labs[1] if len(labs) > 1 else "") + " (esce a 90°)", fontsize=7)
    ax.set_xticks([]); ax.set_yticks([])

# (e) K5 isosurface
ax = fig.add_subplot(2, 3, 5, projection="3d")
tf = os.path.join(SCR, "knot_k5_tris.npy")
if os.path.exists(tf):
    T = np.load(tf); k5 = J("knot_k5.json")
    nrm = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]); nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-12
    lit = 0.35 + 0.65 * np.abs(nrm @ np.array([0.4, -0.5, 0.77]))       # two-sided lambert (triangle orientation is arbitrary)
    fc = np.outer(lit, matplotlib.colors.to_rgb("tab:orange")).clip(0, 1)
    ax.add_collection3d(Poly3DCollection(T, facecolors=fc, edgecolors=(0, 0, 0, 0.15), linewidths=0.2))
    lim = np.abs(T).max() * 1.1
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_title(f"K5 solo il nodo: {k5['triangles']} triangoli per {k5['balls']} palline ({100*k5['ratio']:.2f} %)", fontsize=9)

# (f) K4
ax = fig.add_subplot(2, 3, 6)
if os.path.exists(os.path.join(OUT, "knot_k4.json")):
    k4 = J("knot_k4.json"); p = k4["anw_radial"]; e = k4["experiment"]
    names = ["|μp/μn|", "μp", "|μn|", "g_A", "F_π"]
    th = [p["mu_ratio"], p["mu_p"], -p["mu_n"], p["g_A"], p["F_pi_MeV"]]
    ex = [e["mu_ratio"], e["mu_p"], -e["mu_n"], e["g_A"], e["F_pi_MeV"]]
    ax.bar(range(5), np.array(th) / np.array(ex), color=["C2", "C1", "C1", "C3", "C3"])
    ax.axhline(1, color="k", lw=0.8); ax.set_xticks(range(5)); ax.set_xticklabels(names)
    ax.set_ylabel("Skyrme / misura"); ax.set_title("K4  ANW (ODE esatta): |μp/μn| −3 %, g_A e F_π −30/−50 %")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "knot_figure.png"), dpi=110)
print("saved")
