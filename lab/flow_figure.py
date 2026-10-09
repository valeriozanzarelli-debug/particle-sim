"""flow_figure.py - lab/out/flow_figure.png from flow_{a,b1,b2,c1,c2,z3}.json/.npz"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
CASES = ["a", "b1", "b2", "c1", "c2", "z3"]
R = {c: json.load(open(os.path.join(OUT, f"flow_{c}.json"))) for c in CASES}
Z = {c: np.load(os.path.join(OUT, f"flow_{c}.npz")) for c in CASES}


def vmean(z, w):
    N = len(z["pos"]); s = np.zeros(N); n = np.zeros(N)
    np.add.at(s, z["E"][:, 0], w); np.add.at(s, z["E"][:, 1], w); np.add.at(n, z["E"][:, 0], 1); np.add.at(n, z["E"][:, 1], 1)
    return s / n


def scat(a, z, val, title, cmap="coolwarm", vmin=None, vmax=None):
    p = z["pos"]
    if p.shape[1] == 3:
        fr = p[:, 2] > 0; p = p[fr][:, :2]; val = val[fr]
    sc = a.scatter(p[:, 0], p[:, 1], c=val, s=22, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar(sc, ax=a, fraction=0.046)
    a.set_aspect("equal"); a.axis("off"); a.set_title(title, fontsize=9)


fig, ax = plt.subplots(3, 4, figsize=(19, 13))
z = Z["a"]
scat(ax[0, 0], z, np.log10(vmean(z, z["w0"])), "(a) hex torus: log10 mean link length BEFORE (40% noise)", "viridis", -1.5, 0.6)
scat(ax[0, 1], z, np.log10(vmean(z, z["w"])), f"(a) AFTER {R['a']['steps']} steps (stop: {R['a']['stop']}): lengths split\n"
     f"(link-length spread {R['a']['w_cv0']:.2f} -> {R['a']['w_cv1']:.2f})", "viridis", -1.5, 0.6)
a = ax[0, 2]
for c in CASES:
    h = R[c]["hist"]
    a.semilogy([q["it"] for q in h], [q["k_std"] for q in h], label=f"{c} ({R[c]['kind']})")
a.set_xlabel("flow step"); a.set_ylabel("spread of Ollivier curvature (std over links)"); a.legend(fontsize=7)
a.set_title("curvature spread during the flow (pre-reg: must fall 10x for a, b)", fontsize=9)
a = ax[0, 3]
for c in CASES:
    h = R[c]["hist"]
    a.plot([q["it"] for q in h], [q["k_mean"] for q in h], label=c)
a.axhline(0, c="k", lw=0.5); a.set_xlabel("flow step"); a.set_ylabel("mean Ollivier curvature"); a.legend(fontsize=7)
a.set_title("mean curvature: sphere (c) vs torus (b)", fontsize=9)
z = Z["z3"]
lim = max(np.nanmax(np.abs(z["d0"])), 1e-3)
scat(ax[1, 0], z, z["d0"], "Z3 pattern: Regge deficit BEFORE (rad)", "coolwarm", -lim, lim)
scat(ax[1, 1], z, z["d1"], f"Z3 pattern: deficit AFTER  -> verdict: {R['z3']['z3']['verdict']}", "coolwarm", -lim, lim)
lk = max(np.abs(z["vc0"]).max(), 1e-3)
scat(ax[1, 2], z, z["vc0"], "Z3 pattern: Ollivier curvature per vertex BEFORE", "coolwarm", -lk, lk)
scat(ax[1, 3], z, z["vc1"], "Z3 pattern: Ollivier curvature per vertex AFTER", "coolwarm", -lk, lk)
for a in ax[1]:
    p = z["pos"]
    for f in (0, 5 * 15 + 5, 10 * 15 + 10):
        a.scatter(p[f, 0], p[f, 1], s=120, facecolors="none", edgecolors="k")
zb = Z["b1"]
scat(ax[2, 0], zb, zb["d1"], "(b1) Delaunay torus: Regge deficit AFTER", "coolwarm",
     -np.nanmax(np.abs(zb["d1"])), np.nanmax(np.abs(zb["d1"])))
zc = Z["c1"]
scat(ax[2, 1], zc, zc["vc1"], "(c1) sphere (top half): Ollivier curvature AFTER", "viridis")
a = ax[2, 2]
x = np.arange(len(CASES))
a.bar(x - 0.2, [R[c]["hol0"]["frac_vertices_pass"] for c in CASES], 0.4, label="before")
a.bar(x + 0.2, [R[c]["hol1"]["frac_vertices_pass"] for c in CASES], 0.4, label="after")
a.set_xticks(x); a.set_xticklabels(CASES); a.set_ylim(0, 1.05); a.legend(fontsize=8)
a.set_ylabel("fraction of vertices where the spinor comes back (|deficit| < 0.05)")
a.set_title("invariant spinors of the whole surface: 0 in every case, before and after\n(2 would need ALL vertices = 1.0)", fontsize=9)
a = ax[2, 3]; a.axis("off")
lines = ["case  steps  spread0 -> spread1 (x)   eye after: loops cav cones(v1/v2)   N_inv   max|deficit|"]
for c in CASES:
    r = R[c]; e = r["eye1"]
    lines.append(f"{c:4s} {r['steps']:5d}  {r['k_std0']:.4f} -> {r['k_std1']:.4f} ({r['spread_ratio']:.1f})   "
                 f"{e['loops']}  {e['cavities']}  {e['cones']}/{e['cones_v2']}   {r['hol1']['N_inv']}   {r['hol1']['max_abs_deficit']:.3f}")
lines += ["", f"z3: C_deficit {R['z3']['z3']['C_deficit0']:.2f} -> {R['z3']['z3']['C_deficit1']:.2f} "
              f"(area {R['z3']['z3']['area_baseline']:.2f}); max|def| {R['z3']['z3']['maxdef0']:.3f} -> {R['z3']['z3']['maxdef1']:.3f}"]
a.text(0, 0.95, "\n".join(lines), va="top", family="monospace", fontsize=7.5)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "flow_figure.png"), dpi=100)
print("saved")
