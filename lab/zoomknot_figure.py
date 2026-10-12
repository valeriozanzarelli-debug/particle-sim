"""zoomknot_figure.py - lab/out/zoomknot_figure.png: (a) Z2 nested grid, z = 0 cross-section coloured by pressure,
(b) p(r) uniform vs nested vs agglomerated vs exact ODE, (c) error vs cost, (d) agglomerate error vs block size."""
import os, sys, json, math, numpy as np, torch
import matplotlib; matplotlib.use("Agg"); import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zoomknot_core as zc

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
PRE = json.load(open(os.path.join(OUT, "zoomknot_prereg.json")))
RES = json.load(open(os.path.join(OUT, "zoomknot_results.json")))
L = json.load(open(os.path.join(OUT, "zoomknot_learn.json")))
O = np.load(os.path.join(OUT, "press_profiles.npz"))
ld = lambda t: json.load(open(os.path.join(OUT, f"zoomknot_{t}.json")))

fig, axs = plt.subplots(2, 2, figsize=(13, 11))
# (a) cross-section
ax = axs[0, 0]
H = zc.Hier(PRE["configs"]["Z2"], L["k2"], L["k4"])
S = torch.load(os.path.join(OUT, "zoomknot_ckpt", "Z2_learned.pt"))
with torch.no_grad():
    full = H.fields(S["us"])
    lim = 3.0
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(H.L)))
    for l, lv in enumerate(H.L):
        h, N = lv["h"], lv["N"]
        T, _, _ = zc.stress_central(full[l], h)
        p = ((T[0, 0] + T[1, 1] + T[2, 2]) / 3)[:, :, N].numpy()
        x = zc.ax(N, h).numpy()
        W = H.W[l][:, :, N].numpy()
        X, Y = np.meshgrid(x, x, indexing="ij")
        m = (W > 0) & (np.abs(X) <= lim) & (np.abs(Y) <= lim)
        sc = ax.scatter(X[m], Y[m], c=p[m], s=max(1.5, 120 * h * h), cmap="RdBu_r", vmin=-4, vmax=12, marker="s",
                        linewidths=0)
        if lv["R"] is not None and lv["R"] <= lim:
            ax.add_patch(Rectangle((-lv["R"], -lv["R"]), 2 * lv["R"], 2 * lv["R"], fill=False, ec="k", lw=0.8))
            ax.text(-lv["R"] + 0.03, lv["R"] - 0.13, f"h={h:g}", fontsize=7)
    plt.colorbar(sc, ax=ax, label="p (Skyrme units)", shrink=0.8)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
ax.set_title("(a) Z2: z = 0 slice, |x|<3 (8 levels h = 0.05...6.4 out to 89.6)\neach square = one ball / agglomerate, colour = pressure")
ax.set_xlabel("x"); ax.set_ylabel("y")
# (b) p(r)
ax = axs[0, 1]
r = O["r"]; ax.plot(r[r < 3], O["p"][r < 3], "k-", lw=2, label="exact ODE")
best = RES["verdicts"].get("Z1c", {}).get("best", "Z1_h0.05_R1.0")
for t, st in (("U0", "C3o-"), (best, "C0s-"), ("Z2", "C1^-"), ("Z2_learned", "C2d-")):
    if not os.path.exists(os.path.join(OUT, f"zoomknot_{t}.json")):
        continue
    m = ld(t); P = m["profile"]
    rr = np.array(P["r"]); pp = np.array(P["p_central"]); s = rr < 3
    ax.plot(rr[s], pp[s], st, ms=3, lw=0.8, label=f"{t}: p(0) = {m['p0_central']:.2f}, {m['sites']} sites")
ax.axhline(0, color="gray", lw=0.5); ax.set_xlabel("r"); ax.set_ylabel("p(r) (shell average, 4th-order central)")
ax.set_title("(b) pressure profile"); ax.legend(fontsize=8)
# (c) error vs cost
ax = axs[1, 0]
for row in RES["rows"]:
    c = "C3" if row["tag"] == "U0" else ("C2" if row["tag"].startswith("Z2") else "C0")
    ax.loglog(row["sites"], abs(row["err_p0"]), "o", color=c, ms=8)
    ax.loglog(row["sites"], abs(row["err_E"]), "s", color=c, mfc="none", ms=8)
    ax.annotate(row["tag"].replace("_h", " h").replace("_R", " R"), (row["sites"], abs(row["err_p0"])), fontsize=7,
                xytext=(4, 4), textcoords="offset points")
ax.axvline(64**3, color="k", ls="--", lw=0.8); ax.text(64**3 * 1.02, 0.3, "64^3 budget", fontsize=8)
ax.axhline(0.05, color="C1", ls=":", lw=0.8); ax.axhline(0.003, color="C4", ls=":", lw=0.8)
ax.text(1.05e5, 0.055, "5% (p0 target)", fontsize=8); ax.text(1.05e5, 0.0033, "0.3% (E target)", fontsize=8)
ax.set_xlim(1e5, 4e5); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks([1.2e5, 1.6e5, 2e5, 2.62e5, 3.4e5]); ax.set_xticklabels(["120k", "160k", "200k", "262k", "340k"])
ax.set_xlabel("sites (all levels, incl. ghosts)"); ax.set_ylabel("|relative error|")
ax.set_title("(c) error vs cost: filled = p(0), open = E\nred = uniform h=0.2, blue = nested (Z1), green = agglomerates (Z2)")
# (d) agglomerate error vs block size
ax = axs[1, 1]
pl = L["per_level"]
b = [d["block_in_finest_balls"] for d in pl]
ax.semilogx(b, [100 * d["err_bare"] for d in pl], "o-", color="C3", label="bare lattice coupling")
ax.semilogx(b, [100 * d["err_learned"] for d in pl], "s-", color="C2", label="learned coupling (held-out lam = 1)")
ax.semilogx(b, [100 * d["err_fine_h4"] for d in pl], "^:", color="gray", label="fine truth used for learning (h/4)")
ax.axhline(0, color="k", lw=0.5)
for d in pl:
    ax.annotate(f"k2={d['k2']:.3f}", (d["block_in_finest_balls"], 100 * d["err_learned"]), fontsize=7,
                xytext=(3, -10), textcoords="offset points")
ax.set_xlabel("agglomerate side in finest balls (h_l / 0.05)"); ax.set_ylabel("static energy error of the shell (%)")
ax.set_title("(d) check: error introduced by each agglomerate level")
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "zoomknot_figure.png"), dpi=110)
print("saved")
