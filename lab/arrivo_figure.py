"""arrivo_figure.py - figure of the arrivo_ study -> lab/out/arrivo_figura.png
(a) the rule in the guide (2-D, spin x): density as soft balls + surfers, at three times
(b) Das-Duerr at L = 100: flux, Bohm spin x (cut-off) and spin z, NR vs rule (hybrid D = 40)
(c) 2-D rule surfers (no spin term by hand), spin x / y / z vs NR, L = 10, omega = 20
(d) detectors at L = 10: flux, ABR kappa = pi, ABR kappa = 2mc/hbar, rule cell detector (1-D D = 20/40/80, 2-D x/y/z)
(e) prediction for the proposed experiment (electron, d = 50 um, L = 5 mm), SI times"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import arrivo_nr as A
import arrivo_rule2d as G
import arrivo_abr as B

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SCR = os.environ.get("ARRIVO_SCRATCH", "/tmp/claude-0/-home-claude-particle-sim/9de435bc-26c2-5d88-b70d-07efb8033e5b/scratchpad")
r3 = json.load(open(os.path.join(OUT, "arrivo_r3.json")))
res = json.load(open(os.path.join(OUT, "arrivo_results.json")))

fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 6, hspace=0.35, wspace=0.6)

# (a) soft balls in the guide
axa = fig.add_subplot(gs[0, 0:3])
D, a, om = 30, 0.5, 20.0
g = G.Guide2D(D, a, om, 10.0, 0.8, "free", zc_dd=20, ell_dd=6)
g.prepare((1, 1))
rng = np.random.default_rng(3)
X0, Y0, Z0 = A.sample_initial(300, om, rng)
Y, Z = Y0 * D, Z0 * D
T = np.full(300, np.nan)
al = np.ones(300, bool)
ts = 0.0
r0, y0f, z0f = g.fields()
snaps = []
for tsnap in (0.0, 1.5, 3.0):
    while ts / D**2 < tsnap - 1e-9:
        r1, y1f, z1f, tm = G.exact_fields(g, 4)
        G._traj2(Y, Z, T, al, ts, tm - ts, a, g.y0, 1e9, r0, y0f, z0f, r1, y1f, z1f)
        r0, y0f, z0f = r1, y1f, z1f
        ts = tm
    rho = r0.copy()
    snaps.append((tsnap, rho, Y.copy() / D, Z.copy() / D))
cols = ["#1f77b4", "#d62728", "#2ca02c"]
for k, (tt, rho, yy, zz) in enumerate(snaps):
    # density as soft balls: sample cells with prob ~ rho, draw translucent discs
    p = rho.ravel() / rho.sum()
    idx = rng.choice(p.size, 1500, p=p)
    iy, iz = np.unravel_index(idx, rho.shape)
    axa.scatter(g.z[iz] / D + rng.normal(0, 0.02, idx.size), g.y[iy] / D + 0.35 * k + rng.normal(0, 0.005, idx.size),
                s=18, alpha=0.06, color=cols[k], edgecolors="none")
    axa.scatter(zz, yy + 0.35 * k, s=4, color="k", alpha=0.7)
    axa.text(10.2, 0.35 * k, "τ = %.1f" % tt, fontsize=9, va="center")
axa.axvline(10, color="gray", lw=1)
axa.set_xlim(-0.2, 11.5)
axa.set_ylim(-0.2, 0.9)
axa.set_xlabel("z / d (guida a piatti, spin lungo x)")
axa.set_yticks([])
axa.set_title("(a) la regola nella guida: densità (palline morbide) e surfisti (punti)", fontsize=10)

# (b) DD L = 100
axb = fig.add_subplot(gs[0, 3:6])
tt = np.array(r3["flux_L100"]["tau"])
axb.plot(tt, r3["flux_L100"]["Pi"], "k-", lw=1, label="flusso Π_qf (NR esatto)")
d = np.load(os.path.join(SCR, "arrivo_r1_NR_A.npz"))
Tx = np.array(res["NR_A"]["x_refined"]["max"])
bins = np.linspace(0, 60, 121)
axb.hist(d["T_alpha=1.5708"], bins=bins, density=True, histtype="step", color="C3", lw=1.5,
         label="Bohm spin x (NR), τ_max = %.1f" % res["NR_A"]["x_refined"]["max"])
axb.hist(d["T_alpha=0.0000"][np.isfinite(d["T_alpha=0.0000"])], bins=bins,
         weights=np.full(np.isfinite(d["T_alpha=0.0000"]).sum(), 1 / d["T_alpha=0.0000"].size / 0.5),
         histtype="step", color="C0", lw=1.5, label="Bohm spin z (NR)")
try:
    Th = np.load(os.path.join(SCR, "r2", "x_L100_D40_a0.125.npz"))["T"]
    axb.hist(Th, bins=bins, density=True, histtype="step", color="C1", ls="--", lw=1.2,
             label="regola (D = 40, a = 0,125) + termine di spin NR")
except FileNotFoundError:
    pass
axb.axvline(42.9, color="C3", ls=":", lw=1)
axb.text(43.3, 0.07, "DD: 42,9", color="C3", fontsize=8)
axb.set_xlabel("τ (unità m d²/ħ)")
axb.set_ylabel("densità del tempo di arrivo")
axb.set_title("(b) Das–Dürr, L = 100 d, ω = 10³", fontsize=10)
axb.legend(fontsize=8)

# (c) 2-D rule surfers
axc = fig.add_subplot(gs[1, 0:2])
dC = np.load(os.path.join(SCR, "arrivo_r1_NR_C.npz"))
bins = np.linspace(0, 8, 41)
axc.hist(dC["T_alpha=1.5708"][:2000], bins=bins, density=False, weights=np.full(2000, 1 / 2000 / 0.2),
         histtype="step", color="C3", lw=1.5, label="NR spin x")
axc.hist(dC["T_alpha=0.0000"][:2000], bins=bins, weights=np.full(2000, 1 / 2000 / 0.2), histtype="step",
         color="C0", lw=1.5, label="NR spin z")
for s, c, ls in (("x", "C3", "--"), ("z", "C0", "--"), ("y", "C2", ":")):
    try:
        T = np.load(os.path.join(SCR, "r4", "free_%s_a0.125.npz" % s))["T"]
        T = T[np.isfinite(T)]
        axc.hist(T, bins=bins, weights=np.full(T.size, 1 / 2000 / 0.2), histtype="step", color=c, ls=ls, lw=1.2,
                 label="regola 2-D spin %s" % s)
    except FileNotFoundError:
        pass
axc.set_xlabel("τ")
axc.set_title("(c) surfisti della regola, L = 10, ω = 20\n(nessun termine di spin aggiunto)", fontsize=10)
axc.legend(fontsize=7)

# (d) detectors
axd = fig.add_subplot(gs[1, 2:4])
fF = r3["flux_L10"]
axd.plot(fF["tau"], fF["F"], "k-", lw=1.5, label="flusso (rivelatore che non riflette)")
ab = r3["ABR_L10"]
axd.plot(ab["pi"]["tau"], ab["pi"]["F"], color="C2", lw=1.2, label="ABR κ = π (Tumulka, Jozani–T.)")
for Dd, c in ((20, "C1"), (40, "C4"), (80, "C5")):
    axd.plot(ab["2D_%d" % Dd]["tau"], ab["2D_%d" % Dd]["F"], color=c, lw=1, ls=":")
    try:
        dd = np.load(os.path.join(SCR, "r2", "det_L10_D%d_a0.0625.npz" % Dd))
        axd.plot(dd["tau"], dd["F"], color=c, lw=1.5, label="cella della regola, d = %d λ_C (… ABR κ = 2mc/ħ)" % Dd)
    except FileNotFoundError:
        pass
for s, ls in (("x", "-"), ("y", "--"), ("z", ":")):
    try:
        dd = np.load(os.path.join(SCR, "r4", "det_%s_a0.25.npz" % s))
        axd.plot(dd["tau"], dd["F"], color="C9", ls=ls, lw=1.2, label="cella 2-D, d = 30 λ_C, spin %s" % s)
    except FileNotFoundError:
        pass
axd.set_xlim(0, 100)
axd.set_xlabel("τ")
axd.set_ylabel("probabilità rivelata entro τ")
axd.set_title("(d) rivelatori, L = 10 d", fontsize=10)
axd.legend(fontsize=6.5, loc="lower right")

# (e) prediction
axe = fig.add_subplot(gs[1, 4:6])
P = r3["PRED"]
tu = P["time_unit_s"]
t100 = np.array(r3["flux_L100"]["tau"])
axe.semilogx(t100 * tu, r3["flux_L100"]["F"], "k-", lw=1.5, label="flusso / ABR ideale (≈ ms)")
Tx = np.sort(d["T_alpha=1.5708"][np.isfinite(d["T_alpha=1.5708"])])
axe.semilogx(Tx * tu, np.arange(1, Tx.size + 1) / Tx.size, color="C3", lw=1.2, label="Bohm spin x (taglio ≈ 0,9 ms)")
mod = B.modal_large_kappa(100.0, P["kappa_DD"], 400000)
ts = np.logspace(-4, 6, 400)
axe.semilogx(ts, B.modal_cdf(mod, ts / tu), color="C5", lw=2, label="cella della regola (κ = 2mc/ħ): specchio")
axe.set_xlabel("tempo dal rilascio [s]")
axe.set_title("(e) previsione: elettrone, d = 50 µm, L = 5 mm", fontsize=10)
axe.legend(fontsize=7, loc="upper left")
axe.set_ylim(0, 1.02)
fig.suptitle("Tempi di arrivo di spin 1/2 sulla regola zig-zag (studio arrivo_)", fontsize=12)
fig.savefig(os.path.join(OUT, "arrivo_figura.png"), dpi=110, bbox_inches="tight")
print("saved")
