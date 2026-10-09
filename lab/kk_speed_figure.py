"""kk_speed_figure.py - lab/out/kk_speed.png from kk_speed_kk.json, kk_speed_twist.json/_tmp.npz, kk_speed_flow.json"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = json.load(open("out/kk_speed_kk.json"))
B = json.load(open("out/kk_speed_twist.json"))
Bt = dict(np.load("out/kk_speed_twist_tmp.npz"))
C = json.load(open("out/kk_speed_flow.json"))
fig, ax = plt.subplots(2, 3, figsize=(15, 8.6), dpi=80)

# (a) walk rest energy vs c_w
a = ax[0, 0]
rb = A["b_rest_energy_vs_cw"]
cw = np.array([r["cw"] for r in rb]); E = np.array([r["E1"] for r in rb])
ff = np.geomspace(0.01, 1.2, 200); R = 32 / (2 * np.pi)
a.loglog(ff, np.sqrt(0.02**2 + (ff / R)**2), "k-", lw=1, label="sqrt(theta^2 + (c_w/R)^2)")
a.loglog(cw, E, "o", color="tab:red", label="walk, n = 1 (N = 32, theta = 0.02)")
a.axhline(0.02, color="gray", ls=":", lw=1, label="flip mass theta")
a.set_xlabel("hop speed along the hidden circle  c_w / c"); a.set_ylabel("rest energy of n = 1 (lattice units)")
a.set_title("1a. slower spin on the circle = lighter n = 1\n(and (N, c_w) gives exactly the spectrum of (N/c_w, 1))", fontsize=10)
a.legend(fontsize=8)

# (b) KK: q^2/(G m^2) vs c_w
a = ax[0, 1]
rc = A["c_kk_numbers"]["rows"]
x = np.array([r["cw"] for r in rc]); m1 = np.array([r["m1_over_me"] for r in rc]); al = np.array([r["alpha_KK"] for r in rc])
a.loglog(x, m1, "o-", color="tab:red", label="mass of n = 1 / m_e")
a.loglog(x, al / 7.2973525693e-3, "s-", color="tab:blue", label="KK charge^2 / e^2  (alpha_KK / alpha)")
a.loglog(x, 4 + 0 * x, "-", color="k", lw=2, label="q^2 / (G m^2) = 4 for every c_w")
a.axhline(A["c_kk_numbers"]["electron_q2_over_Gm2"], color="tab:green", ls="--", label="electron: e^2/(G m_e^2) = 4.2e42")
a.set_xlabel("c_w / c  (gauge coupling fixed by gravity, KK)"); a.set_ylim(1e-46, 1e46)
a.set_title("1b. pure KK: mass and charge fall TOGETHER\nq/m never moves (c_w = the radion)", fontsize=10)
a.legend(fontsize=7.5, loc=(0.02, 0.5))

# (c) dynamic twist with its own K: momentum and position of the left packet
a = ax[0, 2]
cols = {"1_1_0.05_1000.0": ("tab:blue", "+1,+1  c_w=0.05 K=1000"), "1_1_0.1_1000.0": ("tab:red", "+1,+1  c_w=0.10 K=1000"),
        "1_1_0.2_1000.0": ("tab:purple", "+1,+1  c_w=0.20 K=1000"), "1_1_0.1_4000.0": ("tab:orange", "+1,+1  c_w=0.10 K=4000"),
        "1_-1_0.1_1000.0": ("tab:green", "+1,-1  c_w=0.10 K=1000")}
a2 = a.twinx()
for key, (col, lab) in cols.items():
    tr = Bt[key]; t = np.arange(tr.shape[0])
    a.plot(t, tr[:, 1], color=col, lw=2, label=lab)
    a2.plot(t, tr[:, 0] - tr[0, 0], color=col, lw=1, ls="--")
a.set_xlabel("tick"); a.set_ylabel("momentum <k> of the left packet (solid)"); a2.set_ylabel("displacement <x> - x0 (dashed)")
a.set_title("2. own-stiffness twist: force set by K only,\ninertia set by c_w only (solid lines overlap)", fontsize=10)
a.legend(fontsize=7.5, loc="lower left")

# (d) 1+1 D chains
a = ax[1, 0]
i = C["i_chains"]
for g, col in ((0.3, "tab:red"), (0.8, "tab:purple")):
    rows = [r for r in i["A_gradient"] if r["g"] == g]
    b = [r["block"] for r in rows]
    a.semilogx(b, [r["v_slow"] for r in rows], "o-", color=col, ms=3, label=f"A: gradient coupling g = {g} (two speeds)")
    a.semilogx(b, [r["v_fast"] for r in rows], "o-", color=col, ms=3)
rows = i["B_mass_mixing"]
a.semilogx([r["block"] for r in rows], [r["v_low"] for r in rows], "s-", color="tab:blue", ms=3, label="B: mass mixing mu = 0.05 (one low band)")
a.axhline(i["B_common_speed"], color="tab:blue", ls=":", lw=1)
a.axhline(1.0, color="gray", lw=0.5); a.axhline(0.5, color="gray", lw=0.5)
a.set_xlabel("block size b (coarse-graining)"); a.set_ylabel("speed at k = pi / b  (c1 = 1, c2 = 0.5)")
a.set_title("3i. 1+1 D: gradient coupling -> speeds never merge;\nmass mixing -> one mode left, at the AVERAGE speed", fontsize=10)
a.legend(fontsize=7.5)

# (e) 3+1 D Yukawa flow
a = ax[1, 1]
for f, col in zip(C["ii_flows"], plt.cm.viridis(np.linspace(0, 0.9, len(C["ii_flows"])))):
    l = np.array(f["l"]); d = np.abs(np.array(f["delta"]))
    a.loglog(l[1:], d[1:], color=col, lw=1.3, label=f"c_s/c_f = {f['cs0']}, g^2 = {f['g2_0']}")
a.axvline(51.5, color="k", ls="--", lw=1); a.text(55, 1.5e-3, "Planck -> m_e\n(51.5 e-folds)", fontsize=8)
a.set_xlabel("l = ln(Lambda / mu)  (towards the infrared)"); a.set_ylabel("|c_s / c_f - 1|")
a.set_ylim(1e-3, 2)
a.set_title("3ii. 3+1 D Yukawa, one loop: speeds DO merge,\nbut only as l^-0.47 (power of a logarithm)", fontsize=10)
a.legend(fontsize=6.5, ncol=1, loc="lower left")

# (f) text summary
a = ax[1, 2]; a.axis("off")
p = C["ii_linear_exponent"]["p"]
r0 = B["runs"][0]
txt = ("Numbers\n\n"
       "1. Pure KK, any c_w:  q^2/(G m^2) = 4.000\n"
       "   electron: 4.2e42.  Light electron (c_w = 1e-21 c)\n"
       "   -> alpha_KK = 7e-45.  c_w is R in disguise.\n\n"
       "2. Own stiffness K (not gravity):\n"
       f"   force  = n1 n2 / (2 K R^2) to {max(abs(r['F_rel_err']) for r in B['runs']):.0e}\n"
       "   inertia = rest energy to 2-6 %\n"
       "   -> q and m decouple, but alpha = 1/137\n"
       "   is now put in by hand (K) and the\n"
       "   needed K is 1e-42 x gravity's.\n"
       "   n = 2 (charge 2e, mass 2 m_e): not seen.\n\n"
       f"3. Speeds merge only as l^-{p:.2f}:\n"
       "   Planck -> m_e shrinks delta by x2-16;\n"
       f"   bound 1e-15 needs ~10^{15/p:.0f} e-folds.\n"
       "   On the compact circle below 1/R there is\n"
       "   no speed left to flow: c_w is a MASS,\n"
       f"   QED changes it by x{C['qed_mass_running_factor']:.2f}, not 1e21.")
a.text(0.0, 1.0, txt, va="top", family="monospace", fontsize=9.2)
fig.tight_layout()
fig.savefig("out/kk_speed.png")
print("wrote out/kk_speed.png")
