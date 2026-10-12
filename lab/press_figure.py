"""press_figure.py - one figure for P1-P5 of the pressure study: lab/out/press_figure.png"""
import os, json, glob, math, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
J = lambda n: json.load(open(os.path.join(OUT, n)))
O = np.load(os.path.join(OUT, "press_profiles.npz")); R = J("press_ode.json")
HBARC = 0.1973269804
ell, pu, eps0 = float(O["ell"]), float(O["pu"]), float(O["eps0"])
r = O["r"]; m = r < 6
rf = r[m] * ell                                  # fm
p_ph, s_ph, T_ph = O["p"][m] * pu, O["s"][m] * pu, O["T00"][m] * pu
D0, Mt, mN = float(O["D0"]), float(O["Mt"]), float(O["mN"])
x = np.linspace(1e-3, 2.5, 600); u = Mt * x / HBARC
pbeg = D0 * Mt**5 * (u - 3) * np.exp(-u) / (192 * np.pi * mN) / HBARC**3      # GeV/fm^3
fig, axs = plt.subplots(2, 2, figsize=(13, 9.5))

ax = axs[0, 0]
ax.plot(rf, rf**2 * p_ph, "C0-", lw=2, label=f"nodo di Skyrme (ANW: e=5.44, F_π=129 MeV, FIT)  D={R['P2']['D']:.2f}")
ax.plot(x, x**2 * pbeg, "C3--", lw=2, label=f"protone, JLab (Burkert+ 2018, tripolo)  D={D0:.2f}")
ax.axhline(0, color="k", lw=0.6)
ax.axvline(R["P2"]["r0_fm"], color="C0", ls=":", lw=1); ax.axvline(R["P2"]["BEG"]["r0_fm"], color="C3", ls=":", lw=1)
ax.set_xlim(0, 2.0); ax.set_xlabel("r (fm)"); ax.set_ylabel("r² p(r)  (GeV/fm)")
ax.set_title("P2  la pressione dentro: spinge al centro, tira sul bordo\n"
             f"zero: nodo {R['P2']['r0_fm']:.2f} fm, protone {R['P2']['BEG']['r0_fm']:.2f} fm;  "
             f"p(0): {R['P2']['p0_Pa']:.1e} vs {R['P2']['BEG']['p0_Pa']:.1e} Pa", fontsize=10)
ax.legend(fontsize=8)

ax = axs[0, 1]
ax.semilogy(rf, T_ph, "k-", label="T$_{00}$ densità di energia (massa)")
ax.semilogy(rf, np.where(p_ph > 0, p_ph, np.nan), "C0-", label="p > 0 (spinge)")
ax.semilogy(rf, np.where(p_ph < 0, -p_ph, np.nan), "C0--", label="−p (p < 0, tira)")
ax.semilogy(rf, s_ph, "C2-", label="s  forza di taglio")
ax.set_xlim(0, 3); ax.set_ylim(1e-4, 5); ax.set_xlabel("r (fm)"); ax.set_ylabel("GeV/fm³")
ax.set_title(f"P1  il nodo esatto (ODE): ∫T00 = M a {R['P1']['a_rel']:.0e},  von Laue ∫r²p = 0 a {R['P1']['vonLaue_rel']:.0e}\n"
             f"D da p e da s: diff. {R['P1']['D_rel_diff']:.0e};  Cebulla 2007 Tab.2 riprodotta entro 0.2 %", fontsize=10)
ax.legend(fontsize=8)

ax = axs[1, 0]
P3 = R["P3"]
lab1 = ["nodo: σ\n(tira)", "nodo: Skyrme\n(spinge)"]
ax.bar([0, 1], [100 * P3["frac_E2"], 100 * P3["frac_E4"]], color=["C0", "C1"])
q = J("press_prereg.json")["references"]["MassDecomp"]["percent"]
names = {"quark_condensate": "massa\nquark", "quark_energy": "energia\nquark", "glue_field_energy": "campo\ngluoni",
         "trace_anomaly": "anomalia\ndi traccia"}
xs = np.arange(3, 3 + len(q))
ax.bar(xs, list(q.values()), color="C7")
ax.set_xticks([0, 1] + list(xs)); ax.set_xticklabels(lab1 + [names[k] for k in q], fontsize=7.5)
for xi, v in zip([0, 1] + list(xs), [100 * P3["frac_E2"], 100 * P3["frac_E4"]] + list(q.values())):
    ax.text(xi, v + 1, f"{v:.0f}%", ha="center", fontsize=9)
ax.axvline(2, color="k", lw=0.6)
ax.set_ylabel("% della massa"); ax.set_ylim(0, 62)
ax.set_title("P3  da dove viene la massa: nodo (Derrick, esatto 50/50) | protone (χQCD 2018, reticolo)\n"
             "confronto solo qualitativo, nessun fit", fontsize=10)

ax = axs[1, 1]
Lt = np.load(os.path.join(OUT, "press_lattice.npz")); LJ = J("press_lattice.json")
ax.plot(r[m], O["p"][m], "k-", lw=1.5, label="ODE esatta")
ax.plot(Lt["r"], Lt["p"], "C3o", ms=4, label=f"reticolo h=0.2, L=48 (derivate centrali, pre-reg.)")
for f in sorted(glob.glob(os.path.join(OUT, "press_posthoc_h*_L*.npz"))):
    d = np.load(f)
    ax.plot(d["r"], d["p_fb"], "s", ms=3, mfc="none", label=f"post-hoc h={float(d['h']):g}: p con le differenze dell'energia")
ax.axhline(0, color="k", lw=0.6); ax.axvline(float(O["r0"]) / 2, color="C7", ls=":", lw=1)
ax.set_xlim(0, 4.7); ax.set_xlabel("r (unità di Skyrme)"); ax.set_ylabel("p (unità di Skyrme)")
ax.set_title(f"P4  il nodo di palline: nucleo (r < r₀/2) off di {100*LJ['core_max_dev_over_p0']:.0f}% di p(0)  → FAIL\n"
             f"r₀ reticolo {LJ['r0_lattice']:.3f} vs ODE {LJ['r0_ode']:.3f}", fontsize=10)
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "press_figure.png"), dpi=110)
print("wrote press_figure.png")
