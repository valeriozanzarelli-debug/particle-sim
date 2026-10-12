"""glue_figure.py - lab/out/glue_figure.png from glue_results.json / glue_verdict.json."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
res = json.load(open(os.path.join(OUT, "glue_results.json")))["results"]
V = json.load(open(os.path.join(OUT, "glue_verdict.json")))
betas = sorted(res.keys(), key=float)
col = {"0++": "#1f77b4", "0++*": "#2ca02c", "2++": "#d62728", "torelon": "#7f7f7f"}
bcol = {"6": "#9467bd", "9": "#ff7f0e", "12": "#17becf"}
lit = {"0++": 4.7367, "0++*": 6.861, "2++": 7.762}

fig, ax = plt.subplots(2, 2, figsize=(12, 9.5))

# (a) effective masses at the finest beta
fb = betas[-1]; r = res[fb]
a = ax[0, 0]
for lab, st in (("0++", r["A1"][0]), ("0++*", r["A1"][1]), ("2++", r["E"][0])):
    me = np.array([np.nan if x is None else x for x in st["meff"]]); de = np.array([np.nan if x is None else x for x in st["meff_err"]])
    t = np.arange(len(me))
    ok = np.isfinite(me) & np.isfinite(de) & (de < 0.5)
    a.errorbar(t[ok] + 0.5, me[ok], de[ok], fmt="o", color=col[lab], capsize=2, label=f"{lab}")
    a.fill_between([st["t_min"], st["t_max"] + 1], st["m"] - st["err"], st["m"] + st["err"], color=col[lab], alpha=0.25)
    s = res[fb]["a_sqrt_sigma"]
    a.axhline(lit[lab] * s, color=col[lab], ls=":", lw=1)
a.set_xlabel("t / a  (punto a t+1/2)"); a.set_ylabel("a m_eff(t)")
a.set_title(f"(a) masse efficaci a beta = {fb} (bande: fit con regola fissa; puntini: letteratura)")
a.set_ylim(0, 1.6); a.set_xlim(0, 9); a.legend(fontsize=9); a.title.set_fontsize(9)

# (b) spectrum m/sqrt(sigma) vs (a sqrt sigma)^2
a = ax[0, 1]
key = {"0++": "m0_s", "0++*": "m0s_s", "2++": "m2_s"}
for lab, k in key.items():
    x = [V["per_beta"][b]["a_sqrt_sigma"][0] ** 2 for b in betas]
    y = [V["per_beta"][b][k][0] for b in betas]; e = [V["per_beta"][b][k][1] for b in betas]
    a.errorbar(x, y, e, fmt="o", color=col[lab], capsize=3, label=lab)
    if "continuum" in V:
        c, ec, sl, _ = V["continuum"][k]
        xx = np.linspace(0, max(x) * 1.05, 20)
        a.plot(xx, c + sl * xx, color=col[lab], lw=1)
        a.errorbar([-0.002], [c], [ec], fmt="s", color=col[lab], mfc="white", capsize=3)
    a.axhline(lit[lab], color=col[lab], ls=":", lw=1)
a.set_xlabel("(a sqrt(sigma))^2   [0 = continuo]"); a.set_ylabel("m / sqrt(sigma)")
a.set_title("(b) spettro (quadrati vuoti: continuo) vs Athenodorou-Teper 2017 (puntini)", fontsize=9)
a.legend(fontsize=9)

# (c) Creutz ratios
a = ax[1, 0]
for b in betas:
    cr = res[b]["creutz"]
    R = np.array(sorted(int(k) for k in cr)); v = np.array([cr[str(k)][0] for k in R]); e = np.array([cr[str(k)][1] for k in R])
    ok = e < 0.5 * np.abs(v)
    sig = res[b]["a_sqrt_sigma"] ** 2
    a.errorbar(R[ok], v[ok] / sig, e[ok] / sig, fmt="o-", color=bcol[b], capsize=2, label=f"beta {b}: chi(R,R)/(sigma a^2)")
    rr = np.linspace(1.8, 6.2, 50)
    a.plot(rr, 1 + np.pi / (24 * (rr - 0.5) ** 2) / sig, color=bcol[b], ls="--", lw=1)
a.axhline(1, color="k", lw=0.8)
a.set_xlabel("R / a"); a.set_ylabel("rapporto di Creutz / tensione di stringa (toreloni)")
a.set_ylim(0, 3); a.legend(fontsize=8)
a.set_title("(c) confinamento: rapporti di Creutz -> sigma (tratteggio: sigma + pi/(24 (R-1/2)^2))", fontsize=9)

# (d) G3: beta a sqrt(sigma) vs 1/beta, and G4 dispersion
a = ax[1, 1]
if "G3" in V:
    bb = np.array([float(b) for b in betas]); y = np.array(V["G3"]["beta_a_sqrt_sigma"]); e = np.array(V["G3"]["err"])
    a.errorbar(1 / bb, y / 4, e / 4, fmt="o", color="k", capsize=3, label="noi: sqrt(sigma)/g^2 = beta a sqrt(sigma) / 4")
    xx = np.linspace(0, 1 / bb.min() * 1.05, 20)
    a.plot(xx, (V["G3"]["sqrt_sigma_over_g2"] * 4 + V["G3"]["slope"] * xx) / 4, "k-", lw=1)
    a.axhline(0.3353, color="r", ls=":", label="Teper 1998: 0.3353")
a.set_xlabel("1 / beta = a g^2 / 4"); a.set_ylabel("sqrt(sigma) / g^2")
txt = "\n".join([f"beta {b}: E(p)^2 - m^2 = c^2 p^2, c^2 = {V['per_beta'][b]['c2'][0]:.2f} +- {V['per_beta'][b]['c2'][1]:.2f}"
                 for b in betas])
a.text(0.02, 0.05, "glueball in moto (G4):\n" + txt, transform=a.transAxes, fontsize=8, va="bottom")
a.legend(fontsize=8, loc="upper left")
a.set_title("(d) l'unica scala e' g^2 (nessuna massa nella regola)", fontsize=9)

fig.suptitle("Glueball in SU(2) 2+1 D: materia fatta solo di 'luce' che si attacca a se stessa", fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "glue_figure.png"), dpi=130)
print("saved")
