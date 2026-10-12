"""h2vib_figure.py - figure of the H2 vibration study -> lab/out/h2vib_figure.png"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from h2vib_run import load_points, Curva, EXACT, OUT
from h2vib_fix import MODELS

r = json.load(open(os.path.join(OUT, "h2vib_run.json")))
fx = json.load(open(os.path.join(OUT, "h2vib_fix.json")))
pts, runs, kv = load_points()
e_inf = 2 * kv["h_E_guided"]
ex = np.array([(R, float(EXACT(R)), de) for R, _, de in pts])
fig, ax = plt.subplots(2, 2, figsize=(12, 9))

a = ax[0, 0]
rr = np.linspace(0.85, 3.0, 400)
c5 = MODELS["rust_nb5_grid"](ex, -1.0); c6 = MODELS["nb6_grid"](ex, -1.0)
a.plot(rr, 1e4 * (c5.v(rr) - EXACT(rr)), "r-", label="vecchio fit (5 esponenziali) sui punti ESATTI")
a.plot(rr, 1e4 * (c6.v(rr) - EXACT(rr)), "g-", label="fit a 6 esponenziali sui punti esatti")
a.plot(pts[:, 0], 0 * pts[:, 0], "k|", ms=12, label="le 18 distanze")
a.set_xlabel("R (bohr)"); a.set_ylabel("fit − curva esatta  (10⁻⁴ hartree)")
a.set_title("B5: il vecchio fit deforma la buca (senza rumore)"); a.legend(fontsize=8)
a.set_xlim(0.85, 3.0); a.axvline(1.4011, color="gray", lw=0.6)

a = ax[0, 1]
res = np.array(r["residuals"]["dE"]); err = np.array(r["residuals"]["err"])
a.errorbar(pts[:, 0], 1e4 * res, 1e4 * err, fmt="ko", ms=4, label="palline (DMC, τ = 0.02) − esatta")
a.axhline(0, color="gray", lw=0.6)
a.set_xlabel("R (bohr)"); a.set_ylabel("10⁻⁴ hartree"); a.set_xscale("log")
a.set_title("B7: punti delle palline − esatto (nessuna curvatura)")
a.legend(fontsize=8)

a = ax[1, 0]
names = list(fx["models"].keys())
we = [fx["models"][k]["dmc"]["we"] for k in names]; sd = [fx["models"][k]["sd_we"] for k in names]
bias = [fx["models"][k]["bias_we"] * 100 for k in names]
y = np.arange(len(names))
a.errorbar(we, y, xerr=sd, fmt="o", color="C0", label="ωe dalle palline (± spread sintetico)")
a.axvline(4401.21, color="k", label="misura 4401.21")
a.axvline(fx["exact_curve"]["we"], color="g", ls="--", label=f"curva BO esatta {fx['exact_curve']['we']:.1f}")
for i, b in enumerate(bias):
    a.text(4480, i, f"distorsione del fit {b:+.2f}%", va="center", fontsize=8)
a.set_yticks(y); a.set_yticklabels(names, fontsize=8); a.set_xlabel("ωe (cm⁻¹)"); a.set_xlim(4250, 4600)
a.set_title("B5/correzione: modelli giudicati SOLO sui dati esatti"); a.legend(fontsize=8, loc="lower left")

a = ax[1, 1]
s = np.array(r["B4"]["syn_we"])
a.hist(s[(s > 4200) & (s < 4500)], bins=50, range=(4200, 4500), color="C1", alpha=0.7, label="vecchio fit su esatto + rumore (2000 copie)")
a.axvline(r["B0"]["python"]["we"], color="r", label="palline, vecchio fit: 4339")
a.axvline(fx["exact_curve"]["we"], color="g", ls="--", label="curva esatta")
a.axvline(4401.21, color="k", label="misura")
a.set_xlabel("ωe (cm⁻¹)"); a.set_title("B4: le palline stanno dove il vecchio fit porta l'esatto")
a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "h2vib_figure.png"), dpi=110)
print("ok")
