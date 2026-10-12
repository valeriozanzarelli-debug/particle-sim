"""box_figure.py - figure of the 'pairs in a box' study -> lab/out/box_figure.png"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
r = json.load(open(os.path.join(OUT, "box_run.json")))
f = json.load(open(os.path.join(OUT, "box_formation.json")))
th = r["th"]
fig, ax = plt.subplots(2, 2, figsize=(12, 9))

a = ax[0, 0]
d = np.array(r["A3_curve"]["d"]); me = np.array(r["A3_curve"]["m_eff"])
a.plot(d * th, me / th, "k.-", label="reticolo (parete di massa), misurata")
dd = np.linspace(2, 80, 300)
a.plot(dd * th, np.sqrt(1 + (np.pi / (dd * th)) ** 2), "r--", label=r"$\sqrt{m^2+(\pi/d)^2}$ (formula del compito)")
a.plot(dd * th, np.sqrt(1 + (np.pi / (2 * dd * th)) ** 2), "b:", label=r"$\sqrt{m^2+(\pi/2d)^2}$")
from scipy.optimize import brentq
pim = [brentq(lambda p: p * x - (np.pi - np.arctan(p / th)), 1e-9, np.pi / x + 1e-9) for x in dd]
a.plot(dd * th, np.sqrt(th**2 + np.array(pim) ** 2) / th, "g-", lw=1, label="Dirac, pareti a massa infinita")
a.axhline(np.sqrt(2), color="gray", lw=0.6); a.axvline(r["A3_dstar_th"], color="gray", lw=0.6)
a.set_xlabel("larghezza d × m  (d in lunghezze Compton)"); a.set_ylabel("m_eff / m")
a.set_ylim(0.9, 6); a.set_xlim(0.2, 8); a.set_xscale("log"); a.legend(fontsize=8)
a.set_title(f"A1b/A3: la massa cresce quando la scatola stringe (d* m = {r['A3_dstar_th']:.2f})")

a = ax[0, 1]
DS = sorted(int(k) for k in r["A1c"])
P = [r["A1c"][str(k)]["P_lat"] for k in DS]; Pp = [r["A1c"][str(k)]["P_pred"] for k in DS]
a.semilogy(np.array(DS) * th, P, "ko-", label="reticolo: coppie per passaggio")
a.semilogy(np.array(DS) * th, Pp, "r--", label=r"$\Sigma_j \exp(-\pi m_j^2/E)$")
a.axhline(r["A1c"][str(DS[0])]["P_1D"], color="b", lw=0.8, label="senza pareti (1+1 D)")
a.set_xlabel("d × m"); a.set_ylabel("coppie (per modo kx che passa)"); a.set_ylim(1e-11, 3)
a.set_title("A1c: campo fisso E = 2 E_crit"); a.legend(fontsize=8)

a = ax[1, 0]
for dval, c in ((8, "C0"), (16, "C1"), (32, "C2"), (128, "C3")):
    p = os.path.join(OUT, f"box_threshold_d{dval}.json")
    if not os.path.exists(p):
        continue
    t = json.load(open(p))
    Om = np.array(t["Omega"]); PP = np.array(t["P"])
    a.plot(Om / th, PP / PP.max(), color=c, label=f"d = {dval}")
    a.axvline(t["two_m_eff"] / th, color=c, ls="--", lw=0.8)
a.axvline(2, color="k", lw=0.8)
a.set_xlabel("frequenza del fotone Ω / m"); a.set_ylabel("assorbimento (normalizzato)")
a.set_title("A2: soglia del fotone (tratteggio = 2 m_eff; nero = 2 m)"); a.legend(fontsize=8)

a = ax[1, 1]
lf = f["ell_f"]
for tag, c, lab in (("free", "k", "senza pareti"), ("strip_d16", "C1", "striscia d = 16 (m -> m_eff)")):
    rows = f["rows"][tag]
    x = np.array([q["ell"] for q in rows]) / lf
    y = np.array([q["ratio"] for q in rows])
    ok = np.array([not q["band_wrap"] for q in rows])
    a.plot(x[ok], y[ok], "o-", color=c, label=lab)
    a.plot(x[~ok], y[~ok], "x", color=c)
a.axhline(1, color="gray", lw=0.6); a.axvline(1, color="gray", lw=0.6)
a.set_xscale("log"); a.set_xlabel("lunghezza del campo ℓ / ℓ_f   (ℓ_f = 2m/E)")
a.set_ylabel("coppie / stima locale (LCFA)")
a.set_title("A4: lunghezza di formazione (x = gradino oltre la banda)"); a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "box_figure.png"), dpi=110)
print("ok")
