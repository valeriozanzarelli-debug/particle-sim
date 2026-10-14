"""mol_figure.py - figure of study mol_: lab/out/mol_figure.png (LiH curve + levels, dipole, electron snapshot;
approximated-tier panel when lab/out/mol_approx.json exists)."""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
EV = 27.211386245988; HC = 219474.6313632

res = json.load(open(os.path.join(OUT, "mol_lih.json")))
an = res["analysis"]; C = res["curve"]
Rs = np.array(sorted(v["R"] for v in C.values()))
Es = np.array([C[f"{R:.3f}"]["E"] for R in Rs]); mus = np.array([C[f"{R:.3f}"]["mu_au"] for R in Rs])
Einf = an["E_Li"] + an["E_H"]
approx = os.path.exists(os.path.join(OUT, "mol_approx.json"))
fig, ax = plt.subplots(1, 4 if approx else 3, figsize=(20 if approx else 15, 4.6))
a = ax[0]
cs = CubicSpline(Rs, Es); x = np.linspace(Rs[0], Rs[-1], 400)
a.plot(x, (cs(x) - Einf) * EV, "k-", lw=1.2)
a.plot(Rs, (Es - Einf) * EV, "o", color="C0", label="ECG 4 elettroni (questo studio)")
a.axhline(-2.515, color="C3", ls=":", label="D_e sperimentale 2,515 eV")
a.axvline(3.015, color="C3", ls="--", lw=0.8, label="R_e sperimentale 3,015 a0")
lev = np.array(an["levels_cm"]) / 8065.54429 + (an["E_min"] - Einf) * EV + (an["D_e_eV"] - an["D0_eV"])
for v, e in enumerate(lev[:5]):
    a.hlines(e, 2.6, 3.6, color="C2", lw=0.8)
a.set_xlabel("distanza Li-H (a0)"); a.set_ylabel("energia rispetto a Li + H (eV)")
a.set_title(f"LiH: R_e={an['R_e']:.3f} a0, D_e={an['D_e_eV']:.3f} eV, omega_e={an['omega_e_cm']:.0f} cm-1", fontsize=9)
a.legend(fontsize=7)
a = ax[1]
a.plot(Rs, np.abs(mus) * 2.541746, "o-", color="C1")
a.axhline(5.882, color="C3", ls=":", label="5,882 D (v=0, sperimentale)")
a.set_xlabel("distanza Li-H (a0)"); a.set_ylabel("momento di dipolo (debye)")
a.set_title(f"dipolo: media su v=0 = {an['mu_v0_D']:.3f} D (Li+ H-)", fontsize=9); a.legend(fontsize=7)
a = ax[2]
p = os.path.join(os.environ.get("MOL_DEST", "/mnt/project-files/debiti/dati/molecole"), "LiH/lih_elettroni.npy")
if os.path.exists(p):
    S = np.load(p)[:6000]
    cols = ["C0", "C3", "C0", "C3"]
    for i in range(4):
        a.scatter(S[:, i, 2], S[:, i, 0], s=0.3, alpha=0.25, color=cols[i])
    a.plot([0], [0], "k*", ms=10); a.plot([res["curve"][f"{Rs[np.argmin(abs(Rs - an['R_e']))]:.3f}"]["R"]], [0], "k^", ms=8)
    a.set_xlim(-3, 7); a.set_ylim(-4, 4); a.set_aspect("equal")
    a.set_xlabel("z (a0)  [Li = stella, H = triangolo]"); a.set_ylabel("x (a0)")
    a.set_title("6000 fotografie dei 4 elettroni (blu spin su, rosso giu')", fontsize=9)
if approx:
    A = json.load(open(os.path.join(OUT, "mol_approx.json")))
    a = ax[3]
    rows = A.get("figure_rows", [])
    if rows:
        names = [r["name"] for r in rows]
        a.barh(names, [r["err_rel_bond"] * 100 for r in rows], color="C4")
        a.axvline(0, color="k", lw=0.8)
        a.set_xlabel("errore della lunghezza di legame (%)"); a.set_title("livello approssimato", fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "mol_figure.png"), dpi=130)
print("saved")
