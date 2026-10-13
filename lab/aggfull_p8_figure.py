"""aggfull_p8_figure.py -> lab/out/aggfull_p8_figure.png"""
import os, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "aggfull_p8_prereg.json")))
E = json.load(open(os.path.join(OUT, "aggfull_p8_err.json")))
J, R = P["config"]["J0_map"], P["config"]["R_glue_a0_over_Z"]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for i, n in enumerate((1, 2)):
    a = ax[i]
    for Z, mk in (("1", "o"), ("3", "s"), ("5", "^")):
        a.plot(R, [100 * E[f"p8_P{Z}_J{j}_{n}s"] for j in J[Z]], mk + "-", label=f"Z = {Z} e {2 * int(Z)}, catena θ/4")
        a.plot([R[-1]], [100 * E[f"th1_P{Z}_J{J[Z][-1]}_{n}s"]], mk, mfc="none", ms=9, label=f"Z = {Z}, catena θ (senza patch)")
    a.axhspan(-1, 1, color="g", alpha=0.1); a.axhline(0, color="gray", lw=0.5); a.axvline(4.06, color="gray", ls=":")
    a.set_xlabel("raggio della colla R (a0/Z)"); a.set_title(f"{n}s: errore α⁴ (coppie alla stessa geometria)")
    a.set_ylim(-13, 4); a.legend(fontsize=6)
ax[0].set_ylabel("errore α⁴ (%)")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "aggfull_p8_figure.png"), dpi=110)
