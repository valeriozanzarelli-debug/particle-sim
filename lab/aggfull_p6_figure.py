"""aggfull_p6_figure.py - phase 6 figure -> lab/out/aggfull_p6_figure.png"""
import os, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
G = json.load(open(os.path.join(OUT, "aggfull_p6_glue_err.json")))
J = np.array([3, 4, 5, 6, 8, 12, 16, 24, 32, 48, 64])
x = (J - 0.5) / 8
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for i, n in enumerate((1, 2)):
    b = np.array([G[f"base_J{j}_{n}s"] for j in J]); l = np.array([G[f"lc_J{j}_{n}s"] for j in J])
    p = 4 * np.exp(-2 * x) if n == 1 else 1.6 * (1 - x / 2) ** 2 * np.exp(-x)
    a = ax[i]
    a.plot(x, 100 * b, "o-", label="T senza accoppiamento (nucleo1s fase 8)")
    a.plot(x, 100 * l, "s-", label="T + accoppiamento imparato (lingua comp. grande)")
    a.plot(x, 100 * (b + p), "k:", label="T + Darwin 'spalmato' esterno (formula)")
    a.axhline(0, color="gray", lw=0.5); a.axhspan(-5, 5, color="g", alpha=0.08)
    a.set_xscale("log"); a.set_xticks([0.3, 0.5, 1, 2, 4, 8]); a.set_xticklabels(["0.3", "0.5", "1", "2", "4", "8"]); a.minorticks_off(); a.set_ylim(-30, 150 if n == 1 else 80)
    a.set_xlabel("raggio della colla R (a0/Z)"); a.set_ylabel(f"errore alfa^4 {n}s (%)")
    a.set_title(f"{n}s: differenza Z=1 / Z=2"); a.legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "aggfull_p6_figure.png"), dpi=110)
