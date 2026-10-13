"""nucleo1s_f8_figure.py -> lab/out/nucleo1s_f8_figure.png : error of the alpha^4 shift (two-Z difference) vs the
radius R_f where the zoom chain hands over to the coarse level; learned level vs exact band; post-hoc points."""
import os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
J = json.load(open(os.path.join(OUT, "nucleo1s_f8_judged.json")))
PH = json.load(open(os.path.join(OUT, "nucleo1s_f8_posthoc.json")))
J0s = [3, 4, 5, 6, 8, 12, 16, 24]
R = [(j + 0.5) / 8 for j in J0s]
fig, ax = plt.subplots(figsize=(7.5, 4.6))
for n, c in ((1, "C0"), (2, "C1"), (3, "C2")):
    ax.plot(R, [100 * J[f"exact_J{j}_{n}s"] for j in J0s], "o-", color=c, label=f"{n}s banda esatta (controllo)")
    ax.plot(R, [100 * J[f"learned_J{j}_{n}s"] for j in J0s], "x--", color=c, ms=8, label=f"{n}s livello imparato (Taylor)")
for n, c in ((1, "C0"), (2, "C1")):
    ax.plot([(j + 0.5) / 8 for j in (32, 48, 64)], [100 * PH[f"err_J{j}_{n}s"] for j in (32, 48, 64)], "s:", color=c, mfc="none",
            label=f"{n}s controllo, dopo (non giudicato)")
ax.axhline(0, color="k", lw=0.8)
ax.axhspan(-5, 5, color="0.9", zorder=0)
ax.set_xscale("log")
ax.set_xlabel("R_f = raggio dove la catena dello zoom passa al livello grosso  [a0/Z]")
ax.set_ylabel("errore dello spostamento alfa^4 (differenza Z=1/Z=2)  [%]")
ax.set_title("Fase 8: il livello imparato va come la banda esatta (H2), la colla no (H1, H3, H4)", fontsize=10)
ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "nucleo1s_f8_figure.png"), dpi=130)
print("saved")
