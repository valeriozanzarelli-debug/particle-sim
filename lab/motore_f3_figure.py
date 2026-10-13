"""motore_f3_figure.py - figure of phase 3 -> lab/out/motore_f3_figure.png"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
L = lambda p: json.load(open(os.path.join(OUT, "motore_f3_%s.json" % p)))
lot, sw, U, T = L("lot"), L("sweep"), L("ultra"), L("truth")
A2 = json.load(open(os.path.join(OUT, "motore_f2_atom.json")))
ex2p = 0.073728
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))

a = ax[0]
lab = ["%d %s" % (r["L"], r["bell"]) for r in lot["calibration"]]
a.bar(range(len(lab)), [max(abs(r["dEkin_roundtrip"]), 1e-16) for r in lot["calibration"]], color="C0", label="calibrazione k0=0.5/a0")
for i, r in enumerate(lot["electron"]):
    a.scatter([len(lab) - 1 + 0.15 * (i + 1)], [abs(r["dEkin_roundtrip"])], color="C3", zorder=3,
              label="esame k0=%g/a0" % r["k0_a0"] if True else None)
a.axhline(1e-6, color="k", ls="--", lw=1, label="soglia 1e-6")
a.axhline(7.0e-5, color="C1", ls=":", lw=1, label="fase 2 (32 sin2, k0=1)")
a.set_yscale("log"); a.set_xticks(range(len(lab))); a.set_xticklabels(lab, rotation=30)
a.set_ylabel("|errore energia cinetica| nel giro")
a.set_title("(a) campana al bordo del blocco: si sceglie la prima che passa")
a.legend(fontsize=7)

a = ax[1]
rows = sorted(sw["rows"].values(), key=lambda r: r["box_a0"])
b = [r["box_a0"] for r in rows]
a.plot(b, [100 * (r["P2p"] / ex2p - 1) for r in rows], "o-", label="Intel a0/8")
a.plot([12], [100 * (A2["ultra"]["kick"]["P2p"] / ex2p - 1)], "s", color="C3", label="Ultra a0/16, scatola 12 (fase 2)")
a.plot([U["box_a0"]], [100 * (U["P2p"] / ex2p - 1)], "*", ms=14, color="C3", label="Ultra, scatola scelta da sola")
a.axhspan(-2, 2, color="C2", alpha=0.12, label="tolleranza Ultra 2 %")
a.axvline(sw["chosen_box_a0"], color="k", ls=":", lw=1)
a.set_xlabel("mezzo lato della scatola (a0)"); a.set_ylabel("P(2p0) dopo il calcio, scarto %")
a.set_title("(b) la scatola nel passaporto: scelta al livello economico")
a.legend(fontsize=7)

a = ax[2]
order = [s for s in ("k0.05", "k0.2", "k0.5", "k1.0", "S1", "S2") if s in T]
x = np.arange(len(order))
a.bar(x - 0.2, [max(T[s]["w_out"], 1e-16) for s in order], 0.4, label="peso fuori base w_out")
a.bar(x + 0.2, [max(T[s]["rho"], 1e-16) for s in order], 0.4, label="rho = w_out/(1-P1s)")
a.axhline(0.01, color="C0", ls="--", lw=1, label="soglia w_out = tol^2")
a.axhline(0.10, color="C1", ls="--", lw=1, label="soglia rho = tol")
for i, s in enumerate(order):
    a.text(i, 3, ("AMMESSA" if T[s]["admitted"] else "FINE") + "\nv_c/v_g=%.2f" % T[s]["v_ratio"], ha="center", fontsize=7)
a.set_yscale("log"); a.set_ylim(1e-16, 30)
a.set_xticks(x); a.set_xticklabels(order)
a.set_title("(c) cartolina: il criterio rimanda al fine gli atomi colpiti")
a.legend(fontsize=7, loc="lower left")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "motore_f3_figure.png"), dpi=110)
print("ok")
