"""locale_s4_figure.py - figure of phase 2 (S4): precursors, outside-cone energy, hydrogen order, Schwinger."""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
O = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
R = json.load(open(os.path.join(O, "locale_s4_run.json")))
D = json.load(open(os.path.join(O, "locale_s4_diag.json")))
fig, ax = plt.subplots(2, 2, figsize=(12, 9))
a = ax[0, 0]
for s, c in (("frozen", "C0"), ("S4", "C3")):
    v = R[f"P6|{s}|smooth|6"]
    d = np.array(v["diag_line_t120"]); r = np.arange(len(d)) * np.sqrt(2)
    a.semilogy(r, np.maximum(d, 1e-40) / d.max(), c, label=f"{s}: diagonale")
    l = np.array(v["axis_line_t120"]); a.semilogy(np.arange(len(l)), np.maximum(l, 1e-40) / l.max(), c + "--", label=f"{s}: asse")
a.axvline(120 + 6, color="k", lw=1); a.text(127, 1e-2, "cono vero\nr = t + R", fontsize=8)
a.axvline(120 * 5 / 3, color="gray", ls=":"); a.text(201, 1e-2, "cono S4\n1,67c", fontsize=8)
a.set_ylim(1e-18, 2); a.set_xlim(60, 250); a.set_xlabel("r (celle)"); a.set_ylabel("energia / picco")
a.set_title("luce 2-D, sorgente liscia R = 6, t = 120: il 'precursore'"); a.legend(fontsize=8)
a = ax[0, 1]
Rs = [3, 6, 12, 24]
for s, c in (("frozen", "C0"), ("S4", "C3")):
    a.loglog(Rs, [R[f"P6|{s}|smooth|{x}"]["outside"]["120"] for x in Rs], c + "o-", label=f"{s} liscia")
    a.scatter([12.0], [R[f"P6|{s}|disk|12"]["outside"]["120"]], s=90, marker="s", facecolors="none", edgecolors=c, zorder=5,
              label=f"{s} disco netto R=12")
a.axhline(1e-6, color="k", ls=":"); a.text(3, 1.5e-6, "soglia P6", fontsize=8)
a.set_xlabel("raggio della sorgente R (celle)"); a.set_ylabel("frazione di energia fuori dal cono vero, t = 120")
a.set_title("quanta luce arriva 'prima' (fuori dal cono)"); a.legend(fontsize=8, loc="lower left")
a = ax[1, 0]
ths = [0.2, 0.15, 0.1, 0.07, 0.06, 0.05, 0.04, 0.03, 0.025]
for s, c in (("frozen", "C0"), ("S4", "C3")):
    for lab, mk in (("1s", "o-"), ("2p1/2", "^--")):
        a.loglog(ths, [abs(D["D1_hydrogen"][f"{s}|{t}"][lab]) for t in ths], c + mk, label=f"{s} {lab}")
t = np.array(ths); a.loglog(t, 3e-3 * (t / 0.1) ** 2, "k:", lw=1); a.loglog(t, 6.5e-6 * (t / 0.1) ** 4, "k-.", lw=1)
a.text(0.026, 2.5e-4, "~θ²", fontsize=8); a.text(0.026, 4e-9, "~θ⁴", fontsize=8)
a.set_xlabel("θ = m a"); a.set_ylabel("|errore relativo sull'energia di legame|")
a.set_title("idrogeno: S4 è 100x meglio, ma gli stati s non scendono come θ⁴"); a.legend(fontsize=8)
a = ax[1, 1]
for s, c in (("frozen", "C0"), ("S4", "C3")):
    for EEc, mk in ((1, "o-"), (2, "s-"), (4, "^-")):
        th = [0.1, 0.05, 0.025, 0.0125]
        a.semilogx(th, [100 * D["D3_schwinger"][f"{s}|{x}|{EEc}"]["vs_continuum"] for x in th], c + mk, label=f"{s} E/Ec={EEc}")
a.axhspan(-5, 5, color="0.9")
a.set_xlabel("θ = m a"); a.set_ylabel("P reticolo / P continuo (stessa spazzata) - 1  [%]")
a.set_title("Schwinger: S4 ha bisogno di θ ≤ 0,05 (periodo lungo)"); a.legend(fontsize=7, ncol=2)
plt.tight_layout(); plt.savefig(os.path.join(O, "locale_s4_figure.png"), dpi=110)
print("ok")
