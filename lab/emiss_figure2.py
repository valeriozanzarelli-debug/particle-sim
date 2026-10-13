"""emiss_figure2.py - figure of phase 2 -> lab/out/emiss_figure2.png"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from emiss_core import OUT

R = json.load(open(os.path.join(OUT, "emiss_results2.json")))
R1 = json.load(open(os.path.join(OUT, "emiss_results.json")))
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
a = ax[0]
for key, lab, col in (("E4_spectrum_pi", "luce pi (E_z)", "C0"), ("E4_spectrum_sigma", "luce sigma- (E_y)", "C3")):
    s = R[key]
    S = np.array(s["S"])
    a.plot(s["om"], S / S.max(), col, label=lab)
a.axvline(R["E4_spectrum_pi"]["w32"], color="k", ls="--", lw=0.8, label="1s->2p3/2")
a.axvline(R["E4_spectrum_pi"]["w12"], color="gray", ls=":", lw=1.2, label="1s->2p1/2")
a.set_xlim(0.070, 0.086)
a.set_xlabel("omega (m c^2/hbar)")
a.set_title("fase 2: con nucleo imparato pi -> 2p3/2, sigma- -> 2p1/2")
a.legend(fontsize=8)

a = ax[1]
q = R["R4_numbers"]
t = np.array(q["t"])
a.plot(t, q["norm"], "C0", label="QE fase 2 (nucleo imparato), K=1")
a.plot(t, np.exp(-q["A_rule"] * q["W2"] * t), "k--", lw=0.8, label="exp(-A_regola W2 t), nessun fit")
s = R["S4_series"]
a.plot(s["t"], s["p"], "C3", label="semiclassico: proiezione sul 2p che batte (luce 4e-33)")
d1 = R1.get("K10_series")
a.set_xlabel("t")
a.set_ylabel("probabilita' eccitata")
a.set_title(f"R4 = {R['R4']['value']:.3f} (fase 1: R1_48 = 0,674)")
a.legend(fontsize=8)

a = ax[2]
d4 = R["D4"]["value"]
v3 = json.load(open(os.path.join(OUT, "emiss_posthoc.json")))["dipoles_48"]
names = list(d4.keys())
x = np.arange(len(names))
a.bar(x - 0.2, [v3[k]["ratio"] for k in names], 0.4, color="C7", label="fase 1 (nucleo nudo, 1s a 4 valli)")
a.bar(x + 0.2, [d4[k]["ratio"] for k in names], 0.4, color="C2", label="fase 2 (nucleo imparato)")
a.axhline(1, color="k", lw=0.8)
a.set_xticks(x)
a.set_xticklabels(names, fontsize=8)
a.set_ylabel("sum |<1s|r|2p>|^2 / Dirac continuo")
a.set_title("dipoli 1s->2p sul 48^3")
a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "emiss_figure2.png"), dpi=110)
print("written emiss_figure2.png")
