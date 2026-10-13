"""emiss_figure3.py - figure of phase 3 -> lab/out/emiss_figure3.png"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from emiss_core import OUT
from emiss_analysis import load

R = json.load(open(os.path.join(OUT, "emiss_results3.json")))
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
a = ax[0]
dq, ds = load("P5_R1_96"), load("P5_S1sc_96")
a.semilogy(dq["ser_t"], dq["ser_eph"] + dq["ser_abs"] + 1e-40, "C0", label="QE: energia data alla luce dal 2p puro")
a.semilogy(ds["ser_t"], ds["ser_eph"] + ds["ser_abs"] + 1e-40, "C3", label="semiclassico: stessa partenza")
a.set_ylim(1e-35, 1e-2)
a.set_xlabel("t")
a.set_title("S5 (96^3): il 2p puro emette solo con la luce quantistica")
a.legend(fontsize=8)
a = ax[1]
w96, wD = R["R5_numbers"]["omega"], 0.0766903
d5 = json.load(open(os.path.join(OUT, "emiss_d5.json")))
vals = [(w96 / wD)**3, d5["ratio"], (w96 / wD)**3 * d5["ratio"], R["R2v5"]["value"]]
labs = ["omega^3\n(regola/Dirac)", "|d|^2\n(regola/Dirac)", "prodotto", "A misurato / W2\n/ Dirac"]
a.bar(range(4), vals, color=["C1", "C1", "C7", "C2"])
for i, v in enumerate(vals):
    a.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=9)
a.axhline(1, color="k", lw=0.8)
a.set_xticks(range(4))
a.set_xticklabels(labs, fontsize=8)
a.set_ylim(0.8, 1.2)
a.set_title("R2v5 passa per compensazione di due errori")
a = ax[2]
x = json.load(open(os.path.join(OUT, "emiss_x5_16000_bh_full.json")))
pk = x["peaks"]
a.stem([p[0] for p in pk], [p[1] for p in pk])
a.set_yscale("log")
for w, lab in ((x["levels"]["w1s"], "1s"), (x["levels"]["w12"], "2p1/2"), (x["levels"]["w32"], "2p3/2"), (1.06848, "scatola 1,0685")):
    a.annotate(lab, (w, 0.3), fontsize=8, rotation=90)
a.set_xlim(0.85, 1.1)
a.set_xlabel("omega del livello")
a.set_ylabel("peso nello stato dopo l'impulso")
a.set_title("righe deboli = 1,0685 - 2p (elettrone libero nella scatola)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "emiss_figure3.png"), dpi=110)
print("written")
