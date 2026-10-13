"""Figura della fase 2 dello studio idrogeno_bohr: convergenza dell'1s ibrido nella grandezza del blocco fine
(regola vecchia e regola T) e parità degli stati. Legge lab/out/idrogeno_bohr_risultati.json.
Uso: python3 lab/idro_bohr_fase2_figura.py"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
R = json.load(open("lab/out/idrogeno_bohr_risultati.json"))["fase2"]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for lab, k, c in [("regola vecchia (72 passi)", "blocco_regola_vecchia_72passi", "#888"), ("regola T (48 passi)", "blocco_regola_T_48passi", "#2a7")]:
    x, y = zip(*R[k]); ax[0].plot(x, y, "o-", color=c, label=lab)
ax[0].axhline(R["grosso_solo_pct"]["T"], color="#2a7", ls=":", lw=1); ax[0].text(3, R["grosso_solo_pct"]["T"] + 0.002, "solo grosso (T)", fontsize=8, color="#2a7")
ax[0].axhline(0.037, color="r", ls="--", lw=1); ax[0].text(3, 0.039, "massa del cammino di base (+0,037 %)", fontsize=8, color="r")
ax[0].axhline(0, color="k", lw=0.6)
ax[0].set_xlabel("lato del blocco fine (celle grosse a0/16), celle fini a0/256"); ax[0].set_ylabel("1s: errore contro Dirac × μ (%)")
ax[0].set_title("ibrido con reflux: il blocco cresce, l'1s si sposta", fontsize=10); ax[0].legend(fontsize=8)
st = ["1s", "2pz", "2px", "2py", "2s"]
v = [R["R4a_T"]["regola_vecchia_parita_1meno"][s] for s in st]
ax[1].bar(np.arange(5) - 0.2, v, 0.4, color="#888", label="regola vecchia")
ax[1].bar(np.arange(5) + 0.2, [max(R["R4a_T"]["parita_1meno_max"], 1e-16)] * 5, 0.4, color="#2a7", label="regola T (max)")
ax[1].set_yscale("log"); ax[1].set_xticks(range(5)); ax[1].set_xticklabels(st)
ax[1].set_ylabel("1 − |parità| (mescolamento s–p)"); ax[1].set_title("parità degli stati (a0/8, ±18 a0)", fontsize=10); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig("lab/out/idrogeno_bohr_fase2_figure.png", dpi=110); print("scritto lab/out/idrogeno_bohr_fase2_figure.png")
