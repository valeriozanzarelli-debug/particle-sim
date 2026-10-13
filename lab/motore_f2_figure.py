"""motore_f2_figure.py - figure of phase 2 (lab/out/motore_f2_figure.png)."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import motore_f2_core as f2
OUT = f2.OUT
L = lambda f: json.load(open(os.path.join(OUT, f)))
A, J, Lo, Po, PH = L("motore_f2_atom.json"), L("motore_f2_judge.json"), L("motore_f2_lot.json"), L("motore_f2_pos.json"), L("motore_f2_posthoc.json")
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
a = ax[0]
labels = ["riga 1s→2p", "P(2p) dopo il calcio", "P(1s) dopo il calcio"]
xs = np.arange(3)
for i, (t, c) in enumerate((("ultra", "tab:red"), ("intel", "tab:blue"))):
    v = [abs(J[t]["line_rel"]), abs(J[t]["P2p_rel"]), abs(J[t]["P1s_rel"])]
    a.bar(xs + (i - 0.5) * 0.35, v, 0.35, color=c, label=f"{t}: cella a0/{16 if t == 'ultra' else 8}")
a.bar([1.35], [abs(PH[1]["P2p_rel"])], 0.2, color="tab:cyan", label="intel, scatola ±18 a0 (dopo)")
a.axhline(0.02, color="tab:red", ls="--", lw=1); a.axhline(0.10, color="tab:blue", ls="--", lw=1)
a.set_yscale("log"); a.set_xticks(xs); a.set_xticklabels(labels)
a.set_ylabel("errore relativo contro l'esatto")
a.set_title("(a) atomo colpito: Ultra e Intel battono alla stessa riga\n(10,202 e 10,156 eV contro 10,204)")
a.legend(fontsize=8)
a = ax[1]
rows = Lo["light"]
a.semilogy([r["N_coarse"] for r in rows], [abs(r["dE_roundtrip"]) for r in rows], "o-", color="tab:green", label="locale a coseni (G3a)")
D1 = L("motore_m3.json")["D1"]
for meth, col in (("cubico", "tab:blue"), ("spettrale", "tab:gray")):
    rr = [r for r in D1 if r["method"] == meth]
    a.semilogy([r["N_coarse"] for r in rr], [abs(r["dE_roundtrip"]) for r in rr], "s--", color=col, label=meth + " (fase 1)")
a.set_xlabel("celle per λ sul livello grosso"); a.set_ylabel("|ΔE/E| andata e ritorno")
a.set_title("(b) cambio di livello LOCALE: energia a 1e-11\nsupporto: blocco ± un blocco vicino (zero esatto oltre)")
a.legend(fontsize=8)
a = ax[2]
p = Po["pure_R"]; q = Po["atomic_packet"]
vals = [q["fake_beat_amp_unprojected"], q["fake_beat_amp_projected"], p["fake_beat_amp_unprojected"]]
a.bar([0, 1, 2], np.maximum(vals, 1e-17), color=["tab:orange", "tab:green", "tab:gray"])
a.set_yscale("log"); a.set_xticks([0, 1, 2])
a.set_xticklabels(["pacchetto atomico\nsenza proiezione", "pacchetto atomico\ncon banda positiva", "R puro fermo\n(segnalato: resta fine)"])
a.set_ylabel("ampiezza del battito finto in P_R")
a.set_title("(c) Intel (B = 192): il battito lento è la Zitterbewegung\nricampionata (ω = %.4f = 2mc² piegato), non le valli" % p["beat_peak_omega"])
fig.suptitle("Motore ottimizzato, fase 2", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "motore_f2_figure.png"), dpi=110); print("ok")
