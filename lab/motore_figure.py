"""motore_figure.py - figure of the study (lab/out/motore_figure.png) from the JSON results."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
L = lambda f: json.load(open(os.path.join(OUT, f)))
A, B1, B2, M2, M3, J, EE = (L("motore_m1_A.json"), L("motore_m1_B1.json"), L("motore_m1_B2.json"), L("motore_m2.json"),
                            L("motore_m3.json"), L("motore_m4_judge.json"), L("motore_m4_EE.json"))
fig, ax = plt.subplots(2, 3, figsize=(17, 10.5))

# (a) error model
a = ax[0, 0]
for c in A["cases"]:
    m = np.linalg.norm(c["meas_dx"]); p = np.linalg.norm(c["spec_dx"])
    a.loglog(p + 1e-4, m + 1e-4, "o", color=plt.cm.viridis(c["ang"] / 45), ms=7)
x = np.logspace(-4, 1, 10)
a.loglog(x, x, "k--", lw=1)
a.set_xlabel("deriva prevista dal modello di Bloch [lunghezze d'onda]")
a.set_ylabel("deriva misurata contro Maxwell continuo [λ]")
a.set_title("(a) A1: il modello d'errore prevede la deriva (16/16)\ncolore = angolo 0° (viola) … 45° (giallo)")

# (b) B1: chosen level vs footprint
a = ax[0, 1]
for tier, col in (("intel", "tab:blue"), ("ultra", "tab:red")):
    rows = B1[tier]["rows"]
    p = [r["footprint"] for r in rows]
    N = [r["N"] for r in rows]
    a.semilogx(p, N, "o-", color=col, label=f"{tier}: celle per λ scelte")
    for r in rows:
        a.annotate("%.2f px" % r["meas_px"] + (" (segn.)" if r["flagged"] else ""), (r["footprint"], r["N"]),
                   textcoords="offset points", xytext=(4, 6 if tier == "intel" else -12), fontsize=7, color=col)
a.set_xlabel("impronta del pixel [λ per pixel]  (sinistra = vicino)")
a.set_ylabel("livello scelto [celle per λ]")
a.set_title("(b) B1: errore sullo schermo ≤ 1 px; costo ÷%.0f (Intel), ÷%.0f (Ultra)"
            % (B1["intel"]["cost_reduction_unflagged"], B1["ultra"]["cost_reduction_unflagged"]))
a.legend(fontsize=8)

# (c) B2
a = ax[0, 2]
for pol, col in (("naive", "tab:orange"), ("passport", "tab:green")):
    h = B2[pol]["history"]
    a.step([x["t_lambda"] for x in h] + [100], [x["N"] for x in h] + [h[-1]["N"]], where="post", color=col,
           label="%s: errore finale %.2f px" % ("ingenua" if pol == "naive" else "passaporto", B2[pol]["final_px"]))
a2 = a.twinx()
t = np.linspace(0, 100, 50)
a2.plot(t, 5 + (0.1 - 5) * t / 100, "k:", lw=1)
a2.set_ylabel("impronta [λ/px] (puntini: la camera si avvicina)")
a.set_yscale("log")
a.set_xlabel("tempo [λ/c]")
a.set_ylabel("celle per λ usate")
a.set_title("(c) B2: l'errore della fisica ha memoria\n(Nanite no). Ingenua 14.7 px, passaporto 0.26 px")
a.legend(fontsize=8, loc="center left")

# (d) cone
a = ax[1, 0]
tab = M2["C1_saving_table"]
a.plot([r["T"] for r in tab], [1 / r["active_fraction"] for r in tab], "s-", color="tab:purple")
a.set_xscale("log", base=2)
a.set_xlabel("orizzonte T [tick] (finestra 48² in un mondo 192²)")
a.set_ylabel("risparmio (celle aggiornate: tutte / cono)")
a.set_title("(d) C1: congelare fuori dal cono di luce è ESATTO\nerrore 0 (al bit), cono 1 cella più stretto: errore %.1f"
            % M2["C1"]["light_err_cone_minus1"])
a.grid(alpha=.3)

# (e) transitions
a = ax[1, 1]
for meth, col in (("spettrale", "tab:green"), ("cubico", "tab:blue"), ("media_blocco", "tab:gray")):
    rows = [r for r in M3["D1"] if r["method"] == meth]
    a.loglog([r["N_coarse"] for r in rows], [abs(r["dE_roundtrip"]) + 1e-16 for r in rows], "o-", color=col,
             label=meth.replace("_", " "))
a.set_xlabel("celle per λ sul livello grosso")
a.set_ylabel("|ΔE/E| nel giro fine→grosso→fine")
a.set_title("(e) D1: cambio di livello della luce\nspettrale = conserva l'energia a 1e-11")
a.legend(fontsize=8)

# (f) dynamics per tier
a = ax[1, 2]
labels = ["velocità 45°", "frequenza 45°", "diffrazione 0°", "cucitura T (30°)", "elettrone v", "elettrone E(1/a0)"]
xs = np.arange(len(labels))
for i, (tier, col) in enumerate((("ultra", "tab:red"), ("intel", "tab:blue"), ("sotto_pavimento", "tab:gray"))):
    j = J[tier]
    e = EE["tiers"][tier]
    vals = [abs(j["EL1_speed"]["45.0"]), abs(j["EL2_freq"]["45.0"]), abs(j["EL3_rows"][0][1] - 1),
            1 - min(r[2] for r in j["EL4_rows"] if r[0] == 30), abs(e["EE1_v_rel"]), abs(e["EE3"]["1.0"]["rel_err_vs_fine"])]
    a.bar(xs + (i - 1) * 0.27, np.maximum(vals, 1e-7), 0.27, color=col,
          label="%s (N=%d, B=%d): luce %s, elettrone %s" % (tier.replace("_", " "), j["N"], e["B"], j["EL_verdict"], j["EE_verdict"]))
a.axhline(0.02, color="tab:red", ls="--", lw=1)
a.axhline(0.10, color="tab:blue", ls="--", lw=1)
a.set_yscale("log")
a.set_xticks(xs); a.set_xticklabels(labels, rotation=25, fontsize=8)
a.set_ylabel("errore relativo contro il riferimento")
a.set_title("(f) La mela in 4K e in 480p: dinamica conservata\ntratteggio = tolleranza Ultra 2 %, Intel 10 %")
a.legend(fontsize=7, loc="upper right", framealpha=0.9); a.set_ylim(1e-5, 2)
fig.suptitle("Motore ottimizzato — lato fisico: errore sullo schermo, cono causale, cambi di livello, dinamica per livello", fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "motore_figure.png"), dpi=110)
print("ok")
