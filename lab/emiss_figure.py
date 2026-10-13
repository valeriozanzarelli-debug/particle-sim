"""emiss_figure.py - figure of the study "Emissione vera" -> lab/out/emiss_figure.png (panels skip missing data)"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from emiss_core import OUT

R = json.load(open(os.path.join(OUT, "emiss_results.json")))
fig, ax = plt.subplots(2, 3, figsize=(17, 10))

# (a) E1 spectra
a = ax[0, 0]
for key, lab, col in (("E1_spectrum_pi", "luce pi (E_z)", "C0"), ("E1_spectrum_sigma", "luce sigma- (E_y)", "C3")):
    if key in R:
        s = R[key]
        S = np.array(s["S"])
        a.plot(s["om"], S / S.max(), col, label=lab)
if "E1_spectrum_pi" in R:
    a.axvline(R["E1_spectrum_pi"]["w32"], color="k", ls="--", lw=0.8, label="1s->2p3/2 (reticolo)")
    a.axvline(R["E1_spectrum_pi"]["w12"], color="gray", ls=":", lw=1.2, label="1s->2p1/2 (reticolo)")
a.set_xlim(0.070, 0.086)
a.set_xlabel("omega (unita' m c^2/hbar)")
a.set_title("E1: righe emesse dopo l'impulso (48^3)")
a.legend(fontsize=8)

# (b) decay: quantum vs semiclassical
a = ax[0, 1]
if "R1_48_numbers" in R:
    import numpy as np
    from emiss_analysis import load
    d = load("P_R1_48")
    n = R["R1_48_numbers"]
    a.plot(d["ser_t"], d["ser_norm"], "C0", label="quantistico (QE), 2p puro, K=1")
    a.plot(d["ser_t"], np.exp(-n["A_rule"] * n["W2"] * d["ser_t"]), "k--", lw=0.8, label="exp(-A_regola W2 t), nessun fit")
    try:
        ds = load("P_S1sc")
        pe = np.array(ds["ser_pexc"])[:, 0]
        a.plot(ds["ser_t"], pe, "C3", label="semiclassico, 2p puro: proiezione che batte (luce emessa 2e-33)")
    except Exception:
        pass
    a.set_ylim(0.93, 1.002)
a.set_xlabel("t")
a.set_ylabel("probabilita' nello stato eccitato")
a.set_title("S1: il 2p puro decade da solo solo con la luce quantistica")
a.legend(fontsize=8)

# (c) K10 / SC10
a = ax[0, 2]
if "K10_series" in R:
    k = R["K10_series"]
    a.plot(k["t"], k["norm"], "C0", label="QE K=10 (2p puro)")
    a.plot(k["t"], np.exp(-k["A"] * np.array(k["t"])), "k--", lw=0.8, label="esponenziale")
PH = json.load(open(os.path.join(OUT, "emiss_posthoc.json")))
if "SC10_series" in R:
    s = R["SC10_series"]
    a.plot(s["t"], s["p"], "C3", label="semiclassico K=10 da (1s+2p)/sqrt2")
    a.plot(s["t"], s["logistic"], "C3", ls=":", label="logistica con A fittato su K10 (esame: FAIL)")
    tt_ = np.array(s["t"]); A0 = PH["SC10_rule"]["A0"]
    a.plot(tt_, 0.5 * np.exp(-A0 * tt_) / (0.5 + 0.5 * np.exp(-A0 * tt_)), "k", ls="-.", lw=0.8, label="logistica con A della regola a t=0 (dopo, dichiarata)")
a.set_xlabel("t")
a.set_title("esponenziale (quantistico) vs logistica (classico)")
a.legend(fontsize=8)

# (d) toy TWA
a = ax[1, 0]
f = os.path.join(OUT, "emiss_twa.json")
if os.path.exists(f):
    tw = json.load(open(f))
    t = np.array(tw["t"])
    a.plot(t, tw["QE"], "k", lw=2, label="QE (Wigner-Weisskopf)")
    a.plot(t, tw["SC_p1.0"], "C3", label="semiclassico da e: fermo")
    a.plot(t, tw["HTWA_e"], "C1", label="Wigner ibrido da e")
    a.plot(t, tw["HTWA_e_rate_eq"], "C1", ls=":", lw=0.8)
    a.plot(t, tw["HTWA_g"], "C2", label="Wigner ibrido da g: pompato!")
    a.plot(t, tw["HTWA_g_rate_eq"], "C2", ls=":", lw=0.8)
    a.plot(t, tw["STWA_e"], "C4", label="Wigner a spin da e")
    a.set_xlabel("tick (giocattolo 1-D)")
    a.set_ylabel("p eccitato")
a.set_title("le opzioni per l'emissione spontanea (giocattolo)")
a.legend(fontsize=7)

# (e) Laboratory beat
a = ax[1, 1]
lb = json.load(open(os.path.join(OUT, "emiss_lab1s.json")))
csv = "/mnt/project-files/ripresa/dati/atomo_accoppiato/serie.csv"
if os.path.exists(csv):
    dd = np.genfromtxt(csv, delimiter=",", names=True)
    a.plot(dd["t"], dd["p_exc"], "C0", lw=0.8, label="Laboratorio dopo la cura dell'1s (export 01:46, con luce)")
fr = lb["free"]
a.plot(fr["t"], fr["p_run"][:len(fr["t"])], "C3", lw=0.8, label="1s del Laboratorio lasciato solo, SENZA luce")
a.set_xlim(0, 1600)
a.set_xlabel("t")
a.set_ylabel("p_exc")
a.set_title(f"il battito e' del 1s sporco: periodo previsto {fr['period_pred']:.0f}, misurato {lb['lab_series']['beat_period']:.0f}")
a.legend(fontsize=8)

# (f) rates and donut
a = ax[1, 2]
lab, v1, v2 = [], [], []
for key, nm in (("R1_48_numbers", "Nd 48^3"), ("R1_96_numbers", "Nd 96^3"), ("P_N1_numbers", "Nd radiale"),
                ("P_H1a_numbers", "H a=0.5"), ("P_H1b_numbers", "H a=0.25")):
    if key in R:
        q = R[key]
        lab.append(nm)
        v1.append(q["A_norm"] / (q["A_rule"] * q["W2"]))
x = np.arange(len(lab))
a.bar(x, v1, 0.6, color="C0", label="A misurato / A previsto dalla regola")
lab.append("Nd 48^3, dipolo\ndel momento (dopo)")
a.bar(len(lab) - 1, np.mean([r_["ratio"] for r_ in PH["P_R1_48_windows"]]), 0.6, color="C7")
for key, nm in (("N1", "N1: Nd vs Dirac"), ("H1", "H1: H vs Dirac")):
    if key in R:
        lab.append(nm)
        a.bar(len(lab) - 1, R[key]["value"], 0.6, color="C2")
a.axhline(1, color="k", lw=0.8)
a.axhspan(0.97, 1.03, color="gray", alpha=0.2)
a.set_xticks(range(len(lab)))
a.set_xticklabels(lab, rotation=30, fontsize=8)
a.set_ylim(0.6, 1.4)
a.set_title("ritmi con alpha vero: blu regola, verde Dirac (banda +-3 %)")
f = os.path.join(OUT, "emiss_ideas.json")
if os.path.exists(f):
    idd = json.load(open(f))
    ins = a.inset_axes([0.68, 0.68, 0.31, 0.31], projection="polar")
    th = np.arccos(np.array(idd["polar_cos"]))
    pr = np.array(idd["polar"]) / idd["fit_sin2"]
    ins.plot(th, pr, "C1.")
    ins.plot(-th, pr, "C1.")
    tt = np.linspace(0, 2 * np.pi, 200)
    ins.plot(tt, np.sin(tt)**2, "k", lw=0.6)
    ins.set_theta_zero_location("N")
    ins.set_xticks([])
    ins.set_yticks([])
    ins.set_title("ciambella sin^2", fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "emiss_figure.png"), dpi=110)
print("written emiss_figure.png")
