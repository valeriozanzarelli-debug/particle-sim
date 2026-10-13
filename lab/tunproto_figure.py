"""tunproto_figure.py - figura del protocollo: previsioni con bande d'errore -> lab/out/tunproto_figura.png"""
import json, os, sys, shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tunproto_core import OUT, spectral_rule, dirac_spectral
from tunproto_prereg import c_tilt, c_plane

PRE = json.load(open(os.path.join(OUT, "tunproto_prereg.json")))
RC = json.load(open(os.path.join(OUT, "tunproto_rC_pulito.json")))
RN = json.load(open(os.path.join(OUT, "tunproto_rC_rumore.json")))
RW = json.load(open(os.path.join(OUT, "tunproto_rD.json")))
C, D = PRE["C_fibra"], PRE["D_guide"]
th, d = C["theta"], C["delta"]
N1 = C["discriminazioni"]["D1"]["N_5sigma_per_punto_per_orologio"]
N2 = C["discriminazioni"]["D2"]["N_5sigma"]

fig, ax = plt.subplots(2, 3, figsize=(17, 9))
SC = json.load(open(os.path.join(OUT, "tunproto_rC_scan.json")))
# (a) fibra: inclinazione contro numero di siti, V0 = E
a = ax[0, 0]
Ns = np.array(C["N_siti"])
for c, col, lab in (("m", "C3", "orologio di massa  τ_mz"), ("V", "C0", "orologio di potenziale  τ_z (Büttiker–Landauer)")):
    y = np.array([RC["d0.15_V2.0_N%d" % N][c]["tz"] for N in Ns]) / th
    e = np.array([np.sqrt(c_tilt(RC["d0.15_V2.0_N%d" % N][c]["sz"], d) / N1) for N in Ns]) / th
    sysn = np.array([abs(RC["d0.15_V2.0_N%d" % N][c]["tz"] - C["val"]["V2.0_N%d" % N]["spett"][c]["tz"]) for N in Ns]) / th
    tot = np.hypot(e, sysn)
    a.fill_between(Ns, y - tot, y + tot, color=col, alpha=0.25, lw=0)
    a.errorbar(Ns, y, e, fmt="o", color=col, capsize=3, label=lab)
xx = np.linspace(0, 9, 10)
a.plot(xx, xx / np.cos(th) + (RC['d0.15_V2.0_N8']['m']['tz'] / th - 8 / np.cos(th)), "--", color="C3", lw=1, label="regola, barriera spessa: pendenza 1/cos θ = %.3f" % (1 / np.cos(th)))
a.plot(xx, xx, ":", color="k", lw=1, label="L/c (1 passo per sito)")
a.set_xlabel("spessore della barriera N (siti = intervalli di tempo)")
a.set_ylabel("tempo letto dall'ancella (passi)")
a.set_title("(a) fibra, θ = 0,5, V0 = E: inclinazione (barre: N = %.0f trasmessi; banda: + taglio del pacchetto)" % N1, fontsize=9)
a.legend(fontsize=7)
a.set_xlim(1, 9)
# (b) fibra: nel piano contro V0/E, N = 8
a = ax[0, 1]
V0s = np.linspace(1.5, 2.5, 21)
for c, col, lab in (("m", "C3", "massa  τ_my"), ("V", "C0", "potenziale  τ_y")):
    ycur = np.array([spectral_rule(th, 2.0, v, 8 * th, 2.0, d, c, nk=201)["ty"] for v in V0s]) / th
    a.plot(V0s / 2.0, ycur, "-", color=col, lw=1, label=lab + " (previsione, pacchetto σ = 4 siti)")
    pts = [RC["d0.15_V%.1f_N8" % v][c] for v in C["V0"]]
    y = np.array([p["ty"] for p in pts]) / th
    e = np.array([np.sqrt(c_plane(p["sx"], p["sy"], d) / N2) for p in pts]) / th
    a.errorbar(np.array(C["V0"]) / 2.0, y, e, fmt="o", color=col, capsize=3)
a.axhline(0, color="k", lw=0.5)
a.axvline(1, color="k", lw=0.5, ls=":")
a.set_xlabel("V0 / E")
a.set_ylabel("angolo nel piano / δ  (passi)")
a.set_title("(b) fibra, N = 8 siti: l'orologio di massa va all'indietro per V0 > E (barre: N = %.0f)" % N2, fontsize=9)
a.legend(fontsize=7)
# (c) guide d'onda: inclinazione contro L (Compton), V0 = E
a = ax[1, 0]
Lc = np.array(D["L_guide"]) / 6.0
for c, col, lab in (("m", "C3", "massa"), ("V", "C0", "potenziale")):
    y = np.array([RW["pulito"]["V3.0_L%d" % Lg][c]["tz_m"] for Lg in D["L_guide"]])
    pr = np.array([D["previsioni_dirac"]["V3.0_L%d" % Lg][c]["tz"] for Lg in D["L_guide"]])
    a.plot(Lc, pr, "s--", mfc="none", color=col, label=lab + ": Dirac continuo")
    sysb = 0.0
    if c == "m":
        di = RW["disordine"]["dispersione"]
        co = abs(RW["contaminazione"]["tz_m_on"] - RW["contaminazione"]["tz_m_pulito"])
        sysb = np.hypot(di, co)
    a.errorbar(Lc, y, sysb if c == "m" else None, fmt="o", color=col, capsize=3,
               label=lab + ": guide (modi accoppiati)" + (" ± disordine ⊕ η=1 %" if c == "m" else ""))
a.plot([0, 5], [0, 5], ":", color="k", lw=1, label="L/c")
a.set_xlabel("spessore L (lunghezze Compton = 6 guide)")
a.set_ylabel("tempo letto (1/σ = %.2f cm)" % (1 / D["sigma_massa_cm"]))
a.set_title("(c) guide d'onda: κ = 4,5 cm⁻¹, σ = 1,5 cm⁻¹, z = 12 cm, polarizzazione = ancella", fontsize=9)
a.legend(fontsize=7)
a.set_xlim(1.5, 4.5)
# (d) significativita' contro numero di trasmessi
a = ax[1, 1]
Nn = np.logspace(1, 6, 200)
for N5, lab, col in ((N1, "D1: pendenza massa ≠ pendenza potenziale (V0 = E)", "C3"),
                     (C["discriminazioni"]["D1b"]["N_5sigma_per_orologio"], "D1b: un solo spessore N = 8", "C1"),
                     (N2, "D2: τ_my < 0 a V0 = 1,2 E", "C2")):
    a.loglog(Nn, 5 * np.sqrt(Nn / N5), color=col, label=lab + " (5σ a N = %.0f)" % N5)
a.axhline(5, color="k", ls=":", lw=1)
a.set_xlabel("particelle trasmesse rivelate (per punto, per orologio)")
a.set_ylabel("significatività (σ)")
a.set_title("(d) quante ripetizioni: solo rumore di proiezione dell'ancella", fontsize=9)
a.legend(fontsize=7)
# (e) rumore tollerabile (post-hoc)
a = ax[0, 2]
clean = np.polyfit(np.array([4, 6, 8]) * th, [RC["d0.15_V2.0_N%d" % N]["m"]["tz"] for N in (4, 6, 8)], 1)[0]
ss = np.array([0.001, 0.002, 0.005, 0.01, 0.02, 0.05])
for kind, col, lab in (("fase_comune", "C0", "fase, uguale nei due rami"), ("fase_indipendente", "C1", "fase, diversa nei due rami"),
                       ("moneta_comune", "C2", "angolo della moneta")):
    a.semilogx(ss, [SC["%s_%g" % (kind, x)]["pendenza_m"] / clean for x in ss], "o-", color=col, label=lab)
a.axhspan(0.9, 1.1, color="k", alpha=0.08)
a.axhline(0, color="k", lw=0.5)
a.set_xlabel("rumore bianco per sito e per passo (rad, rms)")
a.set_ylabel("pendenza di τ_mz / valore pulito")
a.set_title("(e) post-hoc: rumore tollerabile in fibra (banda grigia: ±10 %)", fontsize=9)
a.legend(fontsize=7)
# (f) perche': il rumore 'scalda' e fa passare sopra la barriera
a = ax[1, 2]
for kind, col in (("fase_comune", "C0"), ("moneta_comune", "C2")):
    a.loglog(ss, [SC["%s_%g" % (kind, x)]["PT_N8"] for x in ss], "o-", color=col, label=kind.replace("_", " "))
a.axhline(RC["d0.15_V2.0_N8"]["m"]["PT"], color="k", ls=":", label="senza rumore (tunnel vero)")
a.set_xlabel("rumore (rad, rms)")
a.set_ylabel("probabilità trasmessa, N = 8, V0 = E")
a.set_title("(f) il rumore che cambia nel tempo non conserva l'energia:\nuna parte passa SOPRA la barriera e diluisce l'orologio", fontsize=9)
a.legend(fontsize=7)
fig.tight_layout()
p = os.path.join(OUT, "tunproto_figura.png")
fig.savefig(p, dpi=130)
os.makedirs("/mnt/project-files/problemi-aperti", exist_ok=True)
shutil.copy(p, "/mnt/project-files/problemi-aperti/tunproto_figura.png")
print("ok", p)
