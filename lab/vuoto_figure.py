"""vuoto_figure.py - figure of the study 'Vuoto che polarizza' -> lab/out/vuoto_figure.png"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from vuoto_core import chi1d_cont, chi_1d, species_weight, pitR_uehling, uehling_V, LOGCOEF, ALPHA, PI
from vuoto_levels import load_pitR, make_interp, load_3d

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
D = load_3d()
fig, ax = plt.subplots(3, 2, figsize=(12, 13.5))

# (a) 1-D
a = ax[0, 0]
qm = np.exp(np.linspace(np.log(0.1), np.log(20), 40))
a.plot(qm, [chi1d_cont(q) for q in qm], "k-", lw=2, label="Dirac continuo 1+1 D (esame)")
w = species_weight()
for th, c in ((0.1, "C1"), (0.025, "C0")):
    a.plot(qm, [chi_1d(q * th, th, w, n=2000, beta=4.0).real for q in qm], "o", ms=3, color=c, label=f"zig-zag, specie fisica, θ={th}")
a.axhline(0, color="C3", ls="--", label="zig-zag, mare intero (esatto 0)")
a.set_xscale("log"); a.set_xlabel("q / m"); a.set_ylabel("χ(q)  (risposta del vuoto 1+1 D)")
a.set_title("(a) 1+1 D: il doppione annulla tutto; la specie fisica = Dirac"); a.legend(fontsize=8)

# (b) 3-D pit(0) vs theta
b = ax[0, 1]
ths = [0.1, 0.05, 0.025, 0.0125]
p0 = []
for th in ths:
    try:
        p0.append(np.mean(list(load_pitR(th)[2].values())))
    except Exception:
        p0.append(np.nan)
p0 = np.array(p0)
x = -np.log(ths)
b.plot(x, p0, "o", color="C0", ms=8, label="zig-zag (media x,y,z), q→0")
ok = np.isfinite(p0)
if ok.sum() >= 2:
    b.plot(x, p0[ok][-1] + LOGCOEF * (x - x[ok][-1]), "k--", label="pendenza 1/(6π²) ↔ 2α/3π")
b.set_xlabel("−ln θ  (θ = m a: reticolo più fine →)"); b.set_ylabel("Π(0)/(4πα)")
b.set_title("(b) rinormalizzazione della carica: il logaritmo"); b.legend(fontsize=8)

# (c) Pi_R vs Uehling, whole sea
c = ax[1, 0]
qq = np.exp(np.linspace(np.log(0.02), np.log(50), 80))
c.plot(qq, [pitR_uehling(q) for q in qq], "k-", lw=2, label="Uehling (esame)")
for th, col in ((0.1, "C1"), (0.0125, "C0")):
    try:
        qs, tab, _ = load_pitR(th)
        c.plot(qs, tab, "o", ms=3, color=col, label=f"zig-zag specie fisica θ={th}")
    except Exception:
        pass
fq, fv = [], []
for k, v in D.items():
    if k.startswith("full|"):
        fq.append(float(k.split("|")[3])); fv.append(v[0])
if fq:
    o = np.argsort(fq); fq = np.array(fq)[o]; fv = np.array(fv)[o]
    c.plot(fq, -fv, "s--", color="C3", ms=4, label="mare intero (tutti i doppioni): −Π(q)/(4πα) > 0, anti-schermo, cresce come ln(1/q)")
c.set_xscale("log"); c.set_yscale("log")
c.set_xlabel("q / m"); c.set_ylabel("Π_R(q)/(4πα)")
c.set_title("(c) polarizzazione rinormalizzata contro Uehling"); c.legend(fontsize=8)

# (d) r-space
d = ax[1, 1]
try:
    L = json.load(open(os.path.join(OUT, "vuoto_levels.json")))
    rr = np.array(L["V3_r"])
    d.plot(rr, np.array(L["V3_ratio_lat_ref"]), "o-", label="V_U zig-zag / V_U Uehling (stesso nucleo smussato)")
    d.axhspan(0.98, 1.02, color="0.9", label="soglia ±2 %")
    H = L["H"]
    F = json.load(open(os.path.join(OUT, "vuoto_posthoc_Hfine.json")))
    e2 = F["2s_0.00625"] + (F["2s_0.00625"] - F["2s_0.0125"]) / 3
    e1 = F["1s_0.00625"] + (F["1s_0.00625"] - F["1s_0.0125"]) / 3
    txt = (f"H 2S: {H['2s']['extrap_MHz']:.2f} MHz pre-reg. (FAIL); a posteriori {e2:.2f}  (rif. {H['H_ref_MHz']['2s']:.2f})\n"
           f"H 1S: {H['1s']['extrap_MHz']:.1f} MHz pre-reg. (FAIL); a posteriori {e1:.1f}  (rif. {H['H_ref_MHz']['1s']:.1f})\n"
           f"Nd⁵⁹⁺ 1s: {L['Nd']['extrap2_eV']:.3f} eV  (rif. {L['Nd']['ref_eV']:.3f})")
    d.text(0.03, 0.04, txt, transform=d.transAxes, fontsize=8, family="monospace")
    d.set_ylim(0.975, 1.025)
except Exception as e:
    d.text(0.1, 0.5, "livelli non ancora calcolati: " + str(e)[:60])
d.set_xscale("log"); d.set_xlabel("r / λ_C"); d.set_ylabel("rapporto")
d.set_title("(d) potenziale di Uehling dal mare e spostamento dei livelli"); d.legend(fontsize=8, loc="upper right")
# (e) phase 2: self-energy loop of the rule - mass shift log and g - 2
e = ax[2, 0]
try:
    SE = json.load(open(os.path.join(OUT, "vuoto_f2_se.json")))
    for lat, mk, lab in ((False, "--", "continuo (controllo formule)"), (True, "o-", "regola: elettrone-walk × luce-walk")):
        ths = [t for t in (0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125) if "dm|%s|%g|0|std" % ("lat" if lat else "cont", t) in SE]
        x = -np.log(ths)
        for comp, col in (("dmC_over_m", "C0"), ("dmT_over_m", "C1"), ("dm_over_m", "C3")):
            y = [SE["dm|%s|%g|0|std" % ("lat" if lat else "cont", t)][comp] / ALPHA for t in ths]
            e.plot(x, y, mk, color=col, ms=4, label=(lab + " " + {"dmC_over_m": "Coulomb", "dmT_over_m": "trasverso", "dm_over_m": "totale"}[comp]) if True else None)
    e.set_xlabel("−ln θ  (= ln(1/(m a)))"); e.set_ylabel("δm/(α m)")
    e.set_title("(e) fase 2: autoenergia dalla regola, pendenze 1/π, 1/2π, 3/2π"); e.legend(fontsize=6.5, ncol=2, loc="lower right")
    X = json.load(open(os.path.join(OUT, "vuoto_f2_exams.json")))
    if "G1" in X:
        g = X["G1"]; g0 = X.get("G0", {})
        txt = "g−2 (α/2π = 1): continuo %.4f, regola %.3f (θ→0)\n" % (g0.get("F2_extrap_over_alpha2pi", np.nan), g["F2_extrap_over_alpha2pi"])
        txt += "regola a θ = " + ", ".join("%s: %.3f" % (k, v["F2_lam0"]) for k, v in g["table"].items())
        e.text(0.02, 0.97, txt, transform=e.transAxes, fontsize=7.5, va="top", family="monospace")
except Exception as ex:
    e.text(0.1, 0.5, "fase 2 non ancora calcolata: " + str(ex)[:60])

# (f) phase 2: Bethe log of the radial rule vs non-relativistic QM at the same photon cutoff; H with the projector
f = ax[2, 1]
try:
    X = json.load(open(os.path.join(OUT, "vuoto_f2_exams.json")))
    KRY = [30, 100, 300, 1000, 3000]
    for n, col in ((1, "C0"), (2, "C1")):
        T = X["L_table"]["%d|0.05" % n]
        f.plot(KRY, [T[str(k)]["L_NR"] for k in KRY], "-", color=col, lw=2, label="%dS QM non relativistica (esame)" % n)
        f.plot(KRY, [T[str(k)]["L_walk"] for k in KRY], "o", color=col, mfc="none", label="%dS regola θ=0,05 (pre-registrata)" % n)
        if "L_table_A2" in X and "%d|0.05" % n in X["L_table_A2"]:
            T2 = X["L_table_A2"]["%d|0.05" % n]
            f.plot(KRY, [T2[str(k)]["L_walk"] for k in KRY], "x", color=col, ms=8, label="%dS regola θ=0,05 (emend. A-L2)" % n)
        f.axhline({1: 2.984128556, 2: 2.811769893}[n], color=col, ls=":", lw=1)
    f.set_ylim(0.5, 4.2); f.set_xscale("log"); f.set_xlabel("taglio sull'energia del fotone K / Ry"); f.set_ylabel("ln k₀(K)")
    H = X.get("H_table")
    if H:
        f.text(0.35, 0.05, "H col proiettore anche sull'elettrone (Uehling):\n1S %.2f MHz (rif. %.2f), 2S %.3f (rif. %.3f)\n1S–2S %+.2f MHz (rif. %+.2f)" % (
            H["1s"]["extrap"], H["1s"]["ref"], H["2s"]["extrap"], H["2s"]["ref"], X["H4"]["value_MHz"], X["H4"]["ref_MHz"]),
            transform=f.transAxes, fontsize=7.5, family="monospace")
    f.set_title("(f) fase 2: logaritmo di Bethe (punteggiate: valori dei libri, K→∞)"); f.legend(fontsize=7, loc="upper left")
except Exception as ex:
    f.text(0.1, 0.5, "fase 2 non ancora calcolata: " + str(ex)[:60])
fig.suptitle("Vuoto che polarizza: il mare di Dirac della regola zig-zag in un potenziale di Coulomb", fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "vuoto_figure.png"), dpi=120)
print("saved")
