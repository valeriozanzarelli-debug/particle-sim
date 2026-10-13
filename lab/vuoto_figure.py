"""vuoto_figure.py - figure of the study 'Vuoto che polarizza' -> lab/out/vuoto_figure.png"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from vuoto_core import chi1d_cont, chi_1d, species_weight, pitR_uehling, uehling_V, LOGCOEF, ALPHA, PI
from vuoto_levels import load_pitR, make_interp

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
D = json.load(open(os.path.join(OUT, "vuoto_3d.json")))
fig, ax = plt.subplots(2, 2, figsize=(12, 9))

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
    c.plot(fq, -(fv - fv[-1]), "s--", color="C3", ms=4, label="mare intero: −(Π(q)−Π(3m)) (log a grande distanza)")
c.set_xscale("log"); c.set_yscale("symlog", linthresh=1e-5)
c.set_xlabel("q / m"); c.set_ylabel("Π_R(q)/(4πα)")
c.set_title("(c) polarizzazione rinormalizzata contro Uehling"); c.legend(fontsize=8)

# (d) r-space
d = ax[1, 1]
try:
    L = json.load(open(os.path.join(OUT, "vuoto_levels.json")))
    rr = np.array(L["V3_r"])
    d.plot(rr, np.array(L["V3_ratio_lat_ref"]), "o-", label="V_U zig-zag / V_U Uehling (stesso nucleo smussato)")
    d.axhspan(0.98, 1.02, color="0.9", label="soglia ±2 %")
    txt = (f"H 2S: {L['H']['2s']['extrap_MHz']:.2f} MHz (rif. {L['H_ref_MHz']['2s'] if 'H_ref_MHz' in L else L['H']['H_ref_MHz']['2s']:.2f})\n"
           if False else "")
    H = L["H"]
    txt = (f"H 2S: {H['2s']['extrap_MHz']:.2f} MHz  (rif. {H['H_ref_MHz']['2s']:.2f})\n"
           f"H 1S: {H['1s']['extrap_MHz']:.1f} MHz  (rif. {H['H_ref_MHz']['1s']:.1f})\n"
           f"Nd⁵⁹⁺ 1s: {L['Nd']['extrap2_eV']:.3f} eV  (rif. {L['Nd']['ref_eV']:.3f})")
    d.text(0.05, 0.05, txt, transform=d.transAxes, fontsize=10, family="monospace")
except Exception as e:
    d.text(0.1, 0.5, "livelli non ancora calcolati: " + str(e)[:60])
d.set_xscale("log"); d.set_xlabel("r / λ_C"); d.set_ylabel("rapporto")
d.set_title("(d) potenziale di Uehling dal mare e spostamento dei livelli"); d.legend(fontsize=8, loc="upper right")
fig.suptitle("Vuoto che polarizza: il mare di Dirac della regola zig-zag in un potenziale di Coulomb", fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "vuoto_figure.png"), dpi=120)
print("saved")
