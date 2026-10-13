"""dilat_figure.py - figura dello studio dilat_ (lab/out/dilat_figura.png)"""
import json, os, sys, shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dilat_core import rate_lattice

HERE = os.path.dirname(os.path.abspath(__file__))
R = lambda n: json.load(open(os.path.join(HERE, "out", "dilat_r_%s.json" % n)))
PRE = json.load(open(os.path.join(HERE, "out", "dilat_prereg.json")))["esami"]

fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))

# (a) i due pacchetti (densita') che si separano, mentre l'orologio tiene il ritmo
d = R("E2d")
rho = np.array(d["rho_sup"])
t = np.array(d["t"][1:])
x = np.arange(rho.shape[1]) * d["dx_rho"]
sel = (x > 100) & (x < 280)
a0 = ax[0, 0]
a0.imshow(rho[:, sel], aspect="auto", origin="lower", extent=[x[sel][0], x[sel][-1], t[0], t[-1]], cmap="magma")
a0.set_xlabel("x  [ħ/mc]")
a0.set_ylabel("t  [ħ/mc²]")
a0.set_title("(a) regola, densità: sovrapposizione p̄=0,5 e 0,8 (Δ=0,3)\nle frange cambiano e i centri si allontanano (%.0f ħ/mc), il ritmo no (%.0e)"
             % (d["separazione_finale"], d["variazione_rel"]), fontsize=10)

# (b) dilatazione classica: errore della regola vs a
e1 = R("E1")
b = ax[0, 1]
A = [1 / 4, 1 / 8, 1 / 16, 1 / 32]
for r in e1["righe"]:
    if r["pbar"] == 0:
        continue
    err = [abs(r["err"][str(a)]) for a in A]
    b.loglog(A, err, "o-", label="p̄ = %.1f mc  (τ/t → %.4f)" % (r["pbar"], r["R_cont"]))
    b.loglog(A, [abs(r["coeff_prev"]) * a * a for a in A], "k:", lw=0.8)
b.set_xlabel("passo reticolare a  [ħ/mc]")
b.set_ylabel("|τ/t (regola) − ⟨m/E⟩| / ⟨m/E⟩")
b.set_title("(b) dilatazione classica: errore = previsione O(a²) (punteggiato)", fontsize=10)
b.legend(fontsize=8)

# (c) correzione quantistica nel regime relativistico
e2 = R("E2")
c = ax[0, 2]
for al, ph, col in [(np.pi / 8, 0.0, "C0"), (np.pi / 8, np.pi, "C1"), (np.pi / 4, np.pi, "C2")]:
    rows = [r for r in e2["E2b"]["righe"] if abs(r["alpha"] - al) < 1e-9 and abs(r["phi"] - ph) < 1e-9]
    s = [r["s"] for r in rows]
    c.plot(s, [r["dq"][str(1 / 32)] for r in rows], "o", color=col, ms=7,
           label="regola a=1/32, θ=%s, φ=%s" % ("π/8" if al < 0.5 else "π/4", "0" if ph == 0 else "π"))
    c.plot(s, [r["dq_cont"] for r in rows], "-", color=col, lw=1)
    c.plot(s, [r["dq_SA_NR"] for r in rows], "--", color=col, lw=1)
c.axhline(0, color="k", lw=0.5)
c.set_xlabel("(p̄₂ − p̄₁)/Δ   (p̄₁ = 0,5 mc, Δ = 0,3 mc)")
c.set_ylabel("Δ_Q = τ/t (sovrapp.) − τ/t (miscuglio)")
c.set_title("(c) correzione quantistica: — continuo esatto, -- Smith–Ahmadi (NR)", fontsize=10)
c.legend(fontsize=8)

# (d) visibilita'
e4 = R("E4")
dd = ax[1, 0]
cv = e4["curva"]
dd.plot(cv["t"], cv["V_sup"], label="sovrapposizione (regola)")
dd.plot(cv["t"], cv["V_mix"], label="miscuglio (regola)")
tv = PRE["E4"]["config"]["t"]
dd.plot(tv, PRE["E4"]["previsione"]["V_sup"], "k+", ms=9, label="continuo (previsione)")
dd.plot(tv, PRE["E4"]["previsione"]["V_mix"], "kx", ms=7)
dd.set_xlabel("t  [ħ/mc²]")
dd.set_ylabel("visibilità |⟨ψ₋|ψ₊⟩| (orologio δ = 0,05 mc²)")
dd.set_title("(d) visibilità dell'orologio: p̄ = 0,5 e 1,0, Δ = 0,3", fontsize=10)
dd.legend(fontsize=8)

# (e) Zitterbewegung dell'orologio (chiralita' R pura)
e6 = R("E6")
ee = ax[1, 1]
tz = np.array(PRE["E6"]["config"]["t"])
prev = np.array(PRE["E6"]["previsione"]["tau_cont"])
rm = PRE["E6"]["previsione"]["ritmo_medio"]
ee.plot(tz, prev - rm * tz, "k-", lw=2, label="continuo di Dirac")
for a in [1 / 8, 1 / 32]:
    ee.plot(tz, np.array(e6["tau"][str(a)]) - rm * tz, "--", label="regola a = 1/%d" % round(1 / a))
ee.set_xlabel("t  [ħ/mc²]")
ee.set_ylabel("τ(t) − ⟨pm/E²⟩ t")
ee.set_title("(e) stato di sola chiralità R: l'orologio trema (ZB)\nritmo medio %.4f = ⟨v m/E⟩ (il ramo − va all'indietro)" % rm,
             fontsize=10)
ee.legend(fontsize=8)

# (f) ritmo su tutta la zona di Brillouin: cos(ka) * sqrt(1 - v^2), doppione all'indietro
ff = ax[1, 2]
for a in [1 / 4, 1 / 8]:
    k = np.linspace(-np.pi / a, np.pi / a, 2001)
    ff.plot(k * a / np.pi, rate_lattice(k, a, 1.0), label="regola a = 1/%d" % round(1 / a))
k = np.linspace(-np.pi / (1 / 8), np.pi / (1 / 8), 2001)
ff.plot(k * (1 / 8) / np.pi, 1 / np.sqrt(1 + k ** 2), "k:", label="continuo m/E")
ff.axhline(0, color="k", lw=0.5)
ff.set_xlabel("k a / π")
ff.set_ylabel("ritmo dell'orologio dτ/dt = √(1−v²)·cos(ka)")
ff.set_title("(f) ritmo esatto della regola: vicino al doppione (ka = ±π)\nl'orologio gira all'indietro", fontsize=10)
ff.legend(fontsize=8)

plt.tight_layout()
p = os.path.join(HERE, "out", "dilat_figura.png")
plt.savefig(p, dpi=110)
os.makedirs("/mnt/project-files/problemi-aperti", exist_ok=True)
shutil.copy(p, "/mnt/project-files/problemi-aperti/dilat_figura.png")
print("scritto", p)
