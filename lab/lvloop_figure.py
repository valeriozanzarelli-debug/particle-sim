"""lvloop_figure.py -> lab/out/lvloop_figura.png
(a) maximum lattice step a_max (Planck units) vs scaling exponent p of a loop-induced Lorentz violation
    delta = alpha * kappa * theta^p (kappa = 1), for each bound, mass-type (solid) and momentum-type (dashed);
(b) symmetry of the interacting rule (lvloop_sym.json): spectral distance for pairs of impurities, and group sizes;
(c) the dressed-electron speed counterterm of vuoto_ phase 2 (read from origin/master, slice 0) vs theta.
Usage: python lvloop_figure.py [path to vuoto_f2_exams.json from origin/master]"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lvloop_prereg as P

here = os.path.dirname(os.path.abspath(__file__))
sym = json.load(open(os.path.join(here, "out", "lvloop_sym.json")))
f2path = sys.argv[1] if len(sys.argv) > 1 else None

fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))
a = ax[0]
tests = [(P.BOUNDS["iso_photon_faster"][0], "gap isotropo, fotone + veloce (sciami, 6e-21)", "C0"),
         (P.BOUNDS["iso_electron_faster"][0], "gap isotropo, elettrone + veloce (Granchio, 2,5e-20)", "C1"),
         (P.BOUNDS["iso_two_sided"][0], "gap isotropo, LEP sincrotrone (5e-15)", "C2"),
         (P.BOUNDS["aniso_electron"][1], "anisotropia elettrone, Yb+ (1e-21)", "C3"),
         (P.BOUNDS["aniso_photon_nonbiref"][0], "anisotropia luce, risonatori (1e-18)", "C4"),
         (P.BOUNDS["photon_biref"][0], "birifrangenza luce (1e-32)", "C5")]
ps = np.array([0.5, 1, 1.5, 2, 2.5, 3])
for b, lab, col in tests:
    am = np.array([P.a_max_mass_type(b["B"], p) for p in ps]) / P.L_PL
    a.semilogy(ps, am, "-o", color=col, ms=3, label=lab)
    if b.get("E_eV"):
        amp = np.array([P.a_max_momentum_type(b["B"], p, b["E_eV"]) for p in ps]) / P.L_PL
        a.semilogy(ps, amp, "--", color=col, lw=1)
a.axhline(1, color="k", lw=1.2); a.text(2.55, 1.6, "Planck", fontsize=8)
a.axhline(4.8e-31 / P.L_PL, color="gray", ls=":"); a.text(0.55, 4.8e-31 / P.L_PL * 1.8, "giro 1 (luce libera)", fontsize=8, color="gray")
a.axhline(0.1 * P.LAMBDA_C / P.L_PL, color="gray", ls="-."); a.text(0.55, 0.1 * P.LAMBDA_C / P.L_PL * 1.8, "Laboratorio", fontsize=8, color="gray")
a.annotate("p = 0: escluso per ogni a\n(serve kappa < B/alpha)", xy=(0.5, 1e-30), fontsize=8, color="C3")
a.set_xlabel("esponente p  (violazione = alpha (m a)^p ; tratteggio: alpha (E a)^p)")
a.set_ylabel("a_max / l_P")
a.set_ylim(1e-32, 1e40); a.set_xlim(0.4, 3.1)
a.legend(fontsize=6.5, loc="upper left"); a.set_title("(a) passo massimo contro esponente di scaling", fontsize=10)

b = ax[1]
pr = sym["pair_rule_n8"]
labs = list(pr.keys()); vals = [max(pr[k], 1e-16) for k in labs]
b.bar(range(len(labs)), vals, color=["C2" if v < 1e-9 else "C3" for v in vals])
b.set_yscale("log"); b.set_xticks(range(len(labs))); b.set_xticklabels(labs, fontsize=8)
b.set_ylabel("distanza fra gli spettri (Hausdorff)")
g = {k: v["n_sym"] for k, v in sym.items() if isinstance(v, dict) and "n_sym" in v}
txt = "simmetrie con V casuale:\nregola: %d (D4h, asse y)\nvariante xyz fisso: %d\nvariante 6 permutazioni: %d (Th)" % (
    g.get("rule_n6_seed11"), g.get("fixed_n6_seed11"), g.get("six_n6_seed11"))
b.text(0.03, 0.97, txt, transform=b.transAxes, va="top", fontsize=8)
b.set_title("(b) coppie di impurezze, regola, n = 8:\nx = z esatto, y diverso", fontsize=10)

c = ax[2]
if f2path and os.path.exists(f2path):
    f2 = json.load(open(f2path))
    th = sorted(float(t) for t in f2["G1"]["table"])
    dv = np.array([f2["G1"]["table"][str(t)]["dv_over_alpha"] for t in th])
    dv0 = np.array([f2["G0"]["table"][str(t)]["dv_over_alpha"] for t in th])
    for j, (nm, col) in enumerate((("x", "C0"), ("y", "C3"), ("z", "C2"))):
        c.plot(th, dv[:, j], "-o", color=col, label="regola, fetta 0, asse " + nm)
    c.plot(th, dv0[:, 0], "k--", label="continuo, stesso proiettore")
    c.plot(th, dv[:, 1] - 0.5 * (dv[:, 0] + dv[:, 2]), "m-s", label="anisotropia y - (x+z)/2")
    c.set_xscale("log"); c.set_xlabel("theta = m a"); c.set_ylabel("delta v / alpha")
    c.legend(fontsize=7)
c.set_title("(c) c dell'elettrone vestito (vuoto_ fase 2):\nnon cala con a", fontsize=10)
plt.tight_layout()
out = os.path.join(here, "out", "lvloop_figura.png")
plt.savefig(out, dpi=130)
print(out)
