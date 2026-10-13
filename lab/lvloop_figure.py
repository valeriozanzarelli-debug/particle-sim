"""lvloop_figure.py -> lab/out/lvloop_figura.png
(a) maximum lattice step a_max (Planck units) vs scaling exponent p of a loop-induced Lorentz violation
    delta = alpha * kappa * theta^p (kappa = 1), for each bound, mass-type (solid) and momentum-type (dashed);
(b) symmetry of the interacting rule (lvloop_sym.json): spectral distance for pairs of impurities, and group sizes;
(c) the dressed-electron limiting speed of vuoto_ phase 3 (exact engine, origin/master 35afcc5) vs theta.
Usage: python lvloop_figure.py <vuoto_f2_exams.json> <vuoto_f3_exams.json>  (both from origin/master via git show)"""
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
f3path = sys.argv[2] if len(sys.argv) > 2 else None
if f3path and os.path.exists(f3path):
    t = json.load(open(f3path))["table"]
    th = [r["theta"] for r in t["R|K12"]["rows"]]
    for key, col, nm in (("R|K12", "C3", "regola, K12"), ("R|K12/3", "C1", "regola, K12/3"), ("C|K12/3", "C2", "ordine a giro, K12/3")):
        rows = t[key]["rows"]
        c.plot(th, [r["iso"] for r in rows], "-o", color=col, label=nm + ": isotropa")
        c.plot(th, [r["aniso"] for r in rows], "--s", color=col, ms=4, label=nm + ": anisotropia y")
    c.plot(th, t["cont|K12"]["dv"], "k:", label="continuo, stesso taglio (isotropa)")
    c.axhline(0, color="gray", lw=0.6)
    c.set_xscale("log"); c.set_xlabel("theta = m a  (5 passi, massa fissa)"); c.set_ylabel("delta v / alpha (contro la luce nuda)")
    c.legend(fontsize=6.5)
    c.set_title("(c) vuoto_ fase 3 (master 35afcc5): c dell'elettrone vestito\nnessuna parte cala con a (p = 0)", fontsize=10)
elif f2path and os.path.exists(f2path):
    f2 = json.load(open(f2path))
    th = sorted(float(x) for x in f2["G1"]["table"])
    dv = np.array([f2["G1"]["table"][str(x)]["dv_over_alpha"] for x in th])
    for j, (nm, col) in enumerate((("x", "C0"), ("y", "C3"), ("z", "C2"))):
        c.plot(th, dv[:, j], "-o", color=col, label="fase 2, asse " + nm)
    c.set_xscale("log"); c.legend(fontsize=7)
plt.tight_layout()
out = os.path.join(here, "out", "lvloop_figura.png")
plt.savefig(out, dpi=130)
print(out)
