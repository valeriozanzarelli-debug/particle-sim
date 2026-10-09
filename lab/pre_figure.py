"""pre_figure.py - one figure for the pre-experiments P1..P5 -> lab/out/pre_figure.png"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pre_p4 import U1

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
J = lambda f: json.load(open(os.path.join(OUT, f)))
p1 = J("pre_p1p5.json"); ph = J("pre_p1_posthoc.json"); p2 = J("pre_p2.json"); p3 = J("pre_p3.json")
p3b = J("pre_p3b.json"); p4 = J("pre_p4.json")
z1 = np.load(os.path.join(OUT, "pre_p1p5.npz")); z3 = np.load(os.path.join(OUT, "pre_p3.npz"))
z3b = np.load(os.path.join(OUT, "pre_p3b.npz"))

fig, ax = plt.subplots(2, 4, figsize=(20, 9.5))
# (a) P1
a = ax[0, 0]; t = z1["t"]
for k, lab in (("P1_uniform_M16_lamw", "M=16, lam=w"), ("P1_uniform_M1_lamw", "M=1 (controllo), lam=w"),
               ("P1_uniform_M1_lam0.0", "M=1, lam=0 (Bell minimo)"), ("P1_half_M1_lam0.0", "M=1, meta' scatola, lam=0"),
               ("P1a_lamw", "partenza |psi|^2")):
    a.semilogy(t, np.maximum(z1[k], 1e-5), label=lab)
a.axhline(p1["P1a_lamw"]["F"], color="k", ls=":", label="pavimento F")
a.set_title("P1 reticolo 24x24: H(t) - anche M=1 si mescola"); a.set_xlabel("tick"); a.legend(fontsize=7)
# (b) P1e refinement
a = ax[0, 1]; R = ph["P1e_refinement"]
th = [0.5, 0.25, 0.125]
for M in (1, 16):
    for lam, mk in (("w", "o-"), ("0.0", "s--")):
        a.loglog(th, [R[f"theta{x}_M{M}_lam{lam}"]["tau_phys"] for x in th], mk, label=f"M={M}, lam={lam}")
a.set_xlabel("theta = m a (reticolo piu' fine ->)"); a.set_ylabel("tau in unita' 1/m"); a.invert_xaxis()
a.set_title("P1e (post-hoc): raffinando, il rumore del salto\nsi spegne solo con lam=0"); a.legend(fontsize=7)
# (c) P2
a = ax[0, 2]; keys = ["shuffle", "noise", "random_graph", "random_state"]
rates = p2["rates"]; x = np.arange(4)
a.bar(x - 0.27, [rates[k]["fp_holes"] for k in keys], 0.27, label="buchi (b1+b2>=1)")
a.bar(x, [rates[k]["fp_lowdim"] for k in keys], 0.27, label="'dimensione' 0.8-3.5")
a.bar(x + 0.27, [np.mean([r["betti"][0] > 1 for r in p2["nulls"][k]]) for k in keys], 0.27, label="b0>1 (grappoli)")
a.axhline(0.1, color="r", ls=":", label="soglia 10%")
a.set_xticks(x); a.set_xticklabels(["mescolata", "rumore", "grafo a caso", "stato a caso"], fontsize=8)
a.set_title("P2 falsi positivi della ricostruzione (input SENZA forma)"); a.legend(fontsize=7)
# (d) P2 robust persistence
a = ax[0, 3]
for i, k in enumerate(keys):
    v = [r["p1_rel"] for r in p2["nulls"][k]]
    a.plot(np.full(len(v), i) + np.random.default_rng(i).uniform(-.15, .15, len(v)), v, "o", ms=4, alpha=.6)
a.plot([4], [p2["positive_control_torus12"]["p1_rel"]], "r*", ms=15, label="toro vero (controllo)")
a.axhline(0.4, color="k", ls=":", label="soglia proposta 0.4 (post-hoc)")
a.set_xticks(range(5)); a.set_xticklabels(["mescolata", "rumore", "grafo", "stato", "toro"], fontsize=8)
a.set_ylabel("buco piu' persistente / diametro"); a.set_title("P2 (post-hoc): un criterio che separa"); a.legend(fontsize=7)
# (e) P3 H(t)
a = ax[1, 0]; t3 = z3["t"]
a.semilogy(t3, z3["plain_uniform_H"], label="senza rinforzo")
a.semilogy(t3, z3["reinf_uniform_H"], label="rinforzo (eta=.05, mu=1)")
a.semilogy(t3, z3["posthoc_strong_uniform_H"], label="forte (eta=.3, mu=2)")
a.semilogy(z3b["t"], z3b["low_strong_uniform_H"], label="forte, budget 0.3")
a.semilogy(t3, z3["reinf_psi2_H"], label="rinforzo, partenza |psi|^2", alpha=.6)
a.axhline(p3["reinf_uniform"]["F"], color="k", ls=":")
a.set_xlim(0, 600); a.set_xlabel("tick"); a.set_title("P3 rinforzo: nessun blocco dell'equilibrio"); a.legend(fontsize=7)
# (f) P3b highways map
a = ax[1, 1]; W = z3b["low_strong_uniform_W3000"]; n = W.shape[1]
for j in (0, 1):
    for xx in range(n):
        for yy in range(n):
            w = W[j][xx, yy]
            if w <= 0:
                continue
            x2, y2 = (xx + 1, yy) if j == 0 else (xx, yy + 1)
            if x2 >= n or y2 >= n:
                continue
            a.plot([xx, x2], [yy, y2], color=plt.cm.viridis(w), lw=0.3 + 4 * w)
a.set_aspect("equal"); a.set_title(f"P3b strade (budget 0.3): Gini {p3b['low_strong_uniform']['gini_end']:.2f}, "
                                   f"top10% {p3b['low_strong_uniform']['top10_end']:.2f}")
a.set_xticks([]); a.set_yticks([])
# (g) P4 dispersion
a = ax[1, 2]
s = np.linspace(0, 1, 300)
path = np.r_[np.c_[np.pi * s, 0 * s], np.c_[np.pi + 0 * s, np.pi * s], np.c_[np.pi * (1 - s), np.pi * (1 - s)]]
for r, c in ((0.0, "C0"), (0.5, "C3")):
    om = np.array([np.sort(-np.angle(np.linalg.eigvals(U1(kx, ky, 0.5, r)))) for kx, ky in path])
    a.plot(om, ".", ms=1.5, color=c)
    a.plot([], [], color=c, label=f"Wilson r={r}")
a.set_xticks([0, 300, 600, 900]); a.set_xticklabels(["(0,0)", "(pi,0)", "(pi,pi)", "(0,0)"])
a.set_ylabel("quasi-energia omega"); a.legend(fontsize=7)
a.set_title("P4 i due fogli: (pi,pi) = copia esatta;\nWilson la spinge via (1 pezzo solo da r>=1)")
# (h) P5
a = ax[1, 3]; P5 = p1["P5"]
a.plot([16, 24, 32], [P5[f"n{k}_M16"]["tau"] for k in (16, 24, 32)], "o-", label="tau vs lato n (M=16)")
a2 = a.twiny()
a2.plot([1, 4, 16, 36], [P5[f"n24_M{m}"]["tau"] for m in (1, 4, 16, 36)], "s--", color="C1", label="tau vs M (n=24)")
a.set_xlabel("lato n"); a2.set_xlabel("modi M"); a.set_ylabel("tau (tick)"); a.set_ylim(0, 30)
a.legend(loc="lower left", fontsize=7); a2.legend(loc="lower right", fontsize=7)
a.set_title("P5 tau ~ 10-16 tick: non cresce con n ne' con M\n(celle 4x4 fisse: mescolamento locale da rumore)", fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "pre_figure.png"), dpi=110)
print("saved")
