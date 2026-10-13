"""vuoto_f3_figure.py - figure of phase 3: lab/out/vuoto_f3_figure.png"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
X = json.load(open(os.path.join(OUT, "vuoto_f3_exams.json")))
T = X["table"]
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
a = ax[0]
for key, col, lab in (("R|K12", "C3", "regola, taglio K12 (fase 2)"), ("R|K12/3", "C0", "regola, taglio K12/3"),
                      ("C|K12/3", "C2", "ordine ciclico (proposta), K12/3")):
    th = [r["theta"] for r in T[key]["rows"]]
    for j, mk in ((0, "o"), (1, "s"), (2, "^")):
        a.plot(th, [r["dv"][j] for r in T[key]["rows"]], mk + "-", color=col, ms=5, mfc="none" if j == 2 else col,
               label=lab + " " + "xyz"[j] if True else None)
for key, ls in (("cont|K12", "--"), ("cont|K12/3", ":")):
    a.plot([r for r in [0.05, 0.025, 0.0125, 0.00625, 0.003125]], T[key]["dv"], "k" + ls, label="continuo, " + key.split("|")[1])
a.set_xscale("log"); a.set_xlabel("θ = m a  (reticolo più fine ←)"); a.set_ylabel("δv_j / α  (c vestito dell'elettrone − c)")
a.set_title("(a) c vestito per asse (motore esatto)", fontsize=10); a.legend(fontsize=6.5, ncol=2)
a.axhline(0, color="gray", lw=0.5)
b = ax[1]
for key, col in (("R|K12", "C3"), ("R|K12/3", "C0"), ("C|K12/3", "C2")):
    th = [r["theta"] for r in T[key]["rows"]]
    b.plot(th, [r["aniso"] for r in T[key]["rows"]], "o-", color=col, label="anisotropia y − (x+z)/2, " + key)
    b.plot(th, [r["iso"] for r in T[key]["rows"]], "s--", color=col, label="parte isotropa, " + key)
b.plot([0.05, 0.025, 0.0125, 0.00625, 0.003125], T["cont|K12"]["dv"], "k:", label="parte isotropa, continuo")
b.set_xscale("log"); b.set_xlabel("θ = m a"); b.set_ylabel("/ α")
b.set_title("(b) anisotropia e parte isotropa: non calano con a", fontsize=10); b.legend(fontsize=6.5)
b.axhline(0, color="gray", lw=0.5)
c = ax[2]
P = X["g2_pairs"]
ths = [0.05, 0.025, 0.0125, 0.00625]
for p, col in zip(["20", "21", "01", "02", "10", "12"], ["C0", "C1", "C2", "C3", "C4", "C5"]):
    c.plot(ths, [P["slice0|%g" % t][p] for t in ths], "o-", color=col, label="coppia %s (asse Q, asse corrente)" % p)
    c.plot([0.0125], [P["slice1|0.0125"][p]], "x", color=col, ms=9)
c.axhline(1, color="k", lw=1)
c.set_xscale("log"); c.set_xlabel("θ = m a"); c.set_ylabel("F₂ / (α/2π)")
c.set_title("(c) g−2 col vertice H_eff: fetta 0 (o) e fetta 1 (x)", fontsize=10); c.legend(fontsize=6.5)
fig.suptitle("Vuoto che polarizza, fase 3: il c vestito dell'elettrone della regola (Collins–Perez–Sudarsky–Urrutia–Vucetich)", fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "vuoto_f3_figure.png"), dpi=110)
print("saved")
