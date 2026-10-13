"""phase 4 figure -> lab/out/vuoto_f4_figure.png"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
E = json.load(open(os.path.join(D, "vuoto_f4_exams.json")))
TH = [0.05, 0.025, 0.0125, 0.00625]
PAIRS = ["20", "21", "01", "02", "10", "12"]
cols = {"R K12": "C0", "R K12/3": "C1", "C K12/3": "C2"}
fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
a = ax[0, 0]
for k, c in cols.items():
    G = E["G"][k]
    a.plot(TH, [G["per_theta"][str(t)]["mean"] for t in TH], "o-", color=c, label=k + " (media)")
    for p in PAIRS:
        a.plot(TH, [G["per_theta"][str(t)]["pairs_kmin0"][p] for t in TH], ".", color=c, alpha=0.4)
a.axhline(1, color="k", lw=0.8)
a.set_xscale("log"); a.set_xlabel("theta = m a"); a.set_ylabel("a / (alpha/2pi)")
a.set_title("(a) g-2, vertice esatto (6 orientazioni)"); a.legend(fontsize=8)
a = ax[0, 1]
for k, c in cols.items():
    G = E["G"][k]
    a.plot(TH, [max(G["per_theta"][str(t)]["spread"], 1e-12) for t in TH], "o-", color=c, label=k)
a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("theta"); a.set_ylabel("dispersione (alpha/2pi)")
a.set_title("(b) dispersione fra orientazioni"); a.legend(fontsize=8)
a = ax[1, 0]
for k, c in cols.items():
    K = E["K"][k]; L = E["L"][k]
    ev = np.array(K["electron_dv_over_alpha"]).mean(1)
    lv = [np.mean([L["per_theta"][str(t)][x]["mean"] for x in "xyz"]) for t in TH]
    a.plot(TH, ev, "o-", color=c, label=k + " elettrone")
    a.plot(TH, lv, "s--", color=c, label=k + " luce")
a.axhline(0, color="k", lw=0.8)
a.set_xscale("log"); a.set_xlabel("theta"); a.set_ylabel("delta c / (alpha c)  (isotropo)")
a.set_title("(c) c vestito: elettrone contro luce"); a.legend(fontsize=7)
a = ax[1, 1]
L = E["L"]["R K12"]["per_theta"]
for i, x in enumerate("xyzd"):
    for j in range(2):
        a.plot(TH, [L[str(t)][x]["dc"][j] for t in TH], "o-" if j == 0 else "x--", color="C%d" % i,
               label="luce lungo %s, pol %d" % (x, j))
ev = np.array(E["K"]["R K12"]["electron_dv_over_alpha"])
for i, x in enumerate("xyz"):
    a.plot(TH, ev[:, i], ":", color="C%d" % i, label="elettrone lungo %s" % x)
a.set_xscale("log"); a.set_xlabel("theta"); a.set_ylabel("delta c / (alpha c)")
a.set_title("(d) regola congelata, taglio K12: per asse e polarizzazione"); a.legend(fontsize=6, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(D, "vuoto_f4_figure.png"), dpi=110)
print("ok")
