"""fw1vuoto_figure.py - figure of docs/FW1_VUOTO.md -> lab/out/fw1vuoto_figura.png (and a copy in
/mnt/project-files/problemi-aperti/ if that folder exists). Run after fw1vuoto_analyze.py."""
import os
import sys
import json
import shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fw1vuoto_analyze as A
from fw1vuoto_core import pitR_uehling, LOGCOEF as L

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "fw1vuoto_figura.png")
TH = A.TH
fig, ax = plt.subplots(2, 2, figsize=(12, 9))

# (a) frozen whole sea vs FW1 whole sea: pit(q) per pi-cell, theta = 0.025
a = ax[0, 0]
qs = [0.05, 0.1, 0.2, 0.4, 1, 3, 10]
for th, ls in ((0.05, "--"), (0.025, "-")):
    y = [A.pit("S", "F0", "cellS", th, "x", q) / L for q in qs]
    a.semilogx(qs, y, "k" + ls, marker="o", ms=4, label=f"congelata, mare intero, θ={th}")
for th, ls in ((0.05, "--"), (0.025, "-")):
    q7 = [0.025, 0.05, 0.2, 0.3, 1, 3, 10]
    y = [A.pit("S", "FW1", "cellG", th, "x", q) / L for q in q7]
    a.semilogx(q7, y, "C0" + ls, marker="s", ms=4, label=f"FW1, mare intero, θ={th}")
a.set_xlabel("q / m")
a.set_ylabel("Π(q) per cella  [unità 1/(6π²)]")
a.set_title("(a) mare intero: congelata (coni a π) contro FW1")
a.legend(fontsize=8)
a.grid(alpha=0.3)

# (b) Uehling shape at theta = 0.0125
b = ax[0, 1]
qe = [0.2, 0.3, 1, 3, 10]
for lab, args, st in (("congelata + proiettore", ("S", "F0", "ball"), "k-o"),
                      ("FW1 + proiettore", ("S", "FW1", "ball"), "C1-^"),
                      ("FW1 mare intero (densità a fine periodo)", ("S", "FW1", "cellG"), "C0-s"),
                      ("FW1 mare intero (media sulle fette, energia)", ("E", "FW1", "cellG"), "C2-d")):
    y = []
    for q in qe:
        try:
            y.append(100 * (A.piR(*args, 0.0125, q) / pitR_uehling(q) - 1))
        except Exception:
            y.append(np.nan)
    b.semilogx(qe, y, st, ms=5, label=lab)
b.axhline(0, color="gray", lw=0.8)
b.axhspan(-0.5, 0.5, color="gray", alpha=0.12)
b.set_xlabel("q / m")
b.set_ylabel("Π_R / Uehling − 1  [%]")
b.set_title("(b) forma di Uehling, θ = m a = 0,0125 (banda: soglia 0,5 %)")
b.legend(fontsize=8)
b.grid(alpha=0.3)

# (c) charge renormalization log
c = ax[1, 0]
x = -np.log(TH)
ex = json.load(open(os.path.join(os.path.dirname(OUT), "fw1vuoto_exams.json")))
p0 = [ex["E3"]["pit0_avg"][str(t)] for t in TH]
c.plot(x, np.array(p0) / L, "C0s-", label="FW1 mare intero: elettrico (fine periodo)")
pe = [np.mean(list(ex["E5"]["pit_e0_slice_avg"][str(t)].values())) for t in TH]
pb = [np.mean(list(ex["E5"]["pit_b0"][str(t)].values())) for t in TH]
c.plot(x, np.array(pe) / L, "C2d-", label="FW1 mare intero: elettrico (media fette)")
c.plot(x, np.array(pb) / L, "C3v-", label="FW1 mare intero: magnetico")
pf = [ex["E4"]["F0_proj_pit0"][str(t)] for t in TH]
c.plot(x, np.array(pf) / L, "ko--", label="congelata + proiettore (Debiti)")
xx = np.linspace(x[0], x[-1], 10)
c.plot(xx, xx - x[0] + p0[0] / L, "C0:", lw=1, label="pendenza 1 (= 2α/3π)")
c.set_xlabel("ln(1/θ) = ln(1/(m a))")
c.set_ylabel("Π(0)  [unità 1/(6π²)]")
c.set_title("(c) logaritmo della carica: una specie sola")
c.legend(fontsize=8)
c.grid(alpha=0.3)

# (d) dressed light speed
d = ax[1, 1]
dc = ex["E5"]["delta_c_over_alpha"]
keys = [k for k in dc[str(TH[0])] if k != "iso"]
for k in keys:
    d.plot(TH, [dc[str(t)][k] for t in TH], marker=".", lw=1, label=k.replace("n", "n=").replace("_e", ", E∥"))
d.plot(TH, [dc[str(t)]["iso"] for t in TH], "k-o", lw=2, label="media isotropa")
d.axhline(0.047, color="gray", ls="--", lw=1)
d.text(0.06, 0.049, "congelata, luce col taglio K12 (Debiti f4)", fontsize=7, color="gray")
d.axhline(0, color="gray", lw=0.6)
d.set_xscale("log")
d.set_xlabel("θ = m a")
d.set_ylabel("(c_luce vestita − c)/(α c)")
d.set_title("(d) luce vestita dal mare intero di FW1 (statico)")
d.legend(fontsize=7, ncol=2)
d.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(OUT, dpi=110)
dst = "/mnt/project-files/problemi-aperti"
if os.path.isdir("/mnt/project-files"):
    os.makedirs(dst, exist_ok=True)
    shutil.copy(OUT, os.path.join(dst, "fw1vuoto_figura.png"))
print("written", OUT)
