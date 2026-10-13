"""nucleo1s_figure.py - figure of the study NUCLEO_1S -> lab/out/nucleo1s_figure.png (from nucleo1s_results.json)"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
R = json.load(open(os.path.join(OUT, "nucleo1s_results.json")))
P = json.load(open(os.path.join(OUT, "nucleo1s_prereg.json")))
fig, ax = plt.subplots(2, 2, figsize=(12, 9))

# (a) the law at Z = 60 and other Z
a_ = ax[0, 0]
cols = {"20": "#4C72B0", "40": "#55A868", "60": "#C44E52", "80": "#8172B2", "92": "#CCB974"}
for Z, c in cols.items():
    d = R["law"][Z]["1,-1"]
    a = np.array(d["a"]); e = np.array(d["err"])
    pz = P["predictions_law"]["values"][Z]
    aa = np.logspace(np.log10(a.min()), np.log10(0.5), 200)
    a_.loglog(a, np.abs(e), "o", color=c, ms=4, label=f"Z = {Z}")
    a_.loglog(aa, np.abs(pz["K"] * aa**pz["two_gamma"] + pz["B"] * aa**2), "-", color=c, lw=1)
a_.set_xlabel("cella a vicino al nucleo [λ_C]"); a_.set_ylabel("|errore relativo 1s|")
a_.set_title("(a) legge prevista K a^{2γ} + B a² (linee, 0 parametri) vs regola (punti)")
a_.legend(fontsize=8)

# (b) cost: cells vs error, uniform vs automatic zoom (Z = 60)
b_ = ax[0, 1]
d = R["law"]["60"]["1,-1"]
rmax = 30 + 9 / (60 * 7.2973525693e-3)
a = np.array(d["a"]); e = np.array(d["err"])
b_.loglog(np.ceil(rmax / a), np.abs(e), "s-", color="gray", label="celle uguali")
zz = [(v["cells"], abs(v["err_1s"]), k) for k, v in R["zoom"].items() if k.startswith("budget")]
for cells, err, k in zz:
    b_.loglog(cells, err, "o", color="#C44E52", ms=8)
    b_.annotate(k.replace("budget_", "budget "), (cells, err), fontsize=8, xytext=(5, 5), textcoords="offset points")
geo = [g for g in R.get("posthoc_geo", []) if g["hmax"] == 4.0 and g["hmin"] == 1e-3]
if geo:
    b_.loglog([g["cells"] for g in geo], [abs(g["err"]) for g in geo], "^-", color="#55A868",
              label="celle ∝ r (idea geometrica, post hoc)")
b_.axhline(0.01, color="k", ls=":", lw=1)
b_.loglog([], [], "o", color="#C44E52", label="zoom automatico")
b_.set_xlabel("numero di celle (radiali)"); b_.set_ylabel("|errore 1s| (Z = 60)")
b_.set_title("(b) stesso errore con molte meno celle"); b_.legend(fontsize=8)

# (c) where the +39.9 % comes from (regions)
c_ = ax[1, 0]
reg = R["decomp"]["regions"]
keys = list(reg.keys())
vals = [R["decomp"]["lab"]] + [reg[k] for k in keys]
labels = ["tutto 0,5"] + [f"fini r<{float(k):.3g}" if float(k) < 100 else "tutto fine" for k in keys]
c_.bar(range(len(vals)), 100 * np.array(vals), color=["#C44E52"] + ["#4C72B0"] * len(keys))
c_.set_xticks(range(len(vals))); c_.set_xticklabels(labels, rotation=35, fontsize=8)
c_.set_ylabel("errore 1s [%] (Z = 60, celle 0,5 fuori)")
c_.set_title("(c) il +39,9 %: quasi tutto nella prima cella")
for i, v in enumerate(vals):
    c_.text(i, 100 * v + 0.5, f"{100 * v:.2f}", ha="center", fontsize=7)

# (d) Z = 1: contact removed by the zoom, what is left is the kinetic mass theta^2/3
d_ = ax[1, 1]
z1 = R["z1"]
th = np.array([0.05, 0.025, 0.0125])
from_doc = {0.2: 1.30e-2, 0.1: 3.08e-3, 0.05: 7.10e-4, 0.025: 1.54e-4, 0.0125: 3.07e-5}   # ZIGZAG.md 5.4 (uniform)
tu = np.array(sorted(from_doc)); eu = np.array([from_doc[t] for t in tu])
d_.semilogx(tu, eu / tu**2, "s-", color="gray", label="celle uguali (ZIGZAG §5.4)")
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_core import dirac_binding, ALPHA
ez = [z1["levels"][str(t)]["1"] / dirac_binding(1, -1, ALPHA) - 1 for t in th]
if True:
    d_.semilogx(th, np.array(ez) / th**2, "o-", color="#C44E52", label="con lo zoom al nucleo")
d_.axhline(1 / 3, color="k", ls=":", lw=1, label="θ²/3 (massa cinetica della regola)")
d_.set_xlabel("θ = cella [λ_C]"); d_.set_ylabel("errore 1s / θ²  (Z = 1)")
d_.set_title("(d) idrogeno: il 'θ² ln θ' è il contatto; lo zoom lo toglie"); d_.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "nucleo1s_figure.png"), dpi=130)
print("saved")
