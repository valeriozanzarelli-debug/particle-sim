"""nucleo1s_3d_figure.py - figure of NUCLEO_1S phase 2 (3-D) -> lab/out/nucleo1s_3d_figure.png"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA, dirac_formula  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
R = json.load(open(os.path.join(OUT, "nucleo1s_3d_results.json")))
J = R["judge"]
Za = 60 * ALPHA
D1 = dirac_formula(1, -1, Za)
fig, ax = plt.subplots(2, 2, figsize=(12.5, 9))

# (a) the Laboratorio 1s: four lines, and the beat seen in the Laboratorio
a_ = ax[0, 0]
L = R["lab"]["(0, 0, 0)"]["lines"]
for x in L:
    a_.vlines(x["omega"], 0, x["weight"], color="#C44E52", lw=3)
lr = R.get("scan_learned_0.5")
if lr:
    for x in lr["lines"]:
        a_.vlines(x["omega"] + 0.0004, 0, x["weight"], color="#55A868", lw=3)
a_.axvline(D1, color="k", ls=":", lw=1)
a_.text(D1, 0.95, " Dirac", fontsize=8)
for pk, (y, c) in zip(R["lab"]["lab_pexc_peaks"], ((0.6, "#4C72B0"), (0.7, "#8172B2"))):
    hi = max(x["omega"] for x in L)
    a_.annotate("", xy=(hi - pk["omega"], y), xytext=(hi, y), arrowprops=dict(arrowstyle="<->", color=c))
    a_.text(hi - pk["omega"], y + 0.02, f"battimento del Lab {pk['period_ticks']:.0f} tick", fontsize=8, color=c)
a_.set_xlim(0.83, 0.915); a_.set_ylim(0, 1.0)
a_.set_xlabel("ω (mc²/ħ)"); a_.set_ylabel("peso nello stato di partenza del Lab")
a_.set_title("(a) l'1s del Laboratorio sono 4 livelli (rosso); nucleo imparato: 1 (verde)")

# (b) lines vs a
b_ = ax[0, 1]
T = J["table"]
aa = sorted(float(a) for a in T)
for a in aa:
    r = T[str(a)]
    for o in r.get("L", []):
        b_.plot(a, 1 - o, "o", color="#C44E52", ms=4)
    if "learned_centroid" in r:
        b_.plot(a, 1 - r["learned_centroid"], "s", color="#55A868", ms=6)
b_.axhline(1 - D1, color="k", ls=":", lw=1, label="Dirac 1s")
b_.plot([], [], "o", color="#C44E52", label="regola nuda: le righe dell'1s")
b_.plot([], [], "s", color="#55A868", label="nucleo imparato (2 numeri)")
b_.set_xlabel("cella a (λ_C)"); b_.set_ylabel("legame (mc²)")
b_.set_title("(b) il multipletto si chiude come a²; il nucleo imparato lo chiude ma lega poco")
b_.legend(fontsize=8)

# (c) error budget
c_ = ax[1, 0]
A = np.array(aa)
c_.plot(A, [100 * T[str(a)]["bare_rel"] for a in aa], "o-", color="#C44E52", label="nudo: centro del multipletto")
c_.plot(A, [100 * T[str(a)]["U_rel"] for a in aa], "^-", color="#4C72B0", label="nudo: coppia cieca al nucleo (U)")
c_.plot(A, [100 * T[str(a)]["valley_part"] for a in aa], "--", color="#C44E52", alpha=0.6, label="parte 'valli' = centro − U")
la = [a for a in aa if "learned_rel" in T[str(a)]]
c_.plot(la, [100 * T[str(a)]["learned_rel"] for a in la], "s-", color="#55A868", label="nucleo imparato")
c_.axhspan(-1, 1, color="0.9")
c_.axhline(0, color="k", lw=0.5)
c_.set_xlabel("cella a (λ_C)"); c_.set_ylabel("errore del legame 1s (%)")
c_.set_title("(c) bilancio: quasi tutto l'errore 3-D è lo split delle valli")
c_.legend(fontsize=7)

# (d) the glue
d_ = ax[1, 1]
gl = R.get("glue")
if gl:
    js = [3, 4, 5, 6]
    for suf, lab, c, mk in (("_phase1", "colla fase 1", "gray", "o"), ("", "colla tau (pre-registrata)", "#C44E52", "s"),
                             ("_cusp", "colla del cuspide (emendamento)", "#55A868", "^")):
        d_.plot(js, [100 * gl["Z1_1s"][f"J0={j}{suf}"] for j in js], mk + "-", color=c, label=lab)
    d_.plot(js, [100 * gl["Z2_1s"][f"J0={j}_cusp"] for j in js], "^:", color="#55A868", alpha=0.6, label="cuspide, He⁺ 1s")
    d_.axhline(0, color="k", lw=0.5)
    d_.set_xlabel("J0 (celle grosse date dalla regola fine)"); d_.set_ylabel("1s Z = 1 − Dirac (%)")
    d_.set_title("(d) idrogeno ibrido: la colla non dipende più da J0")
    d_.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "nucleo1s_3d_figure.png"), dpi=130)
print("saved")
