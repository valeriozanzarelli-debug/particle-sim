"""rifl_r5.py - il limite superiore rigoroso "per tick" e il suo controllo.
Teorema (Cauchy-Schwarz): per ogni stato della banda positiva con k in W,
   |L_n(1)|^2 = |int_W f e^{-iwn} uL e^{ik} dk/sqrt(2pi)|^2 <= int|f|^2 * a(theta,W),  a = int_W |uL|^2 dk/2pi
quindi  Delta_N <= sum_n |L_n(1)|^2 <= N a(theta, W).
Continuo di Dirac: |L(0,t)|^2 <= int_0^inf L(p)^2 dp/2pi = m/(4 pi)  =>  Delta(T) <= m c^2 T/(4 pi hbar) = 1/(pi eps^2).
Qui: a(theta,W) per la regola, controllo Delta <= N a su tutti i dati Q2, rapporto di saturazione per Dirac.
uso: python lab/rifl_r5.py   (< 1 min)"""
import sys, os, json
import numpy as np
from scipy.integrate import quad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rifl_core import *

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
r3 = json.load(open(os.path.join(OUT, "rifl_r3.json")))
r1 = json.load(open(os.path.join(OUT, "rifl_r1.json")))


def a_of(theta, kmax):
    f = lambda k: np.abs(walk_modes(np.array([k]), theta)[2][0])**2 / (2 * np.pi)
    return quad(f, 0, kmax, limit=400, epsabs=1e-15, epsrel=1e-12, points=[min(theta, kmax / 2)])[0]


res = {"teorema": __doc__.split("\n")[1:6]}
tab = {}
for th in [1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.3, 0.6, 1.0, 1.3, 1.5]:
    ap, af = a_of(th, np.pi / 2), a_of(th, np.pi)
    l1p = walk_backflow_gram(th, 1, np.pi / 2)
    l1f = walk_backflow_gram(th, 1, np.pi)
    tab[f"{th:g}"] = {"a_fisica": ap, "a_intera": af, "a_fisica_su_theta/4pi": ap / (th / (4 * np.pi)),
                      "Lambda1_fisica": l1p, "Lambda1_intera": l1f, "Lambda1/a fisica": l1p / ap}
    print(th, tab[f"{th:g}"], flush=True)
res["per_theta"] = tab
viol = []
for th, rows in r3["Q2_theta"].items():
    ap, af = a_of(float(th), np.pi / 2), a_of(float(th), np.pi)
    for N, v in rows.items():
        if v["fisica"] > int(N) * ap * (1 + 1e-9) or v["intera"] > int(N) * af * (1 + 1e-9):
            viol.append((th, N))
res["violazioni_Delta_le_N_a"] = viol
sat = {}
q1 = r3.get("Q1_eps", {})
for e, v in q1.items():
    sat[e] = v["c_rbf"] * np.pi * float(e)**2
for e, v in r1["E1"]["per_eps"].items():
    sat.setdefault(f"{float(e):g}", v["c_rbf"] * np.pi * float(e)**2)
res["Dirac_saturazione_c_su_limite"] = dict(sorted(sat.items(), key=lambda z: float(z[0])))
res["Dirac_limite_rispettato"] = bool(all(s <= 1 + 1e-6 for s in sat.values()))
json.dump(res, open(os.path.join(OUT, "rifl_r5.json"), "w"), indent=1)
print("violazioni:", viol, "Dirac rispetta il limite:", res["Dirac_limite_rispettato"])
print(res["Dirac_saturazione_c_su_limite"])
