"""fw1atomo_figura.py - figure of the study fw1atomo_: lab/out/fw1atomo_figura.png (and a copy in
/mnt/project-files/problemi-aperti/).  Reads lab/out/fw1atomo_p1.json, _p2.json, _p3.json."""
import os, sys, json, shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
import zz_core as zz

OUT = os.path.join(os.path.dirname(__file__), "out")
P1 = json.load(open(os.path.join(OUT, "fw1atomo_p1.json")))
P2 = json.load(open(os.path.join(OUT, "fw1atomo_p2.json")))
P3 = json.load(open(os.path.join(OUT, "fw1atomo_p3.json")))
ZA = 60 * zz.ALPHA
OMD = zz.dirac_formula(1, -1, ZA)


def pair(cfgname, lo=0.84, hi=0.93):
    lv = [l for l in P1[cfgname]["levels"] if l["resid"] < 1e-5 and max(l["w_start"]) >= 0.01 and lo <= l["omega"] <= hi]
    return sorted(lv, key=lambda l: l["omega"])


fig, ax = plt.subplots(2, 2, figsize=(12, 9))
# (a) theta = 0,5, 48^3: levels per period-sublattice
a = ax[0, 0]
f0 = [l for l in P2["F0"]["ritz"]["levels"] if l["resid"] < 1e-6 and l["w"] >= 0.01]
for c in range(8):
    for l in f0:
        a.scatter(c - 0.15, l["omega"], c=[[1 - l["val_Gamma"], 0.2, l["val_Gamma"]]], marker="s", s=40)
for l in P2["FW1"]["ritzc"]:
    if l["resid"] < 1e-6 and 0.84 <= l["omega"] <= 0.93:
        up = l["spin"][2] > 0
        a.scatter(l["cls"] + 0.15, l["omega"], c=[[1 - l["val_Gamma"], 0.2, l["val_Gamma"]]], marker="^" if up else "v", s=60)
a.axhline(OMD, color="k", lw=0.8, ls="--")
a.text(7.6, OMD + 0.001, "Dirac", ha="right", fontsize=8)
a.set_xticks(range(8))
a.set_xlabel("sottoreticolo di periodo (i%2, j%2, k%2)")
a.set_ylabel("livello esatto  omega  [mc^2/hbar]")
a.set_title("(a) 1s, Z = 60, theta = 0,5, 48^3: quadrati congelata, triangoli FW1\n(colore: blu = valle vera, rosso = gemelli; punta = segno di S_z)", fontsize=9)

# (b) spin splitting inside one sublattice vs theta and vs h
b = ax[0, 1]
for c, mk in ((0, "o"), (1, "s")):
    xs, ys = [], []
    for name, aa in (("b_FW1_c%d" % c, 0.5), ("d_FW1_a0.400_c%d" % c, 0.4), ("d_FW1_a0.333_c%d" % c, 1 / 3), ("d_FW1_a0.250_c%d" % c, 0.25)):
        if name in P1:
            lv = pair(name)
            if len(lv) >= 2:
                xs.append(aa)
                ys.append(abs(lv[-1]["omega"] - lv[-2]["omega"]))
    b.plot(xs, ys, mk + "-", label=f"in theta (h = 0,4), sottoreticolo {c}")
    xs, ys = [], []
    for name, hh in (("c_h0.2_c%d" % c, 0.2), ("c_h0.3_c%d" % c, 0.3), ("b_FW1_c%d" % c, 0.4)):
        if name in P1:
            lv = pair(name)
            if len(lv) >= 2:
                xs.append(hh)
                ys.append(abs(lv[-1]["omega"] - lv[-2]["omega"]))
    b.plot(xs, ys, mk + ":", color="gray", label=f"in h (theta = 0,5), sottoreticolo {c}")
b.set_yscale("log")
b.set_xlabel("theta = a / lambda_C   (oppure h, linee grigie)")
b.set_ylabel("distanza fra le due righe dello stesso sottoreticolo")
b.axhline(8e-3, color="k", lw=0.6, ls="--")
b.text(0.25, 8.4e-3, "risoluzione di senzadoppi_ r7", fontsize=8)
b.legend(fontsize=7)
b.set_title("(b) le due righe = spin spaccato dalla T rotta: dipendenza da theta e da h", fontsize=9)

# (c) survival
c_ = ax[1, 0]
for rule, col in (("F0", "k"), ("FW1", "tab:blue")):
    p = np.array(P2[rule]["surv"]["p_exc_every_5"])
    t = np.arange(len(p)) * 5 * 1.0
    c_.plot(t, p, color=col, lw=1, label="congelata" if rule == "F0" else "FW1")
c_.set_xlabel("t  [lambda_C / c]")
c_.set_ylabel("p_exc = 1 - |<1s radiale|psi(t)>|^2")
c_.legend(fontsize=8)
c_.set_title("(c) battimento (48^3, theta = 0,5, partenza: 1s radiale del Lab)", fontsize=9)

# (d) 3-D atom: binding vs a^(2 gamma)
d = ax[1, 1]
gam = np.sqrt(1 - ZA**2)
BD = P3["B_dirac"]
for key, lab, mk in (("ritz_FW1_dirac_%s_c0", "FW1 (sottoreticolo 0)", "o"), ("ritz_FW1_dirac_%s_c1", "FW1 (sottoreticolo 1)", "x")):
    xs, ys = [], []
    for aa, e in P3["per_a"].items():
        k = key % aa
        if k in e and "main_B" in e[k]:
            xs.append(float(aa)**(2 * gam))
            ys.append(e[k]["main_B"] / BD - 1)
    d.plot(xs, np.array(ys) * 100, mk, label=lab)
xs, ys = [], []
for aa, e in P3["per_a"].items():
    k = "ritz_F0_dirac_%s_c0" % aa
    if k in e:
        for g in e[k]["groups"]:
            for om in g["omega"]:
                xs.append(float(aa)**(2 * gam))
                ys.append((1 - om) / BD - 1)
d.plot(xs, np.array(ys) * 100, "s", mfc="none", color="gray", label="congelata: tutti i livelli del multipletto")
xs = [float(aa)**(2 * gam) for aa in P3["per_a"]]
d.plot(xs, [(e["radial_B"] / BD - 1) * 100 for e in P3["per_a"].values()], "k+", ms=10, label="regola radiale (stesso a)")
X2 = os.path.join(OUT, "fw1atomo_x2.json")
if os.path.exists(X2):
    x2 = json.load(open(X2))
    for key, mk, lab in (("b20", "D", "FW1, scatola ~20 (x2)"), ("b24", "P", "FW1, scatola 24 (x2)")):
        pts = [(float(k.split("b")[0]), v["err"]) for k, v in x2.items() if (("b" not in k and key == "b20") or k.endswith(key))]
        if pts:
            d.plot([p[0]**(2 * gam) for p in pts], [100 * p[1] for p in pts], mk, mfc="none", ms=8, label=lab)
if "extrap_FW1_a2gamma" in P3:
    d.plot([0], [P3["extrap_FW1_a2gamma"]["rel"] * 100], "r*", ms=12, label="FW1 estrapolato (a^2gamma, scatola 16)")
d.axhline(0, color="k", lw=0.6)
d.set_xlabel("a^(2 gamma)   (gamma = 0,899)")
d.set_ylabel("errore sul legame dell'1s  [%]")
d.legend(fontsize=7)
d.set_title("(d) atomo 3-D, Z = 60: legame contro Dirac", fontsize=9)
fig.tight_layout()
png = os.path.join(OUT, "fw1atomo_figura.png")
fig.savefig(png, dpi=110)
os.makedirs("/mnt/project-files/problemi-aperti", exist_ok=True)
shutil.copy(png, "/mnt/project-files/problemi-aperti/fw1atomo_figura.png")
print("written", png)
