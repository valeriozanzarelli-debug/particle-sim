"""aggmol_h2_figure.py - figure of phase 2 (H2 with two electrons): lab/out/aggmol_h2_figure.png
(a) Born-Oppenheimer curve with the vibrational levels; (b) exams vs thresholds (|relative or absolute error| / threshold)."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
r = json.load(open(os.path.join(OUT, "aggmol_h2.json")))
EH_EV = 27.211386245988
cur = r["curve"]
Rs = np.array(sorted(v["R"] for v in cur.values()))
Et = np.array([cur[f"{R:.2f}"]["E_tot"] for R in Rs])
EH = -1.0000000010921748 / 2
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
a = ax[0]
x = np.linspace(Rs[0], Rs[-1], 600)
cs = CubicSpline(Rs, Et)
a.plot(x, (cs(x) - 2 * EH) * EH_EV, color="tab:blue", label="curva (ECG a due centri, 2 elettroni)")
a.plot(Rs, (Et - 2 * EH) * EH_EV, "o", color="tab:blue", ms=4)
for i, lv in enumerate(r["vib"]["levels_Eh"]):
    e = (lv - 2 * EH) * EH_EV
    xs = x[(cs(x) - 2 * EH) * EH_EV <= e]
    a.hlines(e, xs.min(), xs.max(), color="tab:red", lw=1)
    a.text(xs.max() + 0.05, e, f"v={i}", fontsize=8, color="tab:red", va="center")
if "H5" in r:
    a.plot([1.4], [(r["H5"]["E_HF_grid"] - 2 * EH) * EH_EV], "s", color="tab:green", label="Hartree-Fock sulla griglia a0/4")
a.axhline(0, color="gray", lw=0.6)
a.set_xlabel("distanza fra i protoni R (a0)"); a.set_ylabel("E - 2E(H) (eV)")
a.set_ylim(-5.0, 1.5); a.set_xlim(0.7, 6.1)
v = r["vib"]
a.set_title(f"H2: R_e {v['R_e']:.4f} a0, D_e {v['D_e_eV']:.4f} eV, omega_e {v['omega_e_cm']:.0f} cm-1", fontsize=10)
a.legend(fontsize=8, loc="lower right")
b = ax[1]
V = r["verdict"]
items = [("H0 E(1.4)", abs(V["H0"]["err_mEh"]) / 0.05),
         ("H1 R_e", abs(V["H1"]["rel"]) / 3e-3),
         ("H2 D_e", abs(V["H2"]["diff_meV"]) / 43.4),
         ("H3 omega_e", abs(V["H3"]["rel"]) / 0.01),
         ("H4 nu01", abs(V["H4"]["rel"]) / 0.01),
         ("H4 D_0", abs(V["H4"]["D0_diff_meV"]) / 43.4)]
if "H5" in r:
    items.append(("H5 HF griglia", abs(r["H5"]["err_mEh"]) / 1.6))
if "G1" in r:
    items.append(("G1 H atomo", abs(r["G1"]["H_atom"]["err"]) / 1e-6))
    items.append(("G1 H2+", abs(r["G1"]["H2plus_R2"]["err"]) / 1e-6))
names = [i[0] for i in items]; vals = [max(i[1], 1e-3) for i in items]
b.barh(names, vals, color=["tab:green" if q <= 1 else "tab:red" for q in vals])
b.axvline(1, color="k", lw=1); b.set_xscale("log")
b.set_xlabel("|errore| / soglia pre-registrata  (< 1 = PASS)")
b.set_title("esami contro valori noti (solo esame, mai ingresso)", fontsize=10)
b.invert_yaxis()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "aggmol_h2_figure.png"), dpi=110)
print("ok")
