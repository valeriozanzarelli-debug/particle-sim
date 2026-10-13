"""supercrit_figure.py - figure of the study: lab/out/supercrit_figura.png"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
L = lambda n: json.load(open(os.path.join(HERE, "out", f"supercrit_{n}.json")))
from supercrit_prereg import GODUNOV
c1, c2, g1, ex = L("c1"), L("c2"), L("g1"), L("exams")
fig = plt.figure(figsize=(14, 9.5))
# (a) point nucleus
ax = fig.add_subplot(2, 3, 1)
rows = sorted([v for k, v in c1.items() if k[0].isdigit()], key=lambda x: x["a"])
A = np.array([r["a"] for r in rows]); Z = np.array([r["Zcr"] for r in rows])
ax.semilogx(A, Z, "o-", color="C0", label="regola, nucleo puntiforme")
ax.axhline(137.036, color="k", ls="--", lw=1, label="Zα = 1 (137,04)")
f = ex["C1_fit"]
aa = np.logspace(-6, -1, 100)
ax.semilogx(aa, f["Z_star"] + f["c"] * aa**f["power"], ":", color="C3",
            label=f"Z* + c a^q: Z* = {f['Z_star']:.2f}, q = {f['power']:.2f}")
ax.set_xlabel("passo del reticolo a [ħ/mc]"); ax.set_ylabel("Z in cui l'1s tocca −mc²")
ax.set_title("(a) nucleo puntiforme: il reticolo si tuffa a 131,4"); ax.legend(fontsize=8)
# (b) levels, shell nucleus
ax = fig.add_subplot(2, 3, 2)
tr = np.array(c2["shell"]["0.00025"]["track"])
ax.plot(tr[:, 0], tr[:, 1], "-", color="C0", label="regola, a = 0,00025 (legato)")
Zr = sorted(int(z) for z in g1)
ax.plot(Zr, [g1[str(z)]["extrap"]["E_res_quad3"] for z in Zr], "o", color="C0", label="regola a→0 (risonanza)")
ax.plot(list(GODUNOV), [v[0] for v in GODUNOV.values()], "x", color="k", label="Godunov et al. 2017")
ax.axhline(-1, color="gray", lw=1); ax.text(200, -0.9, "mare di Dirac (−mc²)", fontsize=8, color="gray")
ax.set_xlim(158, 234); ax.set_xlabel("Z"); ax.set_ylabel("energia dell'1s [mc²]")
ax.set_title("(b) nucleo esteso R = 12 fm: Z_cr = 177,463"); ax.legend(fontsize=8)
# (c) widths
ax = fig.add_subplot(2, 3, 3)
pw = L("posthoc_wide")
Gw = [g1[str(z)]["extrap"]["Gamma_quad3"] if str(z) not in pw or "walk_win3_extrap" not in pw[str(z)]
      else pw[str(z)]["walk_win3_extrap"]["Gamma"] for z in Zr]
ax.semilogy(Zr, np.array(Gw) * 511e3, "o", color="C0", label="regola (fase), a→0")
zg = [z for z in GODUNOV if GODUNOV[z][1] > 0]
ax.semilogy(zg, [GODUNOV[z][1] * 511e3 for z in zg], "x", color="k", label="Godunov et al. 2017")
g2 = L("g2"); pd = L("posthoc_decay")
ax.semilogy([int(z) for z in pd], [pd[z]["Gamma_fit_floor"] * 511e3 for z in pd], "s", mfc="none", color="C3",
            label="regola (decadimento nel tempo), a = 0,002")
ax.set_xlabel("Z"); ax.set_ylabel("larghezza Γ [eV]"); ax.set_title("(c) larghezza della risonanza")
ax.legend(fontsize=8)
# (d) the positron leaving: 2-D slice of the density, Z = 200
d = np.load(os.path.join(HERE, "out", "supercrit_g2_200.npz"))
r, snaps, dt = d["r"], d["snaps"], float(d["snap_dt"])
x = np.linspace(-45, 45, 500)
X, Y = np.meshgrid(x, x)
R = np.sqrt(X**2 + Y**2)
for j, k in enumerate((0, 3, 6, 11)):
    ax = fig.add_axes([0.04 + j * 0.24, 0.05, 0.21, 0.36])
    img = np.interp(R, r, snaps[k] / (r[1] - r[0]))     # weight per unit radius (shows the outgoing shell)
    vmax = np.sqrt(np.max(snaps[:, r > 2.0]) / (r[1] - r[0]))
    ax.imshow(np.sqrt(img), extent=[-45, 45, -45, 45], cmap="magma", origin="lower", vmin=0, vmax=vmax)
    ax.add_patch(plt.Circle((0, 0), 40, fill=False, ls=":", color="w", lw=0.6))
    ax.set_xticks([]); ax.set_yticks([])
    t = (k + 1) * dt
    ax.set_title(f"(d) Z = 200, t = {t:.0f} ħ/mc²  (Γt = {t * g2['200']['Gamma_phase']:.1f})", fontsize=9)
    if j == 0:
        ax.text(-43, -43, "carica che esce = positrone (+e)\nal centro resta il vuoto carico (−e)\n(colore = peso per raggio, scala comune;\n cerchio = inizio della spugna)", color="w", fontsize=7)
out = os.path.join(HERE, "out", "supercrit_figura.png")
fig.savefig(out, dpi=110)
print("written", out)
