"""zstar_figure.py - figure of ZSTAR: lab/out/zstar_figura.png (+ copy in /mnt/project-files/problemi-aperti/).
(a) the lattice zero mode near the nucleus: |g(rho)| for zeta below, at and above zeta* (r^+gamma vs r^-gamma);
(b) Z_cr(a): walk (points) vs matching formula (lines) for the base rule and variants;
(c) Z* of the formula vs offset s and reflection phase eta (where the false dive sits; none for s > s_c)."""
import json, os, shutil, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from zstar_core import scale_free_coins, _fwd
from zz_core import ALPHA

HERE = os.path.dirname(os.path.abspath(__file__))
PRE = json.load(open(os.path.join(HERE, "out", "zstar_prereg.json")))
SGN = json.load(open(os.path.join(HERE, "out", "zstar_runs_sign.json")))
C1 = json.load(open(os.path.join(HERE, "out", "supercrit_c1.json")))
PR = PRE["predictions"]
ZS = PR["base"]["zeta_roots"][0]

fig, ax = plt.subplots(1, 3, figsize=(16, 5))
# (a)
for z, col, lab in [(0.93, "#7f8c8d", "ζ = 0,93 (Z = 127)"), (ZS, "#c0392b", f"ζ* = {ZS:.5f} (Z* = {ZS / ALPHA:.2f})"),
                    (0.985, "#2e86c1", "ζ = 0,985 (Z = 135)")]:
    rho, M = scale_free_coins(z, -1, 0.5, 200000)
    cp, cm, lsc = _fwd(M, -1.0 + 0j, 1.0 + 0j)
    g = np.abs(cp + cm) / np.sqrt(2) * np.exp(lsc)
    ii = np.unique(np.logspace(0, np.log10(len(rho) - 1), 400).astype(int)); ii = np.concatenate([[0], ii])
    ax[0].loglog(rho[ii], g[ii], color=col, label=lab)
    gam = np.sqrt(1 - z * z)
    if abs(z - ZS) < 1e-9:
        c0 = g[1000] * rho[1000]**gam
        ax[0].loglog(rho[100::5000], 1.3 * c0 * rho[100::5000]**(-gam), "k:", lw=1, label="∝ r^−γ (×1,3)")
ax[0].set_xlabel("r / a  (celle dal nucleo)")
ax[0].set_ylabel("|g|  (soluzione regolare della regola, E = 0, m = 0)")
ax[0].set_title("(a) il modo di contatto: a ζ* la regola sceglie solo r^−γ")
ax[0].legend(fontsize=8)
# (b)
A = PRE["configs"]["a"]
cols = {"base": "#c0392b", "R+1": "#8e44ad", "R-i": "#d35400", "S0.25": "#16a085", "K+1": "#2c3e50",
        "S1": "#2e86c1", "S0.866": "#5dade2"}
for v, col in cols.items():
    aa = A + [1e-6, 1e-9]
    zm = [PR[v]["Zcr_matching"][f"{a:.6g}"] for a in aa]
    ax[1].semilogx(aa, zm, "-", color=col, lw=1)
    AW = [3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]
    zw = [SGN.get(f"{v}|{a:.6g}", {}).get("Zcr") for a in AW]
    zw = [np.nan if z is None else z for z in zw]
    ax[1].semilogx(AW, zw, "o", color=col, ms=5, label=v if v != "base" else "base (s = ½, refl0 = −1)")
ax[1].axhline(1 / ALPHA, color="k", lw=0.8, ls="--")
ax[1].text(2e-9, 137.6, "Z = 1/α = 137,04", fontsize=8)
ax[1].set_xlabel("passo a [ħ/mc]")
ax[1].set_ylabel("Z_cr(a): l'1s tocca −mc²")
ax[1].set_title("(b) regola (punti) e formula di raccordo (linee)")
ax[1].legend(fontsize=7, ncol=2)
ax[1].set_ylim(100, 160)
# (c)
ss = sorted(PRE["s_scan_kappa-1"].items(), key=lambda x: float(x[0]))
sx = [float(k) for k, v in ss if v]
sy = [v[-1] / ALPHA for k, v in ss if v]
ax[2].plot(sx, sy, "o-", color="#16a085", label="Z*(s), refl0 = −1")
ax[2].axvline(PRE["s_c_kappa-1"], color="#16a085", ls=":", lw=1)
ax[2].text(PRE["s_c_kappa-1"] + 0.01, 102, f"s_c ≈ {PRE['s_c_kappa-1']:.3f}:\noltre, niente Z* finto", fontsize=8)
ee = sorted(PRE["eta_scan_kappa-1_s0.5"].items(), key=lambda x: float(x[0]))
ex_, ey = [], []
for k, v in ee:
    for z in v:
        ex_.append(float(k) / 2); ey.append(z / ALPHA)
ax2 = ax[2].twiny()
ax2.plot(ex_, ey, "s", color="#8e44ad", ms=4, label="Z*(η), s = ½")
ax2.set_xlabel("fase della riflessione η / 2π  (refl0 = e^{iη})", color="#8e44ad")
ax[2].axhline(1 / ALPHA, color="k", lw=0.8, ls="--")
ax[2].set_xlabel("offset s delle celle: r_i = (i + s) a", color="#16a085")
ax[2].set_ylabel("Z* (formula, a → 0)")
ax[2].set_ylim(20, 140)
ax[2].set_title("(c) Z* dipende dai dettagli al centro: è un artefatto di contorno", fontsize=10)
h1, l1 = ax[2].get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax[2].legend(h1 + h2, l1 + l2, fontsize=8, loc="lower left")
fig.tight_layout()
out = os.path.join(HERE, "out", "zstar_figura.png")
fig.savefig(out, dpi=110)
os.makedirs("/mnt/project-files/problemi-aperti", exist_ok=True)
shutil.copy(out, "/mnt/project-files/problemi-aperti/zstar_figura.png")
print("wrote", out)
