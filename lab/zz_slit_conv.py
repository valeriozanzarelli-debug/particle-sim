"""zz_slit_conv.py - lattice convergence of the double slit made of light walkers (added after the
first double-slit run failed, to find out why; information, no pass/fail).

Cause found before this run (dbg, see docs/ZIGZAG.md): on the lattice the two-tick light walk has
cos(omega dt) = cos(k_x a) cos(k_z a) in the x-z plane, so a beam along z spreads with the transverse
curvature  d^2 omega / d k_x^2 = cot(k_z a)/a  instead of 1/k_z: diffraction is too weak by the factor
(k_z a) cot(k_z a) ~ 1 - (k a)^2/3  (lambda = 8a: 0.785 -> beams 20 % too narrow).
Same geometry as zz_slit.py but screen at z = 40 lambda (|x| < 24), so lambda/a = 8, 16, 24 fit in the
time budget. Reference: exact 2-D Rayleigh-Sommerfeld flux at z = 40.
Run: timeout 590 python lab/zz_slit_conv.py
"""
import json, os, time
import numpy as np
from zz_slit import double_slit, OUT

t0 = time.time()
res = {}
pats = {}
for per_lam in (8, 16, 24):   # (the first run reached 8 and 16; 24 was cut by the 590 s limit)
    if per_lam == 24 and time.time() - t0 > 150:
        print("lambda/a = 24 skipped (time)", flush=True)
        continue
    r, p, ref = double_slit(per_lam, Ls=40.0, xhalf=31.0, half=24.0)
    k = 2 * np.pi / per_lam
    r["diffraction_factor_kcotk"] = float(k / np.tan(k))
    res[str(per_lam)] = r
    pats[per_lam] = p
    print(f"lambda/a = {per_lam}: corr(fine bins) with exact RS flux at z = 40: {r['corr_fine']:.4f}, coarse {r['corr_coarse']:.4f}"
          f"  [(ka)cot(ka) = {k / np.tan(k):.4f}]  ({r['seconds']:.0f} s)", flush=True)
    json.dump(res, open(os.path.join(OUT, "zz_slit_conv.json"), "w"), indent=1, default=float)
json.dump(res, open(os.path.join(OUT, "zz_slit_conv.json"), "w"), indent=1, default=float)
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xs = np.linspace(-24, 24, 97)[:-1] + 0.25
    fig, ax = plt.subplots(figsize=(8, 3.4))
    ax.plot(xs, ref, "k-", lw=1.4, label="exact Rayleigh-Sommerfeld flux, z = 40")
    for pl, c in ((8, "r"), (16, "tab:orange"), (24, "tab:green")):
        if pl in pats:
            ax.plot(xs, pats[pl], "-", color=c, lw=0.9, label=f"light walk, lambda = {pl}a")
    ax.set_xlabel("x / lambda"); ax.legend(fontsize=7); ax.set_title("double slit: lattice convergence")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "zz_slit_conv.png"), dpi=90)
except Exception as ex:
    print("figure skipped:", ex)
print(f"total {time.time() - t0:.0f} s")
