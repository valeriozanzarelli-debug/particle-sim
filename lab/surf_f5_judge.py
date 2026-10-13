"""surf_f5_judge.py -> lab/out/surf_f5_results.json + lab/out/surf_f5_figure.png : verdicts of H1-H4 (surf_f5_prereg.json)."""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "surf_f5_prereg.json")))
run = json.load(open(os.path.join(OUT, "surf_f4_run_0.20.json")))
prep = json.load(open(os.path.join(OUT, "surf_f4_prep_0.20.json")))
win = sorted((int(k), w) for k, w in run["windows"].items())
v02 = float(np.mean([w["vphi_over_Za_sin"] for _, w in win[:2]]))
d02 = 1 - v02
pts = {float(k): v for k, v in P["fit_phase4"]["points"].items()}; pts[0.2] = d02
g2 = 2 * np.sqrt(1 - (60 / 137.035999084)**2)
x = np.array(sorted(pts)); y = np.array([pts[a] for a in x])
A = np.vstack([np.ones(len(x)), x**g2]).T
c, *_ = np.linalg.lstsq(A, y, rcond=None)
R2 = float(1 - np.sum((A @ c - y)**2) / np.sum((y - y.mean())**2))
others = {}
for name, p in (("a^2", 2.0), ("a", 1.0)):
    B = np.vstack([np.ones(len(x)), x**p]).T; cc, *_ = np.linalg.lstsq(B, y, rcond=None)
    others[name] = dict(d0=float(cc[0]), c=float(cc[1]), R2=float(1 - np.sum((B @ cc - y)**2) / np.sum((y - y.mean())**2)))
R = {}
R["H1"] = dict(value=dict(measured=d02, predicted=P["prediction"]["deficit"], diff=d02 - P["prediction"]["deficit"]),
               ok=abs(d02 - P["prediction"]["deficit"]) <= 0.005)
R["H2"] = dict(value=dict(d0=float(c[0]), c=float(c[1])), ok=abs(c[0]) <= 0.01)
R["H3"] = dict(value=R2, ok=R2 >= 0.999)
mf = min(run["rec"]["fid"] + [run["final_fid"]]); mz = max(abs(z) for z in run["rec"]["z"] + [run["final_z"]])
ws = [w["vphi_over_Za_sin"] for _, w in win]
R["H4"] = dict(value=dict(prepared=prep["prepared"], min_fid=mf, windows=ws, max_abs_z=mz),
               ok=prep["prepared"]["valley0"] >= 0.99 and prep["prepared"]["sigma_z"] >= 0.9 and mf >= 0.9
               and len(ws) == 3 and all(0.75 <= q <= 1.25 for q in ws) and mz <= 4)
n = sum(r["ok"] for r in R.values())
json.dump(dict(exams=R, n_pass=n, n_fail=len(R) - n, other_laws=others, v02=v02), open(os.path.join(OUT, "surf_f5_results.json"), "w"), indent=1, default=float)
for k, r in R.items():
    print(k, "PASS" if r["ok"] else "FAIL", r["value"])
print("other laws", others); print("PASS", n, "FAIL", len(R) - n)
fig, ax = plt.subplots(figsize=(7, 4.8))
xx = np.linspace(0, 0.55, 100)
ax.plot(x[x > 0.2], y[x > 0.2], "ko", ms=8, label="fase 4 (a = 0,5 / 0,4 / 0,3 / 0,25)")
ax.plot([0.2], [d02], "r*", ms=16, label="fase 5: a = 0,2 (misurato)")
ax.errorbar([0.2], [P["prediction"]["deficit"]], yerr=[0.005], fmt="s", color="tab:blue", mfc="none", capsize=4,
            label="previsione pre-registrata ± 0,005")
ax.plot(xx, c[0] + c[1] * xx**g2, "r-", label=f"d0 + c a^(2 gamma), 5 punti: d0 = {c[0]:.4f}")
ax.axhspan(-0.01, 0.01, color="tab:green", alpha=0.15, label="|d0| <= 1 %")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("a (lambda_C)"); ax.set_ylabel("quanto manca a Dirac: 1 - v_phi/(Z alpha c sin)")
ax.set_title("Fase 5: la legge del contatto per il giro dei surfisti"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "surf_f5_figure.png"), dpi=110)
