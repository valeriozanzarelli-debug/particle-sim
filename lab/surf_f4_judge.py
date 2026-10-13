"""surf_f4_judge.py -> lab/out/surf_f4_results.json + lab/out/surf_f4_figure.png : verdicts of G1-G6 (surf_f4_prereg.json)."""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
AS = (0.5, 0.4, 0.3, 0.25)
RUN = {a: json.load(open(os.path.join(OUT, f"surf_f4_run_{a:.2f}.json"))) for a in AS}
PREP = {a: json.load(open(os.path.join(OUT, f"surf_f4_prep_{a:.2f}.json"))) for a in AS}
R = {}
R["G1"] = dict(value={a: PREP[a]["prepared"] for a in AS},
               ok=all(PREP[a]["prepared"]["valley0"] >= 0.99 and PREP[a]["prepared"]["sigma_z"] >= 0.9 for a in AS))
g2 = {a: dict(min_fid=min(RUN[a]["rec"]["fid"] + [RUN[a]["final_fid"]]), min_valley0=min(RUN[a]["rec"]["valley0"] + [RUN[a]["final_valley0"]])) for a in AS}
R["G2"] = dict(value=g2, ok=all(v["min_fid"] >= 0.9 and v["min_valley0"] >= 0.9 for v in g2.values()))
win = {a: {int(k): w for k, w in RUN[a]["windows"].items()} for a in AS}
g3 = {a: [round(w["vphi_over_Za_sin"], 4) for k, w in sorted(win[a].items())] for a in AS}
R["G3"] = dict(value=g3, ok=all(len(v) == 5 and all(0.75 <= x <= 1.25 for x in v) for v in g3.values()))
g4 = {a: max(abs(z) for z in RUN[a]["rec"]["z"] + [RUN[a]["final_z"]]) for a in AS}
R["G4"] = dict(value=g4, ok=all(v <= 4 for v in g4.values()))
v = {a: float(np.mean([w["vphi_over_Za_sin"] for k, w in sorted(win[a].items())][:2])) for a in AS}
d = {a: 1 - v[a] for a in AS}
R["G5"] = dict(value=dict(v=v, deficit=d), ok=d[0.5] > d[0.4] > d[0.3] > d[0.25])
x = np.array(AS); y = np.array([d[a] for a in AS])
A1 = np.vstack([np.ones(4), x]).T; c1, *_ = np.linalg.lstsq(A1, y, rcond=None)
A2 = np.vstack([np.ones(4), x**2]).T; c2, *_ = np.linalg.lstsq(A2, y, rcond=None)
r2 = lambda A, c: float(1 - np.sum((A @ c - y)**2) / np.sum((y - y.mean())**2))
pw = None
if np.all(y > 0):
    cp, *_ = np.linalg.lstsq(np.vstack([np.ones(4), np.log(x)]).T, np.log(y), rcond=None)
    pw = dict(C=float(np.exp(cp[0])), p=float(cp[1]), R2=float(1 - np.sum((cp[0] + cp[1] * np.log(x) - np.log(y))**2) / np.sum((np.log(y) - np.log(y).mean())**2)))
R["G6"] = dict(value=dict(linear=dict(d0=float(c1[0]), c1=float(c1[1]), R2=r2(A1, c1)), quadratic=dict(d0=float(c2[0]), c2=float(c2[1]), R2=r2(A2, c2)),
                          power=pw), ok=abs(c1[0]) <= 0.05)
n = sum(r["ok"] for r in R.values())
json.dump(dict(exams=R, n_pass=n, n_fail=len(R) - n, prep=PREP), open(os.path.join(OUT, "surf_f4_results.json"), "w"), indent=1, default=float)
for k, r in R.items():
    print(k, "PASS" if r["ok"] else "FAIL", r["value"] if k != "G1" else "")
print("PASS", n, "FAIL", len(R) - n)

fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
cols = dict(zip(AS, ("tab:blue", "tab:orange", "tab:green", "tab:purple")))
for a in AS:
    r = RUN[a]["rec"]; t = np.array(r["t"]) * a
    ax[0].plot(t, r["flux_vphi"], color=cols[a], lw=1, label=f"a = {a}: corrente dello stato")
    w = sorted(win[a].items())
    ax[0].plot([k * a for k, _ in w], [q["vphi_over_Za_sin"] for _, q in w], "o", color=cols[a], label=f"a = {a}: surfisti")
    ax[1].plot(t, r["fid"], color=cols[a], label=f"a = {a}")
ax[0].axhspan(0.75, 1.25, color="tab:green", alpha=0.12); ax[0].axhline(1, color="k", lw=0.6, ls="--")
ax[0].set_xlabel("tempo (lambda_C / c)"); ax[0].set_title("rotazione v_phi / (Z alpha c sin theta)"); ax[0].legend(fontsize=6)
ax[1].set_xlabel("tempo (lambda_C / c)"); ax[1].set_title("fedelta' |<psi0|psi(t)>|^2"); ax[1].legend(fontsize=7)
xx = np.linspace(0, 0.55, 50)
ax[2].plot(x, y, "ko", ms=8, label="mancanza d(a) = 1 - v(a)")
ax[2].plot(xx, c1[0] + c1[1] * xx, "b-", label=f"d0 + c a (d0 = {c1[0]:.3f})")
ax[2].plot(xx, c2[0] + c2[1] * xx**2, "r--", label=f"d0 + c a^2 (d0 = {c2[0]:.3f})")
ax[2].axhline(0, color="k", lw=0.5); ax[2].axhspan(-0.05, 0.05, color="tab:green", alpha=0.12)
ax[2].set_xlabel("a (lambda_C)"); ax[2].set_title("quanto manca a Dirac, contro il passo"); ax[2].legend(fontsize=7)
fig.suptitle("Fase 4: rotazione dell'1s sulla valle k = 0 con il nucleo trasferito (Z = 60)")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "surf_f4_figure.png"), dpi=110)
