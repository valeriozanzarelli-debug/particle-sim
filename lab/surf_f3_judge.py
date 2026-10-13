"""surf_f3_judge.py -> lab/out/surf_f3_results.json + lab/out/surf_f3_figure.png : verdicts of F1-F6
(lab/out/surf_f3_prereg.json) and the phase-3 figure. Amendment C1 (surf_f3_amend.json): the 3990 window was never
measured (loop ends at t = 4000 before the window closes); F3/F5 judged on the 4 measured windows."""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "surf_f3_prereg.json")))
D = {k: json.load(open(os.path.join(OUT, f"surf_f3_run_{k}.json"))) for k in ("learned", "bare")}
PR = {k: json.load(open(os.path.join(OUT, f"surf_f3_prep_{k}.json"))) for k in ("learned", "bare")}
L, B = D["learned"], D["bare"]
win = lambda d: {int(t): w["vphi_over_Za_sin"] for t, w in d["windows"].items()}
R = {}
R["F1"] = dict(value=PR["learned"]["prepared"], ok=PR["learned"]["prepared"]["valley0"] >= 0.99 and PR["learned"]["prepared"]["sigma_z"] >= 0.9)
mf, mv = min(L["rec"]["fid"] + [L["final_fid"]]), min(L["rec"]["valley0"] + [L["final_valley0"]])
R["F2"] = dict(value=dict(min_fidelity=mf, min_valley0=mv), ok=mf >= 0.9 and mv >= 0.9)
wl = win(L)
R["F3"] = dict(value=wl, ok=len(wl) >= 4 and all(0.75 <= v <= 1.25 for v in wl.values()), note="4 of 5 windows measured (C1)")
mz = max(abs(z) for z in L["rec"]["z"] + [L["final_z"]])
R["F4"] = dict(value=mz, ok=mz <= 4)
wb = win(B)
out_b = sum(not (0.75 <= v <= 1.25) for v in wb.values())
R["F5"] = dict(value=dict(min_valley0=min(B["rec"]["valley0"]), windows=wb, windows_outside=out_b),
               ok=min(B["rec"]["valley0"]) < 0.5 or out_b >= 2)
R["F6"] = dict(value=L["rec"]["flux_vphi"][0], ok=0.75 <= L["rec"]["flux_vphi"][0] <= 1.25)
n = sum(r["ok"] for r in R.values())
json.dump(dict(exams=R, n_pass=n, n_fail=len(R) - n, prep=PR), open(os.path.join(OUT, "surf_f3_results.json"), "w"), indent=1, default=float)
for k, r in R.items():
    print(k, "PASS" if r["ok"] else "FAIL", r["value"])
print("PASS", n, "FAIL", len(R) - n)

fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
for k, c in (("learned", "tab:blue"), ("bare", "tab:red")):
    r = D[k]["rec"]; t = np.array(r["t"])
    lab = "nucleo imparato" if k == "learned" else "nucleo nudo (controllo)"
    ax[0].plot(t, r["valley0"], color=c, label=lab)
    ax[1].plot(t, r["sigma_z"], color=c, label=lab)
    ax[2].plot(t, r["flux_vphi"], color=c, lw=1, label=lab + ": corrente dello stato")
    w = win(D[k]); ax[2].plot(list(w), list(w.values()), "o", color=c, ms=8, label=lab + ": surfisti")
ax[0].set_title("peso nella valle k = 0"); ax[1].set_title("<Sigma_z>")
ax[2].axhspan(0.75, 1.25, color="tab:green", alpha=0.15, label="finestra Dirac (pre-registrata)")
ax[2].axhline(0, color="k", lw=0.5)
ax[2].set_title("rotazione v_phi / (Z alpha c sin theta)")
for A in ax:
    A.set_xlabel("tick"); A.legend(fontsize=7)
fig.suptitle("Fase 3: 1s polarizzato sulla valle k = 0 (48^3, Z = 60, a = 0,5)")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "surf_f3_figure.png"), dpi=110)
