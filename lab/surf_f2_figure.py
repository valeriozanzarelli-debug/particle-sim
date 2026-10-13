"""surf_f2_figure.py -> lab/out/surf_f2_figure.png + lab/out/surf_f2_results.json : verdicts of exams V1-V7, S1-S5
(lab/out/surf_f2_prereg.json) from the measured files, and the phase-2 figure."""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "surf_f2_prereg.json")))
V = {(n, M): json.load(open(os.path.join(OUT, f"surf_f2_val_{n}_{M}.json"))) for n in (48, 96) for M in (1, 4, 16, 25)}
Q = json.load(open(os.path.join(OUT, "surf_f2_valli.json")))
OV = json.load(open(os.path.join(OUT, "surf_f2_overlap.json")))
PH = json.load(open(os.path.join(OUT, "surf_f2_posthoc.json"))) if os.path.exists(os.path.join(OUT, "surf_f2_posthoc.json")) else {}
cont = P["part1"]["exams"]["V6"]["continuum"]
R = {}


def f(n, M, s):
    return V[(n, M)]["fits"][s]


v1 = max(max(V[k]["H"]["born"]) / (3 * V[k]["floor"]) for k in V)
R["V1"] = dict(value=f"max H(born)/3F = {v1:.2f}", ok=v1 <= 1)
v2 = {f"n{n}_{s}": f(n, 1, s)["ratio"] for n in (48, 96) for s in ("rho2", "half")}
R["V2"] = dict(value=v2, ok=all(x >= 0.9 for x in v2.values()))
v3 = {f"n{n}_M{M}_{s}": f(n, M, s)["ratio"] for n in (48, 96) for M in (16, 25) for s in ("rho2", "uniform")}
R["V3"] = dict(value=v3, ok=all(x < 0.1 for x in v3.values()))
t4 = [f(48, M, "rho2")["tau"] for M in (1, 4, 16, 25)]
R["V4"] = dict(value=t4, ok=t4[0] > t4[1] > t4[2] >= t4[3])
v5 = {f"n{n}_M{M}": f(n, M, "rho2")["exp_R2"] for n in (48, 96) for M in (16, 25)}
R["V5"] = dict(value=v5, ok=all(x >= 0.9 for x in v5.values()),
               power_R2={f"n{n}_M{M}": f(n, M, "rho2")["pow_R2"] for n in (48, 96) for M in (16, 25)})
v6 = {f"n{n}_M{M}": (f(n, M, "rho2")["tau_sigma"], cont[str(M)]["tau_sigma"]) for n in (48, 96) for M in (16, 25)}
R["V6"] = dict(value=v6, ok=all(0.5 <= a / b <= 2 for a, b in v6.values()))
a, b = f(48, 16, "rho2")["tau_phys"], f(96, 16, "rho2")["tau_phys"]
R["V7"] = dict(value=(a, b), ok=abs(a - b) / max(a, b) <= 0.3)
ref = (0.8404, 0.8863, 0.8993, 0.8996)
lv = [L["omega"] for L in Q["ritz_levels"]]
R["S1"] = dict(value=lv, ok=len(lv) == 4 and all(abs(x - y) <= 1e-3 for x, y in zip(lv, ref)))
st = Q["states"]
w = {k: v["weight"] for k, v in st[Q["S5_level"]]["valleys"].items()}
R["S2"] = dict(value=w, E5_state={k: v["weight"] for k, v in st["E5_listened"]["valleys"].items()},
               ok=all(abs(x - 0.25) <= 0.05 for x in w.values()))
S3 = {}
for name in (Q["S5_level"], "E5_listened"):
    d = st[name]
    own = {k: (v["sigma_z_lab"] if k == "0" else -v["sigma_z_lab"]) for k, v in d["valleys"].items()}
    S3[name] = dict(sum_err=abs(d["sum_w_sigma"] - d["sigma_z"]), own_frame=own,
                    ok=abs(d["sum_w_sigma"] - d["sigma_z"]) <= 0.05 and
                    all(abs(own[k]) >= 0.9 for k, v in d["valleys"].items() if v["weight"] > 0.1))
R["S3"] = dict(value=S3, ok=all(x["ok"] for x in S3.values()))
sm = st.get("smooth_nucleus_1s")
s4 = dict(valley0=sm["valleys"]["0"]["weight"], sigma_z=sm["sigma_z"], vphi=Q["S4"]["vphi_over_Za_sin"])
R["S4"] = dict(value=s4, ok=s4["valley0"] >= 0.9 and s4["sigma_z"] >= 0.9 and 0.75 <= s4["vphi"] <= 1.25)
s5 = Q["S5"]
R["S5"] = dict(value=dict(level=Q["S5_level"], vphi=s5["vphi_over_Za_sin"], max_abs_z=s5["max_abs_z"], drift=s5["drift"]),
               ok=0.75 <= s5["vphi_over_Za_sin"] <= 1.25)
npass = sum(r["ok"] for r in R.values())
summary = dict(exams=R, n_pass=npass, n_fail=len(R) - npass, posthoc_blocks=PH, overlap=OV)
json.dump(summary, open(os.path.join(OUT, "surf_f2_results.json"), "w"), indent=1, default=float)
for k, r in R.items():
    print(k, "PASS" if r["ok"] else "FAIL", r["value"] if k not in ("S3",) else {n: (x["own_frame"], x["sum_err"]) for n, x in r["value"].items()})
print("PASS", npass, "FAIL", len(R) - npass)

# ---------------- figure
fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
cols = {1: "k", 4: "tab:orange", 16: "tab:blue", 25: "tab:green"}
for i, n in enumerate((48, 96)):
    A = ax[0, i]
    for M in (1, 4, 16, 25):
        d = V[(n, M)]
        t = np.array(d["t"]) * d["theta"]
        y = (np.array(d["H"]["rho2"]) - d["floor"]) / (d["H"]["rho2"][0] - d["floor"])
        A.semilogy(t, np.maximum(y, 1e-3), color=cols[M], label=f"M = {M}")
        yb = (np.array(d["H"]["born"]) - d["floor"]) / (d["H"]["rho2"][0] - d["floor"])
        if M == 16:
            A.semilogy(t, np.maximum(np.abs(yb), 1e-3), color="gray", ls=":", label="partenza |psi|^2 (M=16)")
    for key, r in PH.items():
        if key.startswith(f"n{n}_"):
            t = np.array(r["t"]) * V[(n, 1)]["theta"]
            y = (np.array(r["H"]) - r["floor"]) / (r["H"][0] - r["floor"])
            Mk = int(key.split("_M")[1].split("_")[0])
            A.semilogy(t, np.maximum(y, 1e-3), color=cols[Mk], ls="--", lw=0.9,
                       label=f"M={Mk}, blocchi {r['block_phys']:.2g} (post-hoc)")
    A.set_ylim(1e-3, 1.5); A.set_xlabel("tempo fisico (t * a)"); A.set_ylabel("(H - F) / (H0 - F)")
    A.set_title(f"Valentini sul reticolo, flusso esatto, n = {n} (partenza |psi|^4)"); A.legend(fontsize=7)
A = ax[0, 2]
lab = []; xs = []; i = 0
for M in (16, 25):
    for n in (48, 96):
        A.bar(i, f(n, M, "rho2")["tau_sigma"], color="tab:blue"); lab.append(f"reticolo\nn={n} M={M}"); i += 1
    A.bar(i, cont[str(M)]["tau_sigma"], color="tab:red"); lab.append(f"continuo\nM={M}"); i += 1
A.set_xticks(range(i)); A.set_xticklabels(lab, fontsize=7); A.set_ylabel("tau * sigma_omega (adimensionale)")
A.set_title("V6: il reticolo rilassa ~10x prima del continuo")
names = ["0.84037", "0.88627", "0.89930", "0.89960", "E5_listened", "smooth_nucleus_1s"]
short = ["0.8404", "0.8863", "0.8993", "0.8996", "E5 (vecchio)", "nucleo spalmato"]
vk = ["0", "xy", "xz", "yz"]; vc = ["tab:blue", "tab:orange", "tab:green", "tab:purple"]
A = ax[1, 0]
bot = np.zeros(len(names))
for v, c in zip(vk, vc):
    h = np.array([st[nm]["valleys"][v]["weight"] for nm in names]); A.bar(range(len(names)), h, bottom=bot, color=c, label=f"valle {v}"); bot += h
A.set_xticks(range(len(names))); A.set_xticklabels(short, fontsize=8, rotation=20); A.set_ylabel("peso della valle")
A.set_title("Pesi delle 4 valli (livelli esatti: 25 % ciascuna)"); A.legend(fontsize=7)
A = ax[1, 1]
wd = 0.18
for j, (v, c) in enumerate(zip(vk, vc)):
    h = [st[nm]["valleys"][v]["sigma_z_lab"] for nm in names]
    A.bar(np.arange(len(names)) + (j - 1.5) * wd, h, wd, color=c, label=f"valle {v}")
A.plot(range(len(names)), [st[nm]["sigma_z"] for nm in names], "k*", ms=12, label="<Sigma_z> totale")
A.axhline(0, color="k", lw=0.5); A.set_xticks(range(len(names))); A.set_xticklabels(short, fontsize=8, rotation=20)
A.set_ylabel("<Sigma_z> (lab) per valle"); A.set_title("Spin per valle; E5 = 49 % 0.8993 + 51 % 0.8996"); A.legend(fontsize=7)
A = ax[1, 2]
x = np.arange(len(names))
A.bar(x - 0.2, [st[nm]["flux_vphi_over_Za_sin"] for nm in names], 0.4, color="k", label="corrente totale")
A.bar(x + 0.2, [st[nm]["valley0_only_flux_vphi"] for nm in names], 0.4, color="tab:blue", label="solo valle k = 0")
A.axhspan(0.75, 1.25, color="tab:green", alpha=0.15, label="finestra E5 (Dirac)")
A.set_xticks(x); A.set_xticklabels(short, fontsize=8, rotation=20); A.set_ylabel("v_phi / (Z alpha sin)")
A.set_title("Rotazione: la valle 0 gira giusta, le altre la cancellano"); A.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "surf_f2_figure.png"), dpi=110)
print("figure written")
