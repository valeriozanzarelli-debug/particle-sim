"""h2v_figure.py - figure of VIDEO H2 -> lab/out/h2v_figure.png (and docs/fig/h2v_figure.png)"""
import os, sys, json, shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h2v_sim import SCR

LAB = os.path.dirname(os.path.abspath(__file__))
es = json.load(open(os.path.join(LAB, "out", "h2v_esami.json")))
fig = plt.figure(figsize=(16, 9.6))
for k, (fn, tit) in enumerate([("h2v_foto_A.png", "A · la molecola: 24 000 palline da |Ψ(r₁,r₂)|²"),
                               ("h2v_foto_B.png", "B · dove sta l'altro elettrone (luce = g(r₂|r_rosa))"),
                               ("h2v_foto_C.png", "C · il legame vibra (v0+v1, periodo 8,0 fs)")]):
    ax = fig.add_subplot(2, 3, k + 1)
    p = os.path.join(SCR, fn)
    if os.path.exists(p):
        ax.imshow(plt.imread(p))
    ax.set_title(tit, fontsize=10); ax.axis("off")
ax = fig.add_subplot(2, 3, 4)
H = np.load(os.path.join(SCR, "h2v_V3_hist.npy"))
zc = np.linspace(-4, 4, 41); zc = 0.5 * (zc[1:] + zc[:-1])
ax.plot(zc, H[1].sum(0) / 0.2, "-", color="k", lw=1.5, label="calcolato (gaussiane esatte)")
ax.plot(zc, H[0].sum(0) / 0.2, "o", ms=3.5, color="tab:blue", label=f"palline ({es['V3']['n_balls']:,} accumulate)".replace(",", " "))
ax.axvline(-0.7, color="tab:red", lw=0.6); ax.axvline(0.7, color="tab:red", lw=0.6)
ax.set_xlabel("z lungo il legame (a0)"); ax.set_ylabel("elettroni per a0")
ax.set_title(f"V3: Pearson {es['V3']['pearson']:.5f}, L1 {es['V3']['L1']:.3f}", fontsize=10); ax.legend(fontsize=8)
ax = fig.add_subplot(2, 3, 5)
p6 = os.path.join(SCR, "h2v_V6.npy")
if os.path.exists(p6) and "V6" in es:
    t, R, Rm = np.load(p6)
    ax.plot(t, Rm, "-", color="k", lw=1, label="<R>(t) calcolato")
    ax.plot(t, R, "o", ms=3, color="tab:orange", label="letto dalle immagini (bagliori)")
    ax.set_xlabel("tempo nel video (s)  [4 s = 8,0 fs]"); ax.set_ylabel("R (a0)")
    ax.set_title(f"V6: periodo {es['V6']['T_fs_frames']:.3f} fs contro {es['V6']['T_fs_from_nu01']:.3f}", fontsize=10)
    ax.legend(fontsize=8)
ax = fig.add_subplot(2, 3, 6)
rows = []
for k in ["V0", "V1", "V2", "V3", "V4", "V4b", "V5", "V6"]:
    if k in es:
        rows.append((k, es[k]["PASS"]))
if "V7_v1" in es:
    rows.append(("V7 (pre-reg.)", es["V7_v1"]["PASS"]))
if "V7" in es:
    rows.append(("V7 rifatto (dopo)", es["V7"]["PASS"]))
ax.barh(range(len(rows)), [1] * len(rows), color=["tab:green" if p else "tab:red" for _, p in rows])
ax.set_yticks(range(len(rows))); ax.set_yticklabels([f"{k}  {'PASS' if p else 'FAIL'}" for k, p in rows])
ax.invert_yaxis(); ax.set_xticks([]); ax.set_title("esami pre-registrati", fontsize=10)
fig.suptitle("Video H₂ — la prima molecola: due protoni (fermi), due elettroni dallo stato ECG di Debiti", fontsize=13)
fig.tight_layout()
out = os.path.join(LAB, "out", "h2v_figure.png")
fig.savefig(out, dpi=110)
os.makedirs(os.path.join(LAB, "..", "docs", "fig"), exist_ok=True)
shutil.copy(out, os.path.join(LAB, "..", "docs", "fig", "h2v_figure.png"))
print("ok")
