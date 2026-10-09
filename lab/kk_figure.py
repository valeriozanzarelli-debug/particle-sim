"""kk_figure.py - lab/out/kk_force.png from the outputs of kk_cyclotron.py and kk_ab.py"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = dict(np.load("out/kk_cyclotron_tmp.npz"))
d.update(dict(np.load("out/kk_ab_tmp.npz")))
J = json.load(open("out/kk_ab.json"))
fig, ax = plt.subplots(1, 3, figsize=(13, 4.2), dpi=80)
img = d["acc_p"] + d["acc_m"]
ax[0].imshow(np.sqrt(img.T), origin="lower", cmap="Greys")
for key, col, lab in (("tr_p", "tab:red", "n = +1 (density peak)"), ("tr_m", "tab:blue", "n = -1 (density peak)"), ("tr_0", "tab:green", "n = 0")):
    tr = d[key]
    ax[0].plot(tr[:, 0], tr[:, 1], color=col, lw=1.5, label=lab)
ax[0].plot([64], [64], "k*", ms=9)
ax[0].set_xlim(0, 127); ax[0].set_ylim(0, 127)
ax[0].set_title("hidden circle twisted uniformly:\nopposite n curl opposite ways, n = 0 straight", fontsize=10)
ax[0].legend(fontsize=8, loc="lower right"); ax[0].set_xlabel("x (cells)"); ax[0].set_ylabel("y (cells)")
ab = J["aharonov_bohm"]
for n, col in ((1, "tab:red"), (-1, "tab:blue"), (0, "tab:green"), (2, "tab:purple")):
    rows = [r for r in ab if r["n"] == n]
    Phi = np.array([r["Phi"] for r in rows]); ph = np.array([r["phase_shift_interference_term"] for r in rows])
    ax[1].plot(Phi, ph, "o", color=col, label=f"n = {n:+d}")
    ff = np.linspace(0, 2 * np.pi, 200)
    ax[1].plot(ff, np.angle(np.exp(1j * n * ff)), ".", color=col, ms=1.2, alpha=0.6)
ax[1].set_xlabel("twist around the hole  Phi (rad)"); ax[1].set_ylabel("phase of the interference term (rad)")
ax[1].set_title("Aharonov-Bohm ring: phase = n Phi (mod 2 pi)\npoints: walk, lines: prediction", fontsize=10)
ax[1].legend(fontsize=8)
x = np.arange(len(d["r1_0"])); w = np.abs(x - 200) < 40
for key, col, lab in (("r1_0", "k", "n=+1, Phi=0"), ("r1_h", "tab:red", "n=+1, Phi=pi/2"),
                      ("rm_h", "tab:blue", "n=-1, Phi=pi/2"), ("r0_h", "tab:green", "n=0, Phi=pi/2")):
    ax[2].plot(x[w], d[key][w], color=col, lw=1.2, label=lab, ls="--" if key == "r0_h" else "-")
ax[2].set_xlabel("x on the ring (cells), meeting point = 200"); ax[2].set_ylabel("density")
ax[2].set_title("fringes where the two halves meet\n(twist on the link at x = 100, or spread: identical)", fontsize=10)
ax[2].legend(fontsize=7)
fig.tight_layout()
fig.savefig("out/kk_force.png")
print("saved")
