"""rifl_figura.py - figura dello studio RIFLUSSO (lab/out/rifl_figura.png).
(a) lo stato di massimo riflusso della regola che si muove (densita' R arancio, L blu) e la probabilita' a sinistra;
(b) c_rbf(eps) di Dirac: nostro continuo, tabella di Ashfaque et al., la regola per a -> 0, il limite 1/(pi eps^2);
(c) la regola a theta finito: Lambda(theta, N) contro eps = sqrt(4/(theta N)), finestra fisica e banda intera;
(d) limite NR della regola: lambda(U) contro 1/U, verso c_bm."""
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rifl_core import *
from zz_core import Walk

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
r1 = json.load(open(os.path.join(OUT, "rifl_r1.json")))
r2 = json.load(open(os.path.join(OUT, "rifl_r2.json")))
r3 = json.load(open(os.path.join(OUT, "rifl_r3.json")))

fig = plt.figure(figsize=(13, 10))

# (a) stato di massimo riflusso sulla regola (anello grande, k sulla griglia dell'anello)
theta, N, Lr = 0.1, 200, 1600
j = np.arange(1, Lr // 4)
k = 2 * np.pi * j / Lr
lam, v = top_eig(walk_kernel(k, np.full(len(k), 2 * np.pi / Lr), theta, N), True)
x = np.arange(Lr) - Lr // 2
_, uR, uL = walk_modes(k, theta)
ph = np.exp(1j * np.outer(x, k))
psi = np.array([ph @ (v * uR), ph @ (v * uL)])
psi /= np.sqrt(np.sum(np.abs(psi)**2))
W = Walk("dirac1", Lr, 1.0, dim=1, theta=theta)
i0 = int(np.where(x == 0)[0][0])
win = slice(i0 - 300, i0 + 301)
nshow = 2 * N
RR, LL, back = [], [], [0.0]
for n in range(nshow):
    RR.append(np.abs(psi[0, win])**2)
    LL.append(np.abs(psi[1, win])**2)
    if n < nshow:
        back.append(back[-1] + abs(psi[1, i0 + 1])**2 - abs(psi[0, i0])**2)
    psi = W.tick(psi)
RR, LL = np.array(RR), np.array(LL)
sc = max(RR.max(), LL.max())
img = np.ones(RR.shape + (3,))
r, l = np.clip(RR / RR.max(), 0, 1) ** 0.6, np.clip(LL / LL.max(), 0, 1) ** 0.6   # L in scala propria
img[..., 0] -= 0.1 * r + 0.85 * l
img[..., 1] -= 0.55 * r + 0.55 * l
img[..., 2] -= 0.95 * r + 0.1 * l
ax = fig.add_subplot(2, 2, 1)
ax.imshow(np.clip(img, 0, 1), origin="lower", aspect="auto", extent=[-300, 300, 0, nshow])
ax.axvline(0.5, color="k", lw=0.6, ls=":")
ax.axhline(N, color="k", lw=0.6, ls="--")
ax.set_xlabel("sito x (taglio tra 0 e 1)")
ax.set_ylabel("tick n")
ax.set_title(f"(a) massimo riflusso sulla regola, θ={theta}, N={N}: |R|² arancio, |L|² blu (scala propria)", fontsize=10)
ax2 = ax.inset_axes([0.62, 0.08, 0.35, 0.3])
ax2.plot(np.arange(len(back)), back, color="k", lw=1)
ax2.axvline(N, color="gray", ls="--", lw=0.6)
ax2.axhline(lam, color="tab:red", lw=0.6)
ax2.set_title("P_sin(n) − P_sin(0)", fontsize=8)
ax2.tick_params(labelsize=7)

# (b) Dirac c_rbf(eps)
ax = fig.add_subplot(2, 2, 2)
q1 = r3["Q1_eps"]
e = np.array(sorted(float(z) for z in q1))
c = np.array([q1[f"{z:g}"]["c_rbf"] for z in e])
e1 = np.array(sorted(float(z) for z in r1["E1"]["per_eps"]))
c1 = np.array([r1["E1"]["per_eps"][f"{z:.2f}"]["c_rbf"] for z in e1])
ee = np.concatenate([e1, e[e > e1.max()]])
cc = np.concatenate([c1, c[e > e1.max()]])
ax.loglog(ee, cc, "o-", color="tab:blue", ms=4, label="Dirac continuo (nostro)")
tab = r1["E1"]["per_eps"]
ax.loglog(e1, [tab[f"{z:.2f}"]["tabella"] for z in e1], "x", color="k", ms=7, label="Ashfaque et al. 2019")
e3 = r2.get("E3", {}).get("per_eps", {})
if e3:
    ax.loglog([float(z) for z in e3], [e3[z]["regola_a0"] for z in e3], "s", mfc="none", color="tab:red", ms=8,
              label="regola zig-zag, a → 0")
eg = np.logspace(-1, 4, 200)
ax.loglog(eg, 1 / (np.pi * eg**2), "--", color="tab:green", label="limite 1/(π ε²) = mc²T/(4πħ)")
ax.axhline(C_BM, color="gray", ls=":", label="c_bm = 0,0384517")
ax.set_ylim(1e-9, 0.06)
ax.set_xlabel("ε = √(4ħ / mc²T)")
ax.set_ylabel("massimo riflusso")
ax.set_title("(b) Dirac: dal non relativistico all'ultrarelativistico")
ax.legend(fontsize=8)

# (c) la regola a theta finito
ax = fig.add_subplot(2, 2, 3)
q2 = r3["Q2_theta"]
cols = plt.cm.viridis(np.linspace(0, 0.9, len(q2)))
for col, (th, rows) in zip(cols, sorted(q2.items(), key=lambda z: float(z[0]))):
    Ns = sorted(int(n) for n in rows)
    epsN = [rows[str(n)]["eps"] for n in Ns]
    ax.semilogx(epsN, [rows[str(n)]["fisica"] for n in Ns], "-", color=col, lw=1.2, label=f"θ={th}")
    ax.semilogx(epsN, [rows[str(n)]["intera"] for n in Ns], ":", color=col, lw=1.0)
ax.semilogx(ee, cc, "k-", lw=2, alpha=0.4, label="Dirac continuo")
ax.axhline(C_BM, color="gray", ls=":")
ax.set_xlim(0.03, 40)
ax.set_xlabel("ε = √(4/(θN))")
ax.set_ylabel("massimo riflusso in N tick")
ax.set_title("(c) la regola a θ finito: finestra fisica (—), banda intera (···)")
ax.legend(fontsize=7, ncol=2)

# (d) limite NR
ax = fig.add_subplot(2, 2, 4)
E0, E2 = r1["E0"], r2["E2"]
ax.plot(1 / np.array(E0["U"]), E0["lambda"], "o-", color="tab:blue", label="NR continuo (Bracken–Melloy)")
ax.plot(1 / np.array(E2["U"]), E2["regola_a0"], "s", mfc="none", color="tab:red", ms=8, label="regola zig-zag, a → 0")
ax.plot([0], [C_BM], "k*", ms=12, label="c_bm = 0,0384517")
ax.plot([0], [E2["c_regola"]], "rx", ms=10, label=f"regola estrapolata {E2['c_regola']:.5f}")
ax.set_xlabel("1 / U   (U = p_max √(T/2m))")
ax.set_ylabel("massimo riflusso")
ax.set_title("(d) limite non relativistico: la regola ritrova c_bm")
ax.legend(fontsize=8)

fig.tight_layout()
fn = os.path.join(OUT, "rifl_figura.png")
fig.savefig(fn, dpi=110)
print("scritto", fn, "lambda(a) =", lam)
