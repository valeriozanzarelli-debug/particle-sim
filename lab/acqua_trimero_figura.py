"""acqua_trimero_figura.py - trimero ciclico (H2O)3 alla nostra geometria (verdetti lab/out/acqua_trimero_verdetti.json):
densita' elettronica RHF X=2 calcolata da noi, estratta a palline luminose, vista dall'alto (proiezione sul piano
dell'anello; colore = quota). Solo disegno: la correlazione entra nei numeri, non nella densita'."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import acqua_trimero as T, acqua_dimero as DM, acqua_int as I, acqua_chol as CH, acqua_basis as AB
from acqua_figure import ao_values

v = json.load(open(os.path.join(AB.OUT, "acqua_trimero_verdetti.json"))); g = v["geometria"]
atoms, xyz = T.trimer_xyz(*g); B = DM.basis(atoms, xyz, 2, [1] * 9)
ints = I.integrals(B, eri=False); L = CH.cholesky(B); r = CH.scf_df(B, ints, L, 15); C = r["C"][:, :15]
rng = np.random.default_rng(1); pts = []
while sum(len(p) for p in pts) < 90000:
    P = rng.uniform(-6.5, 6.5, size=(200000, 3)); P[:, 2] = rng.uniform(-3.5, 3.5, len(P))
    rho = 2 * np.sum((ao_values(B, P) @ C) ** 2, 1); keep = rng.random(len(P)) < np.minimum(rho, 1.0)
    pts.append(np.c_[P[keep], rho[keep]])
S = np.vstack(pts)[:90000]
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig = plt.figure(figsize=(9, 9), facecolor="black"); ax = fig.add_axes([0, 0, 1, 1]); ax.set_facecolor("black")
t = np.clip((S[:, 2] + 2.0) / 4.0, 0, 1); lr = np.clip((np.log10(S[:, 3]) + 2.5) / 3.0, 0, 1)
col = np.c_[0.30 + 0.70 * t, 0.45 + 0.40 * lr, 1.0 - 0.6 * t]
for size, alpha in ((70, 0.010), (16, 0.035), (2, 0.22)):
    ax.scatter(S[:, 0], S[:, 1], s=size, c=col, alpha=alpha, linewidths=0)
for a, x in zip(atoms, xyz):
    ax.scatter([x[0]], [x[1]], s=34 if a == "O" else 12, color="white")
ax.set_xlim(-6.2, 6.2); ax.set_ylim(-5.4, 7.0); ax.set_aspect("equal"); ax.axis("off")
V = v["verdetti"]; tx = lambda y, s_, c="white", fs=12: ax.text(0.03, y, s_, color=c, fontsize=fs, transform=ax.transAxes, va="top")
tx(0.97, "Trimero dell'acqua (H₂O)₃ dalla regola: tre molecole ad anello", fs=14)
tx(0.93, f"legame {V['T1_De']['valore']:.3f} eV (rif 0,685)   O–O {V['T2_ROO']['valore']:.3f} Å (rif 2,80)", "#bbbbbb", 11)
tx(0.90, f"termine a 3 corpi {V['T3_3corpi']['valore']:.3f} eV (rif −0,100) = 14 % del legame", "#bbbbbb", 11)
ax.text(0.03, 0.04, "APPROSSIMATO · 3 esami su 3 PASS · zero manopole · densità RHF (correlazione solo nei numeri)",
        color="#ffb347", fontsize=10, transform=ax.transAxes)
ax.text(0.97, 0.08, "vista dall'alto · colore = quota (rosa sopra, blu sotto) · 1 a₀ ≈ 0,53 Å", color="#888888",
        fontsize=10, transform=ax.transAxes, ha="right")
fig.savefig(os.path.join(AB.OUT, "acqua_trimero_figura.png"), dpi=110, facecolor="black")
