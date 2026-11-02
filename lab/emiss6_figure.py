"""emiss6_figure.py - figure of phase 4 (docs/EMISSIONE.md §13) -> lab/out/emiss6_figure.png"""
import os, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from emiss6_core import CK, OUT, VEC
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
d = np.load(os.path.join(CK, "emiss6_q60_0.025.npz"))
r = json.load(open(os.path.join(OUT, "emiss6_res_q60_0.025.json")))
t, n = d["ser_t"], d["ser_norm"]
ax[0].semilogy(t * r["A_fit"], n, color="#c040a0", lw=2, label="norma (atomo ancora eccitato)")
ax[0].semilogy(t * r["A_fit"], np.exp(-r["A_fit"] * t), "k:", label="exp(-A t)")
ax[0].set_xlabel("t · A"); ax[0].set_title("(a) Z=60, 2p puro, K=1: decadimento"); ax[0].legend()
J = d["J"] @ np.conj(VEC); nf = 8 * len(J)
S = np.abs(np.fft.fft(J, nf))**2; fr = -2 * np.pi * np.fft.fftfreq(nf, d=r["b"])
sel = np.abs(fr - r["line_centre"]) < 6 * r["fwhm"]
ax[1].plot((fr[sel] - r["omega"]) / r["A_fit"], S[sel] / S[sel].max(), ".", color="#3070d0", ms=3, label="spettro della luce emessa")
x = np.linspace(-6, 6, 400) * r["fwhm"] / r["A_fit"]
c = (r["line_centre"] - r["omega"]) / r["A_fit"]
ax[1].plot(x, 1 / (1 + ((x - c) / (r["fwhm"] / r["A_fit"] / 2))**2), "k-", lw=1, label=f"Lorentz, FWHM/A = {r['fwhm_over_A']:.4f}")
ax[1].set_xlabel("(ω − ω_regola)/A"); ax[1].set_title("(b) riga naturale"); ax[1].legend()
dr = np.load(os.path.join(CK, "emiss6_rabi60_0.025.npz")); m = int(dr["irec"][0])
ax[2].plot(dr["out_t"][:m], dr["out_pe"][:m], color="#e08020", label="P(2p3/2)")
ax[2].plot(dr["out_t"][:m], dr["out_pg"][:m], color="#30a060", label="P(1s)")
ax[2].set_xlabel("t (ħ/mc²)"); ax[2].set_title("(c) assorbimento: Rabi, Z=60"); ax[2].legend()
plt.tight_layout(); plt.savefig(os.path.join(OUT, "emiss6_figure.png"), dpi=110)
