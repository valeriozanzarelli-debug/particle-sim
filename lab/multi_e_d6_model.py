"""multi_e_d6_model.py - DIAGNOSI POST HOC di D6: quanto sbaglia il MODELLO del fit h0(1 + c r + d r^2) su r < 0.25 a0
su una funzione nota. Densita' di coppia esatta del campo medio (prodotto di orbitali, cuspide c = 0 per costruzione):
  h(u) = (2 pi / u) int r rho(r) int_{|r-u|}^{r+u} s rho(s) ds dr,  rho = |phi|^2 normalizzata
mediata sui gusci da 0.01 a0 e fittata come i dati (pesata). Il c che esce e' il bias del modello."""
import numpy as np
from multi_e_core import OUT, load_json, save_json
from multi_e_d6_diag import fitw, hb, nb, vol
import os
z = np.load(os.path.join(OUT, "multi_e_orb_He.npz"), allow_pickle=True)
r = z["r"]; f = z["P"][0] / r
f /= np.sqrt(np.trapezoid(4 * np.pi * r ** 2 * f ** 2, r))
rho = f ** 2
rg = np.linspace(0, 12, 120001); rh = np.interp(rg, r, rho)
cum = np.concatenate([[0], np.cumsum(0.5 * (rg[1:] * rh[1:] + rg[:-1] * rh[:-1]) * np.diff(rg))])   # int_0^x s rho(s) ds
def h(u):
    inner = np.interp(rg + u, rg, cum, right=cum[-1]) - np.interp(np.abs(rg - u), rg, cum)
    return 2 * np.pi / u * np.trapezoid(rg * rh * inner, rg)
I = np.zeros(nb)
for i in range(nb):
    us = np.linspace(i * hb, (i + 1) * hb, 21)[1:-1] if i else np.linspace(1e-4, hb, 21)[1:-1]
    us = np.linspace(i * hb + hb / 40, (i + 1) * hb - hb / 40, 20)
    I[i] = np.mean([4 * np.pi * u * u * h(u) for u in us])      # I(r) per unita' di r
c, h0 = fitw(I)
print("h(0) diretto", np.trapezoid(4 * np.pi * rg ** 2 * rh ** 2, rg), "fit h0", h0, "c del modello (vero 0):", c)
res = load_json("multi_e_dmc.json"); d = res.get("D6_diag_posthoc", {})
d["model_bias_on_HF_pair_density"] = dict(c_fit=float(c), h0_fit=float(h0), note="funzione con c = 0 esatto: il c del fit e' il bias del modello")
res["D6_diag_posthoc"] = d; save_json("multi_e_dmc.json", res)
