"""search_monopole.py - promozione a 3+1 D della regola vincente "luce + struttura con vuoto al polo"
(docs/RICERCA_REGOLE.md, parte C) e numeri del vincolo C1.

In 3+1 D l'analogo del modello sigma O(3) con gauge (vuoto al polo: il fotone resta senza massa) e' il modello
di Georgi-Glashow SU(2) + Higgs aggiunto; il suo solitone topologico e' il monopolo di 't Hooft-Polyakov
('t Hooft, Nucl. Phys. B 79, 276 (1974); Polyakov, JETP Lett. 20, 194 (1974)).  Ansatz a riccio:
  Phi^a = v h(r) r^a,  A_i^a = eps_aij r^j (1 - K(r)) / (g r),  xi = g v r
  E = (4 pi v / g) f(beta),  f = int dxi [ K'^2 + (K^2-1)^2/(2 xi^2) + xi^2 h'^2/2 + K^2 h^2 + (beta/4) xi^2 (h^2-1)^2 ],
  beta = lambda/g^2  (potenziale lambda/4 (Phi^2 - v^2)^2; m_W = g v, m_H = sqrt(2 lambda) v).
  Limite BPS beta = 0: K = xi/sinh xi, h = coth xi - 1/xi, f = 1.  Per beta -> infinito f -> 1.787 (Kirkman-Zachos 1981).
Uso: python lab/search_monopole.py
"""
import json
import os

import numpy as np
from scipy.optimize import minimize

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
ALPHA = 1 / 137.035999
M_E, M_P, M_MU, M_TAU = 0.51099895, 938.27208816, 105.6583755, 1776.86


def mono_f(beta, xmax=30.0, n=900):
    """minimizza il funzionale radiale (differenze finite, L-BFGS di torch con gradiente esatto)."""
    import torch
    torch.set_default_dtype(torch.float64)
    torch.set_num_threads(1)
    xg = torch.linspace(0, xmax, n + 1)
    dx = float(xg[1] - xg[0])
    xi = xg[1:-1]
    K0 = xi / torch.sinh(xi)
    h0 = 1 / torch.tanh(xi) - 1 / xi
    zK = K0.clone().requires_grad_(True)
    zh = h0.clone().requires_grad_(True)
    one, zero = torch.ones(1), torch.zeros(1)

    def dens(zK, zh):
        K = torch.cat([one, zK, zero])
        h = torch.cat([zero, zh, one])
        Km, hm = 0.5 * (K[1:] + K[:-1]), 0.5 * (h[1:] + h[:-1])
        xm = 0.5 * (xg[1:] + xg[:-1])
        dK, dh = (K[1:] - K[:-1]) / dx, (h[1:] - h[:-1]) / dx
        return (dK ** 2 + (Km ** 2 - 1) ** 2 / (2 * xm ** 2) + xm ** 2 * dh ** 2 / 2 + Km ** 2 * hm ** 2
                + beta / 4 * xm ** 2 * (hm ** 2 - 1) ** 2), xm
    opt = torch.optim.LBFGS([zK, zh], lr=1, max_iter=3000, tolerance_grad=1e-10, tolerance_change=1e-14,
                            history_size=50, line_search_fn="strong_wolfe")

    def closure():
        opt.zero_grad()
        e, _ = dens(zK, zh)
        E = e.sum() * dx
        E.backward()
        return E
    opt.step(closure)
    with torch.no_grad():
        e, xm = dens(zK, zh)
        E = float(e.sum() * dx)
        cum = torch.cumsum(e, 0) * dx
        r_half = float(np.interp(0.5 * float(cum[-1]), cum.numpy(), xm.numpy()))
    # coda coulombiana magnetica 1/(2 xi^2) oltre il taglio; per beta=0 la soluzione BPS e la condizione h(xmax)=1
    # si compensano (valore esatto f=1), quindi nessuna correzione
    tail = 0.0 if beta == 0 else 1 / (2 * xmax)
    return dict(beta=beta, mH_over_mW=float(np.sqrt(2 * beta)), f=E, f_cut_correction=tail,
                r_half_energy_over_mW=r_half)


def main():
    out = {}
    rows = []
    for beta in [0.0, 0.5, 2.0, 10.0, 100.0]:
        r = mono_f(beta)
        rows.append(r)
        print(f"beta = lambda/g^2 = {beta:7.2f}: f = {r['f']:.4f}  (taglio a xi=30: coda 1/(2 xi) = {r['f_cut_correction']:.4f}"
              f" da aggiungere)  m_H/m_W = {r['mH_over_mW']:.3f}  raggio a meta' energia = {r['r_half_energy_over_mW']:.3f}/m_W", flush=True)
    for r in rows:
        r["f_corr"] = r["f"] + r["f_cut_correction"]
    out["monopolo"] = rows
    # predizioni: se il fotone della regola e' il fotone vero, g = e (il W ha carica e), e^2 = 4 pi alpha
    pred = []
    for r in rows:
        MW_ratio = r["f_corr"] / ALPHA  # M/m_W = 4 pi f / e^2 = f / alpha
        for name, M in [("elettrone", M_E), ("protone", M_P)]:
            mW = M / MW_ratio
            pred.append(dict(beta=r["beta"], scala=name, M_over_mW=MW_ratio, m_W_MeV=mW,
                             raggio_meta_energia_fm=r["r_half_energy_over_mW"] * 197.3269804 / mW,
                             carica_magnetica_in_unita_e="g_m = 4 pi / e = 1/(2 alpha) e = 68.5 e"))
    out["predizioni"] = pred
    # look-elsewhere: f(beta) in [1, 1.787] => M/m_W in [137.0, 245]: quanti bersagli cadono nell'intervallo?
    lo, hi = rows[0]["f_corr"] / ALPHA, rows[-1]["f_corr"] / ALPHA
    out["intervallo_M_su_mW"] = [lo, hi]
    out["bersagli_nell_intervallo"] = {"mu/e 206.77": bool(lo < 206.768 < hi)}
    out["nota_look_elsewhere"] = ("con beta libero, M/m_W copre con continuita' [%.1f, %.1f]: 'indovinare' 206.77 "
                                  "scegliendo beta e' un fit, non una predizione" % (lo, hi))
    # C1: PVLAS e scala E*
    A_e = 2 * ALPHA ** 2 / 45  # coeff. (E^2-B^2)^2 in unita' m_e^-4 (Heaviside-Lorentz, hbar=c=1)
    dn_qed_per_T2 = 3.96e-24  # Delta n / B^2 della QED (PVLAS 2020, Ejlli et al.)
    bound_2sigma = (12 + 2 * 17) * 1e-23  # |Delta n|/B^2 < 4.6e-22 T^-2 (2 sigma), PVLAS
    # un termine extra c1 (E^2-B^2)^2 da' Delta n = -4 c1 B^2 ; la QED ha Delta n = (c2 - 4 c1) B^2 = 3 A_e B^2 (c2 = 7 A_e)
    c1_max_over_qed = bound_2sigma / (4 * dn_qed_per_T2 / 3)
    # E*: campo a cui il termine extra vale l'1% di Maxwell: 2 c1 E^2 = 0.01, E in unita' di E_c = m^2/e
    # 2 c1 E^2 = 2 (c1/A_e) A_e (E/E_c)^2 / (4 pi alpha)
    Ec = 1.32e18  # V/m
    ratio_coeff = 2 * c1_max_over_qed * A_e / (4 * np.pi * ALPHA)
    Estar = np.sqrt(0.01 / ratio_coeff) * Ec
    out["C1"] = dict(PVLAS_Delta_n_over_B2="(+12 +- 17) e-23 T^-2 (Ejlli et al., Phys. Rep. 871 (2020) 1)",
                     QED_Delta_n_over_B2=dn_qed_per_T2, c1_extra_max_su_c1_QED_2sigma=c1_max_over_qed,
                     E_star_minimo_V_per_m=float(Estar),
                     laser_piu_intensi_V_per_m=8.7e14,
                     nota="un termine non lineare della luce compatibile con PVLAS e' < 1% di Maxwell sotto E* ~ 1e18 V/m")
    # contributo di un bosone vettore carico di massa m_W al g-2 dell'elettrone (ordine di grandezza)
    out["g-2_stima"] = dict(formula="delta a_e ~ (alpha/pi)^2 (m_e/m_W)^2 / 45 (polarizzazione del vuoto, ordine di grandezza)",
                            m_W_MeV_se_M_protone=pred[1]["m_W_MeV"],
                            delta_a_e=(ALPHA / np.pi) ** 2 * (M_E / pred[1]["m_W_MeV"]) ** 2 / 45,
                            precisione_a_e=1.3e-13)
    with open(os.path.join(OUT, "search_monopole.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "monopolo"}, indent=1))


if __name__ == "__main__":
    main()
