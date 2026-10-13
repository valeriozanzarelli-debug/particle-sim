"""multi_e_d6.py - esame D6 (pre-registrato): cuspide di Kato elettrone-elettrone vista dai camminatori.
He, guida con u = 0 (nessuna cuspide messa a mano nella guida), tau = 0.01, intracula a bin 0.01 a0.
  I_mix  = istogramma pesato dei camminatori DMC (campiona psi_T * phi)
  I_VMC  = istogramma della VMC con la stessa guida (campiona psi_T^2)
  I_est  = 2 I_mix - I_VMC  (stimatore estrapolato, errore del secondo ordine in phi - psi_T)
  h(r)   = I_est / (4 pi r^2): densita' di coppia mediata sulla sfera; fit h0 (1 + c r + d r^2) su r < 0.25 a0,
           con il modello mediato esattamente su ogni guscio (peso r^2), quindi niente errore di 'centro del bin'.
Kato per la densita' di coppia: c = 2 mu q1 q2 = m_kin (mu = m/2 per due elettroni). Esame: c = m_kin +- 0.10.
Errore su c: jackknife su 20 gruppi di pezzi (2 a.u. ciascuno), indipendenti fra DMC e VMC.
Lunghezza delle corse (dichiarata prima di guardare l'istogramma): 400 a.u. di produzione ciascuna.
Run: python lab/multi_e_d6.py"""
import os
import numpy as np
from multi_e_core import OUT, load_json, save_json

JS = "multi_e_dmc.json"
NG = 20
RMAX = 0.25


def blocks(mode, tau=0.01):
    z = np.load(os.path.join(OUT, f"multi_e_dmc_He_u0_{mode}_tau{tau:g}.npz"), allow_pickle=True)
    HB = np.array(z["HB"])
    return HB[:, :-1], HB[:, -1]                     # istogrammi per pezzo, peso totale per pezzo


def density(H, W, hbin):
    """I(r) per unita' di r, normalizzata a una coppia: somma(H) / (somma(W) * hbin)."""
    return H.sum(0) / (W.sum() * hbin)


def fit_c(I, hbin, nb):
    edges = np.arange(nb + 1) * hbin
    lo, hi = edges[:-1], edges[1:]
    vol = 4 * np.pi * (hi ** 3 - lo ** 3) / 3
    h = I[:nb] * hbin / vol                          # densita' di coppia media nel guscio
    # media nel guscio di r^n con peso r^2:  3 (hi^(n+3) - lo^(n+3)) / ((n+3) (hi^3 - lo^3))
    mom = lambda n: 3 * (hi ** (n + 3) - lo ** (n + 3)) / ((n + 3) * (hi ** 3 - lo ** 3))
    X = np.stack([mom(0), mom(1), mom(2)], 1)
    a, *_ = np.linalg.lstsq(X, h, rcond=None)
    return a[1] / a[0], a


if __name__ == "__main__":
    import sys
    tau = float(sys.argv[1]) if len(sys.argv) > 1 else 0.01          # 0.01 = esame; altri tau = diagnosi post hoc
    res = load_json(JS)
    hbin = 0.01
    nb = int(round(RMAX / hbin))
    Hm, Wm = blocks("dmc", tau); Hv, Wv = blocks("vmc")     # la VMC non dipende da tau (Metropolis esatto)
    Im, Iv = density(Hm, Wm, hbin), density(Hv, Wv, hbin)
    Ie = 2 * Im - Iv
    c_est, a_est = fit_c(Ie, hbin, nb)
    c_mix, _ = fit_c(Im, hbin, nb)
    c_vmc, _ = fit_c(Iv, hbin, nb)
    # jackknife: gruppi contigui di pezzi (la correlazione fra pezzi da 2 a.u. e' trascurabile, T_corr ~ 1 a.u.)
    gm = np.array_split(np.arange(len(Wm)), NG); gv = np.array_split(np.arange(len(Wv)), NG)
    jk = []
    for g in range(NG):
        km = np.setdiff1d(np.arange(len(Wm)), gm[g]); kv = np.setdiff1d(np.arange(len(Wv)), gv[g])
        jk.append(fit_c(2 * density(Hm[km], Wm[km], hbin) - density(Hv[kv], Wv[kv], hbin), hbin, nb)[0])
    jk = np.array(jk)
    s_c = float(np.sqrt((NG - 1) / NG * np.sum((jk - jk.mean()) ** 2)))
    m_kin = res[f"He_u0_dmc_tau{tau:g}"]["m_kin"]
    out = dict(c_extrapolated=float(c_est), sigma_c=s_c, c_mixed=float(c_mix), c_vmc=float(c_vmc), target=m_kin, tol=0.10,
               PASS=bool(abs(c_est - m_kin) <= 0.10), h0=float(a_est[0]), d=float(a_est[2] / a_est[0]),
               chunks_dmc=int(len(Wm)), chunks_vmc=int(len(Wv)),
               I_est_first_bins=Ie[:nb].tolist(), I_mix_first_bins=Im[:nb].tolist(), I_vmc_first_bins=Iv[:nb].tolist())
    out["tau"] = tau
    res["D6" if tau == 0.01 else f"D6_diag_posthoc_tau{tau:g}"] = out
    save_json(JS, res)
    print({k: v for k, v in out.items() if "bins" not in k})
