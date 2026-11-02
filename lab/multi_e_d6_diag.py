"""multi_e_d6_diag.py - DIAGNOSI POST HOC (dichiarata, non rigiudica D6) del FAIL di D6.
Tre cause da separare:
  (1) buco dei camminatori a r12 -> 0 a tau finito: con la guida u = 0 l'energia locale contiene +1/r12 senza
      compensazione e il peso primitivo e^{-tau E_L} svuota la zona entro la lunghezza di diffusione sqrt(2 D tau)
      -> si guarda come cambia con tau (0.01, 0.0025, 0.000625);
  (2) il fit pre-registrato non e' pesato: i primi gusci (volume ~ r^2 dr) hanno pochissimi conteggi e pesano come gli
      altri -> qui fit pesato con l'inverso della varianza di Poisson (peso ~ volume del guscio);
  (3) lo stimatore estrapolato 2 I_mix - I_VMC sbaglia al secondo ordine in (phi - psi_T); al contatto la guida u = 0
      ha h(0) quasi doppio del vero, quindi l'errore non e' piccolo. Identita': se il misto ha cuspide m/2 e la VMC 0,
      c_est = m h_mix(0) / (2 h_mix(0) - h_VMC(0)) != m. Invece il MISTO ha cuspide esatta m/2 (psi_T phi, psi_T senza
      cuspide), e quella si puo' controllare senza estrapolare.
Run: python lab/multi_e_d6_diag.py 0.01 0.0025 0.000625"""
import sys
import numpy as np
from multi_e_d6 import blocks, density
from multi_e_core import load_json, save_json

hb = 0.01; nb = 25
e = np.arange(nb + 1) * hb; lo, hi = e[:-1], e[1:]
vol = 4 * np.pi * (hi ** 3 - lo ** 3) / 3
mom = lambda n: 3 * (hi ** (n + 3) - lo ** (n + 3)) / ((n + 3) * (hi ** 3 - lo ** 3))
X = np.stack([mom(0), mom(1), mom(2)], 1)


def fitw(I):
    h = I[:nb] * hb / vol; w = np.sqrt(vol)
    a, *_ = np.linalg.lstsq(X * w[:, None], h * w, rcond=None)
    return a[1] / a[0], a[0]


def jk(fun, Hm, Wm, Hv, Wv, NG=20):
    gm = np.array_split(np.arange(len(Wm)), NG); gv = np.array_split(np.arange(len(Wv)), NG)
    v = np.array([fun(np.delete(Hm, gm[i], 0), np.delete(Wm, gm[i]), np.delete(Hv, gv[i], 0), np.delete(Wv, gv[i])) for i in range(NG)])
    return float(np.sqrt((NG - 1) / NG * np.sum((v - v.mean()) ** 2)))


res = load_json("multi_e_dmc.json")
m = res["He_u0_dmc_tau0.01"]["m_kin"]
Hv, Wv = blocks("vmc")
out = res.get("D6_diag_posthoc", {})
for tau in [float(t) for t in sys.argv[1:]]:
    Hm, Wm = blocks("dmc", tau)
    Im, Iv = density(Hm, Wm, hb), density(Hv, Wv, hb)
    c_mix, h_mix = fitw(Im); c_v, h_v = fitw(Iv); c_e, h_e = fitw(2 * Im - Iv)
    rec = dict(tau=tau, diffusion_length=float(np.sqrt(2 * res["He_u0_dmc_tau0.01"]["D"] * tau)),
               c_mix_weighted=float(c_mix), sigma_c_mix=jk(lambda a, b, c, d: fitw(density(a, b, hb))[0], Hm, Wm, Hv, Wv),
               c_est_weighted=float(c_e), sigma_c_est=jk(lambda a, b, c, d: fitw(2 * density(a, b, hb) - density(c, d, hb))[0], Hm, Wm, Hv, Wv),
               c_vmc_weighted=float(c_v), h0_mix=float(h_mix), h0_vmc=float(h_v), h0_est=float(h_e),
               c_est_predicted_by_identity=float((2 * h_mix * m / 2 - h_v * c_v) / (2 * h_mix - h_v)),
               first_bin_mix_over_plateau=float((Im[0] * hb / vol[0]) / np.mean(Im[5:15] * hb / vol[5:15])),
               target_mix=m / 2, target_est=m)
    out[f"tau{tau:g}"] = rec
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rec.items()}, flush=True)
res["D6_diag_posthoc"] = out
save_json("multi_e_dmc.json", res)
