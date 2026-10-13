"""emiss_analysis.py - measured rates from the qe runs (shared by the dev checks and the final judge)."""
import numpy as np
from zz_core import ALPHA
from emiss_core import Light3D, CK
import os


def load(name):
    return np.load(os.path.join(CK, f"run_{name}.npz"))


def qe_rates(name, omega, b, t_from_periods=3.0, nL=128, sigma=1.0, transverse=True):
    d = load(name)
    t = d["ser_t"]; nrm = d["ser_norm"]; eph = d["ser_eph"] + d["ser_abs"]
    D = d["ser_Dc"]
    D0 = D[0]
    d2 = float(np.sum(np.abs(D0)**2))
    A_rule = 4 / 3 * ALPHA * omega**3 * d2
    L = Light3D(nL, b, sponge=0, sigma=sigma)
    W2 = L.form_factor2(omega)
    t0 = t_from_periods * 2 * np.pi / omega
    m = t >= t0
    # rate from the norm of the vacuum sector: -d ln(norm)/dt (linear fit of ln norm)
    A_norm = -np.polyfit(t[m], np.log(nrm[m]), 1)[0]
    # rate from the photons: N = 2 (energy in box + absorbed) / omega; dN/dt / norm
    N = 2 * eph / omega
    A_ph = np.polyfit(t[m], N[m], 1)[0] / np.mean(nrm[m])
    return dict(A_norm=A_norm, A_photons=A_ph, A_rule=A_rule, W2=W2, A_rule_W2=A_rule * W2,
                ratio_norm=A_norm / (A_rule * W2), ratio_photons=A_ph / (A_rule * W2), d2=d2, t_end=float(t[-1]),
                depletion=float(1 - nrm[-1]), A_norm_corr=A_norm / W2)
