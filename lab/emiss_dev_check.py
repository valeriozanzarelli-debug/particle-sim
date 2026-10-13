"""DEVELOPMENT: compare the decay rate from the norm of the vacuum sector with the rate from the photons (energy x 2/w)"""
import sys
import numpy as np
from zz_core import ALPHA
from emiss_core import Light3D
for nm in sys.argv[1:]:
    d = np.load(f'out/emiss_ck/run_{nm}.npz')
    t = d['ser_t']; nrm = d['ser_norm']; ph = d['ser_eph'] + d['ser_abs']
    w = 0.978799 - 0.899452
    m = t > 200
    sl = -np.polyfit(t[m], nrm[m], 1)[0]
    slp = np.polyfit(t[m], 2 * ph[m] / w, 1)[0]
    D = d['ser_Dc'][-1]
    d2 = np.sum(np.abs(D)**2)
    W2 = Light3D(128, 4.0, sigma=1.0).form_factor2(w)
    A = 4 / 3 * ALPHA * w**3 * d2
    print(f"{nm}: rate(norm) {sl:.4e}  rate(photons) {slp:.4e}  ratio {sl / slp:.4f};  A_rule*W2 {A * W2:.4e}  norm/AW2 {sl / (A * W2):.4f}  photons/AW2 {slp / (A * W2):.4f}")
