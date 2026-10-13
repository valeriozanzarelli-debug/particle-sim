"""emiss_x5c.py - declared extra diagnosis for X5: WHERE in space does the weak dipole line (~0.090) live?
Free evolution of the P4_E1pi final state; demodulate the density at the main line and at the weak line
(rho_w(x) = sum_t h(t) e^{i w t} rho(x, t) / sum h) and split the z-dipole of each line by radius.
Run: python lab/emiss_x5c.py -> lab/out/emiss_x5c.json"""
import json, os, time
import numpy as np
from emiss_core import CK, OUT
from emiss_v4_dev import make_atom, S_LEARNED

t0 = time.time()
d = np.load(os.path.join(CK, "run_P4_E1pi.npz"))
psi = d["psi"].copy(); tk = int(d["t_ticks"])
D = make_atom(n=48, s=S_LEARNED)
W = {"main_3/2": 0.07937, "main_1/2": 0.07651, "weak_0.0897": 0.08966, "weak_0.0925": 0.09250}
NT = 6000
acc = {k: np.zeros(psi.shape[1:], complex) for k in W}
hs = 0.0
for t in range(NT):
    if (tk + t) % 2 == 0:
        tt = t * D.a
        h = np.sin(np.pi * t / NT)**2
        rho = np.sum(np.abs(psi)**2, 0)
        for k, w in W.items():
            acc[k] += h * np.exp(1j * w * tt) * rho
        hs += h
    D.tick(psi, tk + t)
r = D.r
Zc = D.x[None, None, :] * np.ones_like(r)
res = {}
for k in W:
    a = acc[k] / hs
    zc = a * Zc
    tot = abs(np.sum(zc))
    bins = [0, 1, 2, 4, 6, 8, 10, 12, 100]
    prof = [float(abs(np.sum(zc[(r >= bins[i]) & (r < bins[i + 1])])) / max(tot, 1e-30)) for i in range(len(bins) - 1)]
    # parity pattern of the demodulated density: weight of the (i+j+k) staggered and 4-cell modulations
    F = np.abs(np.fft.fftn(a))**2
    n = a.shape[0]
    q = np.abs(2 * np.pi * np.fft.fftfreq(n))
    smooth = F[np.ix_(q < np.pi / 4, q < np.pi / 4, q < np.pi / 4)].sum() / F.sum()
    res[k] = dict(omega=W[k], z_amp=float(tot), radial_share=dict(zip([f"{bins[i]}-{bins[i + 1]}" for i in range(len(bins) - 1)], prof)),
                  smooth_k_fraction=float(smooth))
    print(f"{k} (w {W[k]}): |z| {tot:.4e}; share of |z| by radius {[round(x, 3) for x in prof]}; smooth-k fraction of the oscillating density {smooth:.3f}", flush=True)
json.dump(res, open(os.path.join(OUT, "emiss_x5c.json"), "w"), indent=1)
print(f"({time.time() - t0:.0f} s)")
