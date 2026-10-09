"""dirac_zbw_scan.py - supplementary (added AFTER dirac_free.py failed Z1d by 0.1 % beyond its tolerance):
how the 2+1 D Zitterbewegung amplitude of the mixed packet depends on the packet width sigma.
Nothing in the rule or in the original tolerance is changed; this only shows where 0.475 comes from.
Expectation for a broad packet: A = (hbar/2mc) <m^3/E^3>-like factor ~ 0.5 (1 - (3/2)<k^2>), <k^2> = 2/(4 sigma^2).
Run: timeout 590 python lab/dirac_zbw_scan.py
"""
import json, os, time
import numpy as np
from scipy.optimize import curve_fit
from dirac_core import Dirac, gaussian_packet

t0 = time.time()
out = []
for sig in (2.0, 4.0, 8.0, 16.0):
    L = max(96.0, 12 * sig)
    D = Dirac(384 if sig < 16 else 512, L, 2)
    psi = gaussian_packet(D, sig, [0, 0], [1, 1])
    pk = D.fft(psi)
    ts = np.arange(0, 4 * np.pi, 0.05)
    mx = np.array([D.mean(D.evolve(pk, t), D.x[0]) for t in ts])
    f = lambda t, c0, c1, A, w, ph: c0 + c1 * t + A * np.sin(w * t + ph)
    p, _ = curve_fit(f, ts, mx, p0=[0, 0, 0.5, 2, 0], maxfev=20000)
    A = abs(p[2])
    pred = 0.5 * (1 - 1.5 * 2 / (4 * sig**2))
    out.append(dict(sigma=sig, A=A, omega=abs(p[3]), A_pred_broad=pred))
    print(f"sigma = {sig:5.1f}: A = {A:.5f} (0.5(1 - 3/(4 sigma^2)) = {pred:.5f}), omega = {abs(p[3]):.5f}", flush=True)
with open(os.path.join(os.path.dirname(__file__), "out", "dirac_zbw_scan.json"), "w") as fh:
    json.dump(dict(rows=out, runtime_s=time.time() - t0), fh, indent=1)
print(f"runtime {time.time() - t0:.1f} s", flush=True)
