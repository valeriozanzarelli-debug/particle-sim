"""fluid_f3_post.py - POST-HOC (not pre-registered) diagnosis of the E3c miss: same driven kink, T = 600 instead of 150,
chain 1400 long. Hypothesis: relativistic approach to terminal speed takes ~ gamma^2/alpha, T = 150 too short. RESULT: rejected (same v at F = 0.2, 0.4).
Writes out/fluid_f3_post.json."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fluid_core import fk_acc, kink_profile
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def driven(h, F, alpha=0.05, T=600.0, Lx=1400.0):
    n = int(Lx / h); x = (np.arange(n) - n / 2) * h
    phi = kink_profile(x, 0.6 * Lx / 2, 1.0) + np.arcsin(F); dphi = np.zeros_like(phi)
    dt = 0.2 * h; a = fk_acc(phi, h, F, alpha, dphi); ts, Xs = [], []; t = 0.0
    for s in range(int(round(T / dt))):
        dphi += 0.5 * dt * a; phi += dt * dphi; a = fk_acc(phi, h, F, alpha, dphi); dphi += 0.5 * dt * a; t += dt
        if s % int(round(1 / dt)) == 0 and t > T - 100:
            c = np.where(np.diff(np.floor((phi - np.pi - np.arcsin(F)) / (2 * np.pi))) != 0)[0]
            ts.append(t); Xs.append(x[c[0]])
    v = float(abs(np.polyfit(ts, Xs, 1)[0])); vms = 1 / np.sqrt(1 + (4 * alpha / (np.pi * F)) ** 2)
    return dict(F=F, v=v, v_MS=float(vms), rel=float(v / vms - 1))


res = [driven(0.25, F) for F in [0.05, 0.1, 0.2, 0.4, 0.8]]
json.dump(res, open(os.path.join(OUT, "fluid_f3_post.json"), "w"), indent=1)
for r in res: print(r)
