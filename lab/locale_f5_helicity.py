"""locale_f5_helicity.py - phase 5, E6 (information, measure made concrete after the run: the 'purity of the + band'
of the prereg is ill-defined because the + band is two-fold degenerate in the doubled walk). Test: a helicity-+ photon
of the frozen light walk (positive frequency, all in F+) evolved with the breathing light walk; maximum weight that leaks
into F- over 2000 periods. Writes lab/out/locale_f5_helicity.json."""
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locale_core import SplitWalk, sdot, zigzag_ticks

X6 = np.kron(np.array([[0, 1], [1, 0]], complex), np.eye(3))


def lwalk(mu):
    fac = []
    for it, tk in enumerate(zigzag_ticks(3)):
        s = -1 if it == 0 else +1
        for d in tk:
            S = sdot(d, 1); M = np.zeros((6, 6), complex); M[:3, :3] = S; M[3:, 3:] = -S
            fac.append(("shift", d, M))
        if mu:
            fac.append(("coin", np.cos(s * mu) * np.eye(6) - 1j * np.sin(s * mu) * X6))
    return SplitWalk(fac, 6, np.eye(3), T_phys=2.0, dim=3)


out = {}
w0 = lwalk(0.0)
for mu in (0.0, 0.1, 0.3, 0.5):
    w = lwalk(mu)
    for N in (40, 20):
        kk = 2 * np.pi / N
        K = np.array([[0.3, 0.5, 0.81]]); K = K / np.linalg.norm(K) * kk
        e0, v0 = w0.phases(K)
        cand = [j for j in range(6) if np.sum(abs(v0[0][:3, j]) ** 2) > 0.99 and e0[0][j] > kk]
        j = min(cand, key=lambda j: abs(e0[0][j] - 2 * kk))
        psi = v0[0][:, j] / np.linalg.norm(v0[0][:, j])
        W = w.bloch(K)[0]
        mx, mean = 0.0, 0.0
        for t in range(2000):
            psi = W @ psi
            f = float(np.sum(abs(psi[3:]) ** 2)); mx = max(mx, f); mean += f / 2000
        out[f"{mu}|N{N}"] = dict(max_leak=mx, mean_leak=mean)
        print(mu, N, f"max leak into F- {mx:.3e}, mean {mean:.3e}", flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "locale_f5_helicity.json"), "w"), indent=1)
