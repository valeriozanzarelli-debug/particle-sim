"""box_threshold.py - A2: the photon threshold in the strip. A weak oscillating field A(t) = a0 sin(Omega t) along x
(a 'photon' of frequency Omega, uniform in space) acts on the sea for T = 600 ticks; the excitation summed over a kx
grid is the absorption spectrum P(Omega). Threshold Omega* = position of its maximum (the 1+1 D van Hove edge just
above the gap). Prediction: Omega* = 2 m_eff(d) within 10%, rising above 2 th as d shrinks.
Usage: python box_threshold.py d [kstep]  (kx grid step, default 0.004; 0.008 / 0.016 used for d = 32 / 128 to fit
the 6-minute budget: the Fejer width 2 pi/T = 0.010 in Omega is far larger than the step's effect on Omega)   (writes lab/out/box_threshold_d{d}.json; d = 128 takes ~4 min)"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from box_core import static_part, masses, floquet, sx_phase, W_CUT, PI

TH, A0, T, DOM = 0.1, 0.01, 600, 0.0025
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def spectrum(d, kstep=0.004):
    M = static_part(TH, d)
    L = M.shape[0] // 2
    m = masses(TH, d, M)
    m0 = float(m[0])
    om = np.arange(1.6 * TH, 2.6 * m0 + 1e-9, DOM)
    kmax = 1.5 * m0
    ks = np.arange(-kmax, kmax + 1e-9, kstep)
    cols, kcol, ups = [], [], []
    for k in ks:
        w, v = floquet(M, k)
        sel = (w < 0) & (w > -(om[-1] - m0 + 0.02))
        for j in np.where(sel)[0]:
            cols.append(v[:, j]); kcol.append(k)
        ups.append(v[:, (w > 0) & (w < W_CUT)])
    psi0 = np.array(cols).T
    kcol = np.array(kcol)
    kidx = np.searchsorted(ks, kcol - 1e-12)
    P = []
    for Om in om:
        psi = psi0.copy()
        for t in range(T):
            A = A0 * np.sin(Om * t)
            psi = M @ (sx_phase(kcol + A, L) * psi)
        A_end = A0 * np.sin(Om * T)
        tot = 0.0
        for i, k in enumerate(ks):
            c = kidx == i
            if not c.any():
                continue
            if abs(A_end) > 1e-12:
                w, v = floquet(M, k + A_end)
                U = v[:, (w > 0) & (w < W_CUT)]
            else:
                U = ups[i]
            tot += float((np.abs(U.conj().T @ psi[:, c]) ** 2).sum())
        P.append(tot * (ks[1] - ks[0]) / (2 * PI))
    P = np.array(P)
    i = int(np.argmax(P))
    # onset: first Omega where P exceeds 10% of the maximum (reported, not graded)
    onset = float(om[np.where(P > 0.1 * P[i])[0][0]])
    return dict(d=d, m_eff=m0, two_m_eff=2 * m0, Omega=om.tolist(), P=P.tolist(), Omega_star=float(om[i]),
                ratio=float(om[i] / (2 * m0)), onset_10pct=onset, ncols=int(psi0.shape[1]), nk=len(ks))


if __name__ == "__main__":
    for d in [int(sys.argv[1])]:
        t0 = time.time()
        ksp = float(sys.argv[2]) if len(sys.argv) > 2 else 0.004
        r = spectrum(d, ksp)
        r["kstep"] = ksp
        r["time_s"] = time.time() - t0
        json.dump(r, open(os.path.join(OUT, f"box_threshold_d{d}.json"), "w"), indent=1)
        print(d, "2m_eff", round(r["two_m_eff"], 4), "Omega*", round(r["Omega_star"], 4), "ratio", round(r["ratio"], 4),
              "onset", round(r["onset_10pct"], 4), r["ncols"], f"{r['time_s']:.0f}s", flush=True)
