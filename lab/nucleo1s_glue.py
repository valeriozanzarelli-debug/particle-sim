"""nucleo1s_glue.py - a stronger glue for the radial hybrid of phase 1 (zig-zag zoom inside, learned agglomerate outside).

Phase-1 glue: ghosts j < J0 = the inner zig-zag solution, coarse equations j >= J0, energy from the coarse equation of the
last ghost row.  Its weakness (found in phase 2): the coarse kinetic operator of the agglomerate is spectral (T(k) = k^2/2m
up to k = pi, kernel with an alternating 1/n^2 tail); applied to the odd extension of u = r psi, which has the Coulomb
cusp (u'' jumps at r = 0), it makes an alternating aliasing error that decays slowly with the distance from the nucleus:
the energy then jumps with the parity of J0 (+0.10 %, -0.02 %, +0.12 % for J0 = 2, 3, 4).

Defect (tau) correction, as in full-approximation multigrid: the inner rule tells the coarse rule how wrong its kinetic
operator is on the true solution near the nucleus,
        tau_j = [K (w u_in)]_j - T_c(w u_in)(r_j),     T_c v = -v'' / (2 m_kin),
with w a C-infinity window (1 up to J0, 0 after J0 + window), and the coarse equations become (H - E) u = tau.
tau is linear in the amplitude of u_in (same E), so the problem stays an eigen-problem in E.  Nothing is fitted.
"""
import os, sys
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_hybrid import CoarseRadial, Hybrid, zz_regular, A0  # noqa: E402
from nucleo1s_core import dirac_binding, ALPHA  # noqa: E402


def planck(x):
    x = np.clip(x, 1e-12, 1 - 1e-12)
    with np.errstate(over="ignore"):
        return 1 / (1 + np.exp(1 / (1 - x) - 1 / x))


def solve_glue(c, Z, nst, J0, window=12, tau=True, rmax=None):
    Za = Z * ALPHA
    Eref = -dirac_binding(nst, -1, Za)
    rmax = rmax or (24 if nst == 1 else 40) * A0 / Z
    Jw = J0 + window
    hy = Hybrid(c, rmax, J0=Jw + 3, Za=Za)
    N, h = hy.N, c.h
    K = hy.Hm - np.diag(c.V(N, Za))
    rj = hy.r
    Rw1, Rw2 = (J0 + 0.5) * h, (Jw + 0.5) * h

    def wf(rr):
        return np.where(rr <= Rw1, 1.0, np.where(rr >= Rw2, 0.0, 1 - planck((rr - Rw1) / (Rw2 - Rw1))))

    def parts(E):
        r, g, f = zz_regular(hy.mesh, -1, Za, 1 + E)
        G = g.real
        inside = rj <= r[-1]
        ui = np.where(inside, CubicSpline(r, G)(np.minimum(rj, r[-1])), 0.0)
        v = wf(rj) * ui
        Tc = -CubicSpline(r, wf(r) * G)(np.minimum(rj, r[-1]), 2) / (2 * c.m_kin) * inside
        return ui, (K @ v - Tc) if tau else np.zeros(N)

    def res(E):
        ui, t = parts(E)
        A = hy.Hm - E * np.eye(N)
        ug = ui[:J0]
        uo = np.linalg.solve(A[J0:, J0:], -A[J0:, :J0] @ ug + t[J0:])
        return (A[J0 - 1, :J0] @ ug + A[J0 - 1, J0:] @ uo - t[J0 - 1]) / np.abs(ug).max()
    Es = np.linspace(Eref * 1.03, Eref * 0.97, 13)
    rs = [res(E) for E in Es]
    br = [i for i in range(12) if np.sign(rs[i]) != np.sign(rs[i + 1])]
    if not br:
        return None
    i = min(br, key=lambda i: abs(0.5 * (Es[i] + Es[i + 1]) - Eref))
    E = brentq(res, Es[i], Es[i + 1], xtol=1e-18, rtol=1e-14)
    return E / Eref - 1


def glue_scan(cfg, log=print):
    c = CoarseRadial(cfg["B"], vmode=cfg["vmode"])
    out = dict(cfg=cfg)
    for Z, nst in cfg["states"]:
        rec = {}
        for J0 in cfg["J0"]:
            rec[f"J0={J0}"] = solve_glue(c, Z, nst, J0, cfg["window"], tau=True)
            rec[f"J0={J0}_phase1"] = solve_glue(c, Z, nst, J0, cfg["window"], tau=False)
        out[f"Z{Z}_{nst}s"] = rec
        log(f"Z{Z} {nst}s", {k: (None if v is None else f"{100 * v:+.4f}%") for k, v in rec.items()})
    return out
