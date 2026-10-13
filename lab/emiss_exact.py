"""emiss_exact.py - EXAM ONLY (never an input of the rule): exact Dirac-Coulomb numbers for the 2p -> 1s decay of a
hydrogen-like ion (point nucleus, infinite nuclear mass, no QED radiative corrections).

  * levels: Sommerfeld-Dirac formula;
  * radial functions: closed form (confluent hypergeometric), checked against the radial Dirac equations;
  * rates:  A_LW  = (4/3) alpha w^3 |<1s|r|2p_j>|^2   (long wavelength, length form, summed over final m)
            A_EX  = (alpha w / 2 pi) Int dOmega_k sum_pol sum_mf |<1s mf| alpha.eps e^{-i k.r} |2p_j mi>|^2
                    (golden rule with the full transverse photon, all multipoles; Coulomb gauge)
            A_NR  = (2/3)^8 alpha^5 Z^4 mc^2/hbar (Schrodinger, the H value 6.2649e8 s^-1 scaled with Z^4;
                    6.2649e8 includes the proton reduced mass).
Run: python lab/emiss_exact.py  ->  lab/out/emiss_exact.json
"""
import json, os
import numpy as np
from scipy.special import hyp1f1, roots_legendre, roots_genlaguerre, sph_harm_y
from zz_core import ALPHA, T_UNIT, MC2_EV
from emiss_core import spin_angle

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def level(n, kappa, Za):
    k = abs(kappa)
    g = np.sqrt(k * k - Za * Za)
    nr = n - k
    E = 1 / np.sqrt(1 + (Za / (nr + g))**2)
    return E, g, nr


def radial(n, kappa, Za, r):
    """G = r g, F = r f with psi = (g Omega_kappa, i f Omega_-kappa), normalised Int (G^2 + F^2) dr = 1."""
    E, g, nr = level(n, kappa, Za)
    lam = np.sqrt(1 - E * E)
    Np = np.sqrt(n * n - 2 * nr * (abs(kappa) - g))
    rho = 2 * lam * r
    A = (Np - kappa) * hyp1f1(-nr, 2 * g + 1, rho)
    B = nr * hyp1f1(1 - nr, 2 * g + 1, rho) if nr > 0 else 0 * rho
    pre = rho**g * np.exp(-rho / 2)
    G = np.sqrt(1 + E) * pre * (A - B)
    F = -np.sqrt(1 - E) * pre * (A + B)
    return E, G, F


def check_ode(n, kappa, Za):
    """residual of dG/dr = -k/r G + (E - V + 1) F,  dF/dr = k/r F - (E - V - 1) G  (V = -Za/r)"""
    r = np.linspace(0.05, 30, 20001) / Za
    E, G, F = radial(n, kappa, Za, r)
    V = -Za / r
    dG, dF = np.gradient(G, r), np.gradient(F, r)
    r1 = dG - (-kappa / r * G + (E - V + 1) * F)
    r2 = dF - (kappa / r * F - (E - V - 1) * G)
    s = np.max(np.abs(G))
    return float(np.max(np.abs(r1[5:-5])) / s * (1 / Za)), float(np.max(np.abs(r2[5:-5])) / s * (1 / Za))


def quad_radial(Za, nq=120):
    """Gauss-Laguerre nodes for Int_0^inf h(r) dr with weight e^{-x}, x = Za r (all states decay at least like e^{-Za r/2})"""
    x, w = roots_genlaguerre(nq, 0.0)
    s = 0.5 * Za                      # slowest decay (2p: e^{-Za r/2}) times 2 for a product of two states
    r = x / s
    return r, w * np.exp(x) / s


def rates(Za, nq=120, nth=40, nph=48):
    out = {}
    r, wr = quad_radial(Za, nq)
    E1, G1, F1 = radial(1, -1, Za, r)
    nrm = np.sum(wr * (G1**2 + F1**2))
    G1, F1 = G1 / np.sqrt(nrm), F1 / np.sqrt(nrm)
    # angular grid
    u, wu = roots_legendre(nth)
    ph = (np.arange(nph) + 0.5) * 2 * np.pi / nph
    wph = 2 * np.pi / nph
    TH, PH = np.meshgrid(np.arccos(u), ph, indexing="ij")
    WA = np.outer(wu, np.full(nph, wph))
    X, Y, Z = np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)
    # Dirac alpha in Dirac representation acting on (upper, lower): alpha_i = [[0, s_i], [s_i, 0]]
    sig = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
    for name, (n2, k2) in (("2p1/2", (2, 1)), ("2p3/2", (2, -2))):
        E2, G2, F2 = radial(n2, k2, Za, r)
        nrm = np.sum(wr * (G2**2 + F2**2))
        G2, F2 = G2 / np.sqrt(nrm), F2 / np.sqrt(nrm)
        w = E2 - E1
        Rlen = float(np.sum(wr * (G1 * G2 + F1 * F2) * r))
        A_LW = 4 / 3 * ALPHA * w**3 * Rlen**2 / 3          # sum over final m of |<Omega|r_hat|Omega>|^2 = 1/3
        # full golden rule, photon along z (isotropy: average over initial m, times 4 pi)
        j2 = abs(k2) - 0.5
        tot = 0.0
        for mi in np.arange(-j2, j2 + 1):
            Om2u = spin_angle(k2, mi, X, Y, Z)
            Om2l = spin_angle(-k2, mi, X, Y, Z)
            for mf in (0.5, -0.5):
                Om1u = spin_angle(-1, mf, X, Y, Z)
                Om1l = spin_angle(1, mf, X, Y, Z)
                for e in (0, 1):            # polarisation x, y
                    s = sig[e]
                    # psi_f^dag alpha.eps psi_i = up_f^dag s lo_i + lo_f^dag s up_i, with lo = i f Omega_-k
                    # angular pieces A1 = Om1u^dag s Om2l (times i F2), A2 = Om1l^dag s Om2u (times -i F1)
                    A1 = np.einsum("a...,ab,b...->...", Om1u.conj(), s, Om2l)
                    A2 = np.einsum("a...,ab,b...->...", Om1l.conj(), s, Om2u)
                    # e^{-i k z} with z = r cos th
                    ph_ = np.exp(-1j * w * np.outer(r, Z.ravel())).reshape((len(r),) + Z.shape)
                    rad1 = (G1 * F2 * 1j)[:, None, None]
                    rad2 = (-1j * F1 * G2)[:, None, None]
                    integ = (rad1 * A1[None] + rad2 * A2[None]) * ph_
                    M = np.sum(wr[:, None, None] * r[:, None, None]**2 / r[:, None, None]**2 * integ * WA[None])
                    tot += abs(M)**2
        A_EX = ALPHA * w / (2 * np.pi) * 4 * np.pi * tot / (2 * j2 + 1)
        out[name] = dict(omega=w, E2=E2, R_length=Rlen, A_LW=A_LW, A_EX=A_EX, A_LW_s=A_LW / T_UNIT, A_EX_s=A_EX / T_UNIT,
                         ratio_EX_LW=A_EX / A_LW)
    out["E1s"] = E1
    return out


if __name__ == "__main__":
    res = {}
    for Z in (1, 60):
        Za = Z * ALPHA
        chk = {f"{n}{k}": check_ode(n, k, Za) for n, k in ((1, -1), (2, 1), (2, -2))}
        rr = rates(Za)
        A_NR = (2 / 3)**8 * ALPHA**5 * Z**4
        rr["A_NR_Z4"] = A_NR
        rr["A_NR_Z4_s"] = A_NR / T_UNIT
        rr["ode_residuals"] = chk
        res[f"Z{Z}"] = rr
        print(f"Z = {Z}: ODE residuals {chk}")
        print(f"   1s E = {rr['E1s']:.8f} (binding {1 - rr['E1s']:.6f})   A_NR (Z^4) = {A_NR:.5e} mc^2/hbar = {A_NR / T_UNIT:.5e} s^-1")
        for nm in ("2p1/2", "2p3/2"):
            q = rr[nm]
            print(f"   {nm}: w = {q['omega']:.6f}, R = {q['R_length']:.5f} lC ({q['R_length'] * Za:.5f} a0), "
                  f"A_LW = {q['A_LW']:.5e} ({q['A_LW'] / A_NR:.4f} x NR), A_EX = {q['A_EX']:.5e} ({q['A_EX'] / A_NR:.4f} x NR), "
                  f"EX/LW = {q['ratio_EX_LW']:.5f}; A_EX = {q['A_EX_s']:.5e} s^-1")
    json.dump(res, open(os.path.join(OUT, "emiss_exact.json"), "w"), indent=1, default=float)
