"""zstar_core.py - PERCHE' Z* = 131,37 (docs/ZSTAR.md).

La regola radiale congelata (zz_core.RadialWalk) vicino al nucleo puntiforme e' un problema SENZA SCALA:
nella cella i (r_i = (i + s) a) gli angoli per tick sono kappa/(i+s) (centrifugo) e Z alpha/(i+s) (fase di
Coulomb), indipendenti da a; massa (m a) ed energia (omega a) sono trascurabili per r << 1. Quindi lontano
dalle prime celle (a << r << 1) la soluzione regolare della regola e' una combinazione delle due soluzioni del
continuo a energia zero, r^{+gamma} e r^{-gamma} (gamma = sqrt(kappa^2 - zeta^2)):

    g(r) ~ A (r/a)^{+gamma} + B (r/a)^{-gamma},    beta(zeta) = B/A   (numero della regola, senza a)

e la regola equivale al Dirac del continuo con la condizione al contorno "r^gamma + tau r^-gamma",
tau = beta a^{2 gamma} (estensione autoaggiunta effettiva). Previsione: il livello si tuffa a -mc^2 quando
tau = tau_cr(zeta) (numero del continuo); per a -> 0 serve beta -> infinito, cioe' A(zeta*) = 0: la regola,
a energia zero, ha un modo che decade come r^{-gamma} (uno "stato di contatto" del reticolo).

Strumenti:
  * zero_mode(): ricorrenza verso fuori a lambda = 1, m = 0 (esatta, senza a), con offset s, refl0 complesso,
    varianti della moneta.
  * subdominant(): la soluzione r^{-gamma} per ricorrenza all'indietro (stabile) da L con lo spinore del continuo.
  * D(zeta): "wronskiano" det[s_reg, s_sub] in una cella vicina al nucleo; zero <=> A = 0. zeta* = radice.
  * beta_fit(): beta dal fit della soluzione regolare (controllo, e per la previsione di Z_cr(a)).
  * tau_cr(): continuo (SOLO CONTROLLO/raccordo): tau per cui il Dirac continuo con m = 1, E = -1 ha lo stato
    di soglia (G(r_end) cambia segno), condizione G = r^g + tau r^-g a r0.
Unita' hbar = c = m_e = 1.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA, U_CHI
from supercrit_core import coins, njit

UD = U_CHI.conj().T


def coin_variant(kappa, a, r, V, mass, variant="frozen"):
    """chi-basis coin per site. 'frozen' = zz_core (exp of the sum). 'split' = e^{-iVa} e^{-i kappa/r sx a}
    e^{-i m sz a} (Trotter, order changed). Both are identical for m = 0."""
    if variant == "frozen":
        return coins(kappa, a, r, V, mass)
    dt = a
    c1, s1 = np.cos(dt * kappa / r), np.sin(dt * kappa / r)
    Mx = np.zeros((len(r), 2, 2), complex)
    Mx[:, 0, 0] = c1; Mx[:, 1, 1] = c1; Mx[:, 0, 1] = -1j * s1; Mx[:, 1, 0] = -1j * s1
    Mz = np.zeros((len(r), 2, 2), complex)
    Mz[:, 0, 0] = np.exp(-1j * mass * dt); Mz[:, 1, 1] = np.exp(1j * mass * dt)
    M = np.einsum("nab,nbc->nac", Mx, Mz) * np.exp(-1j * V * dt)[:, None, None]
    return np.einsum("ab,nbc,cd->nad", UD, M, U_CHI)


def scale_free_coins(zeta, kappa, s, L, variant="frozen", vfun=None):
    """coins of the a-free problem: a = 1, m = 0, sites rho = i + s. vfun(rho, zeta) overrides V (default -zeta/rho)."""
    rho = np.arange(L) + s
    V = -zeta / rho if vfun is None else vfun(rho, zeta)
    return rho, coin_variant(kappa, 1.0, rho, V, 1e-300, variant)


@njit(cache=True)
def _fwd(M, refl0, lam):
    """outward recursion of W chi = lam chi (W = M S). state s_i = (chi+_{i-1}, chi-_i); returns chip, chim
    (chi+_i, chi-_i) and log-scale; normalised chi-_0 = 1 (rescaling tracked)."""
    n = M.shape[0]
    chip = np.zeros(n, np.complex128)
    chim = np.zeros(n + 1, np.complex128)
    lsc = np.zeros(n + 1)
    chim[0] = 1.0
    prev = refl0 * chim[0]
    acc = 0.0
    for i in range(n):
        # A = lam M^dag ; (p_i, q_i) = A chi_i with p_i = chi+_{i-1}, q_i = chi-_{i+1}
        a00 = lam * np.conj(M[i, 0, 0]); a01 = lam * np.conj(M[i, 1, 0])
        a10 = lam * np.conj(M[i, 0, 1]); a11 = lam * np.conj(M[i, 1, 1])
        chip[i] = (prev - a01 * chim[i]) / a00
        chim[i + 1] = a10 * chip[i] + a11 * chim[i]
        prev = chip[i]
        sc = abs(chip[i]) + abs(chim[i + 1])
        if sc > 1e50:
            chip[:i + 1] /= sc
            chim[:i + 2] /= sc
            prev /= sc
            acc += np.log(sc)
        lsc[i + 1] = acc
    return chip, chim[:n], lsc[:n]


@njit(cache=True)
def _bwd(M, lam, xp_last, xm_last, i_stop):
    """inward recursion from the last site: given (chi+_{L-1}, chi-_{L-1}) at site L-1 obtain chi at sites < L-1.
    chi_i = M_i (p_i, q_i), p_i = chi+_{i-1}, q_i = chi-_{i+1}  ->  (p_i, q_i) = lam M_i^dag chi_i,
    chi-_{i+1} known, chi+_{i} known -> solve for chi-_i, chi+_{i-1}."""
    n = M.shape[0]
    chip = np.zeros(n, np.complex128)
    chim = np.zeros(n, np.complex128)
    chip[n - 1] = xp_last
    chim[n - 1] = xm_last
    for i in range(n - 1, i_stop, -1):
        # at site i: q_i = chi-_{i+1} unknown only for i = n-1 (use asymptotic: chi- continuity); we instead use the
        # pair relation of site i: chi_i = M_i (p_i, q_i) -> (p_i, q_i) = M_i^dag chi_i / lam* ... (unitary M)
        a00 = np.conj(M[i, 0, 0]); a01 = np.conj(M[i, 1, 0])
        a10 = np.conj(M[i, 0, 1]); a11 = np.conj(M[i, 1, 1])
        p = lam * (a00 * chip[i] + a01 * chim[i])       # = chi+_{i-1}
        # q = lam * (a10 chip[i] + a11 chim[i]) = chi-_{i+1} (consistency, not needed)
        chip[i - 1] = p
        # chi-_{i-1}: from site i-1 relation: q_{i-1} = chi-_i = lam (a10' chip[i-1] + a11' chim[i-1])
        b10 = np.conj(M[i - 1, 0, 1]); b11 = np.conj(M[i - 1, 1, 1])
        chim[i - 1] = (chim[i] / lam - b10 * chip[i - 1]) / b11
        sc = abs(chip[i - 1]) + abs(chim[i - 1])
        if sc > 1e50:
            chip[i - 1:] /= sc
            chim[i - 1:] /= sc
    return chip, chim


def cont_spinor(zeta, kappa, sign):
    """continuum zero-energy massless solution r^{s}, s = sign*gamma: (G, F) = (1, (s + kappa)/zeta)"""
    gam = np.sqrt(kappa**2 - zeta**2)
    sg = sign * gam
    return np.array([1.0, (sg + kappa) / zeta]), sg


def gf_to_chi(G, F):
    """(g, f) -> (chi+, chi-): g = (p+q)/sqrt2, f = i (p-q)/sqrt2"""
    return np.array([(G - 1j * F) / np.sqrt(2), (G + 1j * F) / np.sqrt(2)])


def D_wronskian(zeta, kappa=-1, s=0.5, refl0=-1.0, L=200000, i_eval=3, variant="frozen", vfun=None,
                lam=1.0, return_parts=False):
    """det[ s_reg(i), s_sub(i) ] at site i_eval, s = (chi+_{i-1}, chi-_i), normalised by |s_reg||s_sub|.
    Zero <=> regular lattice solution has no r^{+gamma} part (A = 0)."""
    rho, M = scale_free_coins(zeta, kappa, s, L, variant, vfun)
    cp, cm, _ = _fwd(M, complex(refl0), complex(lam))
    e, sg = cont_spinor(zeta, kappa, -1)
    x = gf_to_chi(e[0], e[1]) * rho[-1]**sg
    # 1/rho correction ignored: its r^{+gamma} contamination dies inward as (rho/L)^{2 gamma}
    bp, bm = _bwd(M, complex(lam), x[0], x[1], max(i_eval - 2, 0))
    i = i_eval
    sr = np.array([cp[i - 1], cm[i]])
    ss = np.array([bp[i - 1], bm[i]])
    d = (sr[0] * ss[1] - sr[1] * ss[0]) / (np.linalg.norm(sr) * np.linalg.norm(ss))
    if return_parts:
        return d, (cp, cm), (bp, bm), rho
    return d


def D_real(zeta, **kw):
    """D with its constant phase removed: unitarity (current conservation) makes D = e^{i phi} x real.
    phi is fixed by the convention: we rotate by the phase of D computed with the r^{-gamma} start replaced by
    r^{+gamma}... simpler: return Re(D * conj(ph)) with ph = D/|D| evaluated at a reference zeta_ref (kw)."""
    ref = kw.pop("zeta_ref", None)
    d = D_wronskian(zeta, **kw)
    if ref is None:
        return d
    d0 = D_wronskian(ref, **kw)
    return (d * np.conj(d0) / abs(d0)).real


def beta_fit(zeta, kappa=-1, s=0.5, refl0=-1.0, L=100000, rmin=40.0, ncorr=3, variant="frozen", vfun=None):
    """regular solution g(rho) = rho^g (A0 + A1/rho + ..) + rho^-g (B0 + B1/rho + ..) fitted on rho in [rmin, L]
    (both g and f, sharing the leading coefficients through the continuum spinors). Returns beta = B0/A0 (real
    part after removing the common phase), A0, B0, fit residual."""
    rho, M = scale_free_coins(zeta, kappa, s, L, variant, vfun)
    cp, cm, lsc = _fwd(M, complex(refl0), 1.0 + 0j)
    # undo the rescaling only within the fit window: work with the window's own scale
    g = (cp + cm) / np.sqrt(2)
    f = 1j * (cp - cm) / np.sqrt(2)
    sel = np.where(rho >= rmin)[0]
    sel = sel[lsc[sel] == lsc[sel[-1]]]          # same scale
    sel = sel[:: max(len(sel) // 4000, 1)]
    r = rho[sel]
    gam = np.sqrt(kappa**2 - zeta**2)
    cols = []
    for k in range(ncorr + 1):
        cols.append(r**gam / r**k)
    for k in range(ncorr + 1):
        cols.append(r**(-gam) / r**k)
    X = np.stack(cols, 1)
    cg, *_ = np.linalg.lstsq(X / r[:, None]**gam, g[sel] / r**gam, rcond=None)
    res = np.abs(X @ cg - g[sel]).max() / np.abs(g[sel]).max()
    A0, B0 = cg[0], cg[ncorr + 1]
    ph = A0 / abs(A0)
    return (B0 / A0).real, A0, B0, res, (B0 / A0).imag


# ----------------------------------------------------------------- continuum matching (CHECK / raccordo)
def _rhs(kappa, E, zeta, m=1.0):
    def f(r, y):
        V = -zeta / r
        G, F = y
        return [-kappa / r * G + (E - V + m) * F, -(E - V - m) * G + kappa / r * F]
    return f


def cont_G_end(zeta, tau, kappa=-1, E=-1.0, r0=1e-7, r_end=40.0):
    """continuum Dirac (m = 1) at energy E with small-r condition G = r^g + tau r^-g (zero-energy spinors);
    returns G(r_end) normalised. For |tau| = inf use the pure r^-g."""
    from scipy.integrate import solve_ivp
    ep, gp = cont_spinor(zeta, kappa, +1)
    em, gm = cont_spinor(zeta, kappa, -1)
    if np.isinf(tau):
        y0 = em * r0**gm
    else:
        y0 = ep * r0**gp + tau * em * r0**gm
    n0 = np.linalg.norm(y0)
    sol = solve_ivp(_rhs(kappa, E, zeta), (r0, r_end), y0 / n0, rtol=1e-11, atol=1e-30, method="DOP853")
    return sol.y[0, -1]


def zeta_cr_of_a(a, kappa, s, refl0, zlo, zhi, beta_fun=None, r_end=40.0, tol=1e-7):
    """prediction of the rule's Z_cr(a): zeta where the continuum with tau = beta(zeta) a^{2 gamma} has its
    threshold state at E = -1 (first sign change of G(r_end) scanning up from zlo)."""
    from scipy.optimize import brentq
    if beta_fun is None:
        beta_fun = lambda z: beta_fit(z, kappa, s, refl0)[0]

    def h(z):
        gam = np.sqrt(kappa**2 - z**2)
        b = beta_fun(z)
        # write the BC as cos(t) r^g + sin(t) r^-g to avoid the pole of beta
        tau = b * a**(2 * gam)
        return cont_G_end(z, tau, kappa, r_end=r_end)
    zs = np.linspace(zlo, zhi, 25)
    hv = [h(z) for z in zs]
    for k in range(len(zs) - 1):
        if np.sign(hv[k]) != np.sign(hv[k + 1]):
            return brentq(h, zs[k], zs[k + 1], xtol=tol)
    return None


if __name__ == "__main__":
    # sanity: nucleo1s_ (NUCLEO_1S.md, origin/master) measured beta(Z=60) = 0.10301 on the same rule
    for Z in (20, 60, 92):
        b = beta_fit(Z * ALPHA)
        print(f"Z={Z}: beta = {b[0]:.5f} (imag {b[4]:.1e}, fit res {b[3]:.1e})")
    for z in (0.95, 0.958, 0.959, 0.97):
        d = D_wronskian(z)
        print(f"zeta={z}: D = {d:.6f}")


# ----------------------------------------------------------------- full prediction of Z_cr(a) by matching
def lattice_AB(zeta, kappa=-1, s=0.5, refl0=-1.0, L=100000, rmin=40.0, ncorr=3, variant="frozen", vfun=None):
    """leading coefficients (A0, B0) of g = A0 rho^gamma + B0 rho^-gamma (complex gamma allowed for zeta > |kappa|)"""
    rho, M = scale_free_coins(zeta, kappa, s, L, variant, vfun)
    cp, cm, lsc = _fwd(M, complex(refl0), 1.0 + 0j)
    g = (cp + cm) / np.sqrt(2)
    sel = np.where(rho >= rmin)[0]
    sel = sel[lsc[sel] == lsc[sel[-1]]]
    sel = sel[:: max(len(sel) // 4000, 1)]
    r = rho[sel]
    gam = np.sqrt(kappa**2 - zeta**2 + 0j)
    if abs(gam) < 1e-4:
        raise ValueError("gamma too close to 0")
    cols = [r**gam / r**k for k in range(ncorr + 1)] + [r**(-gam) / r**k for k in range(ncorr + 1)]
    X = np.stack(cols, 1)
    sc = np.abs(r**gam)
    c, *_ = np.linalg.lstsq(X / sc[:, None], g[sel] / sc, rcond=None)
    return c[0], c[ncorr + 1], gam


def matching_G_end(zeta, a, kappa=-1, s=0.5, refl0=-1.0, r0=1e-8, r_end=40.0, ref=None, variant="frozen",
                   vfun=None):
    """continuum Dirac (m = 1, E = -1) started at r0 with the lattice's own small-r solution
    A0 (r/a)^g e+ + B0 (r/a)^-g e-; returns (G(r_end) normalised, phase-half-angle used) -- real solution."""
    from scipy.integrate import solve_ivp
    A0, B0, gam = lattice_AB(zeta, kappa, s, refl0, variant=variant, vfun=vfun)
    ep = np.array([1.0, (gam + kappa) / zeta])
    em = np.array([1.0, (-gam + kappa) / zeta])
    cA = A0 * (r0 / a)**gam
    cB = B0 * (r0 / a)**(-gam)
    y0 = cA * ep + cB * em
    if abs(gam.imag) < 1e-12:          # zeta < |kappa|: common phase of A0, B0 (constant ~ i, sign fixed by ref)
        ph = (B0 if abs(B0) > abs(A0) else A0)
        ph = ph / abs(ph)
        if ref is not None and (ph * np.conj(ref)).real < 0:
            ph = -ph
        yr = (y0 * np.conj(ph))
        half = ph
    else:                               # zeta > |kappa|: |cA| = |cB| (current), real solution = e^{i(argA+argB)/2} x real
        half = np.exp(0.5j * (np.angle(cA) + np.angle(cB)))
        if ref is not None and (half * np.conj(ref)).real < 0:
            half = -half
        yr = y0 * np.conj(half)
    imag_err = np.abs(yr.imag).max() / np.abs(yr).max()
    yr = yr.real / np.abs(yr).max()
    sol = solve_ivp(_rhs(kappa, -1.0, zeta), (r0, r_end), yr, rtol=1e-11, atol=1e-30, method="DOP853")
    G = sol.y[0, -1] / np.abs(sol.y[:, -1]).max()
    return np.sign(G) * min(abs(G), 1.0), half, imag_err


def predict_Zcr(a, kappa=-1, s=0.5, refl0=-1.0, zlo=None, zhi=None, n=60, variant="frozen", vfun=None, tol=1e-6):
    """first zeta (scanning upward) where the matched continuum has its threshold state at E = -1.
    returns zeta_cr (or None), max imag error of the real-structure check"""
    k = abs(kappa)
    zlo = 0.3 * k if zlo is None else zlo
    zhi = 1.25 * k if zhi is None else zhi
    zs = np.linspace(zlo, zhi, n)
    zs = zs[np.abs(zs - k) > 2e-3]
    vals, ref, errs = [], None, []
    halves = []
    for z in zs:
        v, ref, e = matching_G_end(z, a, kappa, s, refl0, ref=ref, variant=variant, vfun=vfun)
        vals.append(v); errs.append(e); halves.append(ref)
    vals = np.array(vals)
    for j in range(len(zs) - 1):
        if np.sign(vals[j]) != np.sign(vals[j + 1]) and abs(vals[j] - vals[j + 1]) < 1.5:
            lo, hi, rlo = zs[j], zs[j + 1], halves[j]
            flo = vals[j]
            if (lo - k) * (hi - k) < 0:
                continue
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                fm, hm, _e = matching_G_end(mid, a, kappa, s, refl0, ref=rlo, variant=variant, vfun=vfun)
                if np.sign(fm) == np.sign(flo):
                    lo, flo, rlo = mid, fm, hm
                else:
                    hi = mid
                if hi - lo < tol:
                    break
            return 0.5 * (lo + hi), max(errs)
    return None, max(errs)
