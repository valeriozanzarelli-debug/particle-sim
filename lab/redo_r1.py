"""redo_r1.py - R1: pairs broken out of the vacuum of the 1+1 D zig-zag walk by a strong field (Schwinger, Sauter).
Pre-registered in lab/out/redo_prereg.json (R1). Writes lab/out/redo_r1.json.  (~1 min)

Walk (one tick): R(x) <- e^{-iA} R(x-1), L(x) <- e^{+iA} L(x+1)  (Peierls phase A), then coin with potential
    R' = e^{-iV}(c R - i s L),  L' = e^{-iV}(c L - i s R),   c = cos th, s = sin th.
k-space: U(k) = C diag(e^{-i(k+A)}, e^{i(k+A)}) on (R, L); cos w = cos th cos k.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, sys, time
import numpy as np
from scipy.integrate import solve_ivp, quad

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
PI = np.pi


def coin(th):
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -1j * s], [-1j * s, c]])


def Uk(k, th, V=0.0):
    return np.exp(-1j * V) * coin(th) @ np.diag([np.exp(-1j * k), np.exp(1j * k)])


def bands(k, th, V=0.0):
    """(lower, upper) Floquet eigenvectors of U(k): quasi-energy w = -angle(lambda) - V in (-pi,0) / (0,pi)"""
    lam, v = np.linalg.eig(Uk(k, th, 0.0))
    w = -np.angle(lam)
    o = np.argsort(w)
    return v[:, o[0]], v[:, o[1]], w[o]


# ------------------------------------------------------------------ uniform field
def lz(th, E, k_start, ticks):
    """exact evolution of one sea plane wave under A(t) = E t; returns upper-band probability at the end"""
    lo, _, _ = bands(k_start, th)
    psi = lo.copy()
    C = coin(th)
    k = k_start
    for n in range(ticks):
        psi = C @ (np.array([np.exp(-1j * k), np.exp(1j * k)]) * psi)
        k += E
    _, up, _ = bands(k, th)
    return float(abs(np.vdot(up, psi)) ** 2)


def crossing(th, E, centre):
    ticks = int(round(PI / E))
    return lz(th, E, centre - PI / 2, ticks)


def dykhne(th, E):
    sc = np.arccosh(1 / np.cos(th))
    I = quad(lambda s: np.arccos(np.clip(np.cos(th) * np.cosh(s), -1, 1)), 0, sc, epsabs=1e-14, epsrel=1e-12)[0]
    return float(np.exp(-4 * I / E))


def real_space_check(th=0.1, N=128, E=0.02, ticks=None):
    """whole sea on a ring of N sites, Peierls phases in real space; compare upper population per k with k-space"""
    ticks = ticks or int(round(PI / E / 2))
    ks = 2 * PI * np.arange(N) / N
    x = np.arange(N)
    C = coin(th)
    err = 0.0
    psi = np.zeros((2, N, N), complex)                 # (comp, x, state)
    for j, k in enumerate(ks):
        lo, _, _ = bands(k, th)
        psi[:, :, j] = lo[:, None] * np.exp(1j * k * x)[None, :] / np.sqrt(N)
    A = 0.0
    for n in range(ticks):
        R = np.exp(-1j * A) * np.roll(psi[0], 1, axis=0)
        L = np.exp(1j * A) * np.roll(psi[1], -1, axis=0)
        psi = np.stack([C[0, 0] * R + C[0, 1] * L, C[1, 0] * R + C[1, 1] * L])
        A += E
    for j, k in enumerate(ks):
        _, up, _ = bands(k + A, th)
        amp = (up.conj()[:, None] * psi[:, :, j]).sum(0) @ np.exp(-1j * k * x) / np.sqrt(N)
        p_real = abs(amp) ** 2
        # k-space
        lo, _, _ = bands(k, th)
        ph = lo.copy(); kk = k
        for n in range(ticks):
            ph = C @ (np.array([np.exp(-1j * kk), np.exp(1j * kk)]) * ph); kk += E
        p_k = abs(np.vdot(up, ph)) ** 2
        err = max(err, abs(p_real - p_k))
    return float(err)


# ------------------------------------------------------------------ Sauter step
def Mx(phi, th):
    c, s = np.cos(th), np.sin(th)
    return np.array([[np.exp(1j * phi) / c, -1j * s / c], [1j * s / c, np.exp(-1j * phi) / c]])


def modes(M, Jz):
    lam, v = np.linalg.eig(M)
    f = np.real(np.einsum("ia,ij,ja->a", v.conj(), Jz, v))
    return lam, v, f


def transmission(T, Vl, Vr, ML, MR):
    """T = total transfer matrix (left end -> right end); ML/MR one-site transfer matrices of the flat ends.
    Incident right-moving from the left. Returns transmitted flux fraction (= pair number in the Klein zone)."""
    Jz = np.diag([1.0, -1.0])
    lamL, vL, fL = modes(ML, Jz)
    lamR, vR, fR = modes(MR, Jz)
    if not (np.allclose(abs(lamL), 1, atol=1e-9) and np.allclose(abs(lamR), 1, atol=1e-9)):
        return float("nan")
    i_in, i_ref = (0, 1) if fL[0] > 0 else (1, 0)
    i_out = 0 if fR[0] > 0 else 1
    # T (v_in + r v_ref) = t v_out
    A = np.stack([T @ vL[:, i_ref], -vR[:, i_out]], 1)
    r, t = np.linalg.solve(A, -T @ vL[:, i_in])
    return float(abs(t) ** 2 * fR[i_out] / fL[i_in])


def lattice_N(w, th, V0, ell):
    X = int(25 * ell + 50)
    xs = np.arange(-X, X + 1)
    V = 0.5 * V0 * (1 + np.tanh(xs / ell))
    T = np.eye(2, dtype=complex)
    for v in V:
        T = Mx(w - v, th) @ T
    return transmission(T, 0.0, V0, Mx(w, th), Mx(w - V0, th))


def continuum_N(w, m, V0, ell):
    X = 25 * ell + 50
    Vf = lambda x: 0.5 * V0 * (1 + np.tanh(x / ell))
    sz, sx = np.diag([1.0, -1.0]), np.array([[0, 1.0], [1.0, 0]])
    def rhs(x, y):
        Y = y[:4] + 1j * y[4:]
        Y = Y.reshape(2, 2)
        K = 1j * sz @ ((w - Vf(x)) * np.eye(2) - m * sx)
        d = (K @ Y).ravel()
        return np.r_[d.real, d.imag]
    y0 = np.r_[np.eye(2).ravel(), np.zeros(4)]
    sol = solve_ivp(rhs, (-X, X), y0, rtol=1e-10, atol=1e-12, method="DOP853")
    Y = (sol.y[:4, -1] + 1j * sol.y[4:, -1]).reshape(2, 2)
    # flat ends: exp(K) over unit length acts like a one-site transfer matrix (same eigenvectors)
    from scipy.linalg import expm
    KL = 1j * sz @ (w * np.eye(2) - m * sx); KR = 1j * sz @ ((w - V0) * np.eye(2) - m * sx)
    return transmission(Y, 0.0, V0, expm(KL), expm(KR))


def sauter_formula(w, m, V0, ell):
    pL = np.sqrt(w**2 - m**2); pR = np.sqrt((V0 - w)**2 - m**2)
    a = PI * ell / 2
    # FIX (declared): first version had this ratio upside down (N > 1); caught by the pre-registered ODE check.
    return float(np.sinh(2 * a * pR) * np.sinh(2 * a * pL) / (np.sinh(a * (V0 + pR - pL)) * np.sinh(a * (V0 - pR + pL))))


def wavepacket(th, V0, ell, sig=150.0):
    """time-dependent lattice check: gaussian packet of left upper-band states centred in the Klein window"""
    w0 = V0 / 2
    k0 = np.arccos(np.cos(w0) / np.cos(th))
    Nx = 4096
    x = np.arange(Nx) - Nx // 2
    x0 = -int(5 * sig + 3 * ell + 100)
    ks = 2 * PI * np.fft.fftfreq(Nx)
    g = np.exp(-(ks - k0) ** 2 * sig ** 2 / 2) * np.exp(-1j * ks * x0)
    psi = np.zeros((2, Nx), complex)
    Tpred_num, norm = 0.0, 0.0
    for j, k in enumerate(ks):
        if abs(g[j]) < 1e-12:
            continue
        lo, up, ww = bands(k, th)
        psi += g[j] * up[:, None] * np.exp(1j * k * x)[None, :]
        wgt = abs(g[j]) ** 2
        Tpred_num += wgt * lattice_N(ww[1], th, V0, ell); norm += wgt
    psi /= np.sqrt((abs(psi) ** 2).sum())
    V = 0.5 * V0 * (1 + np.tanh(x / ell))
    C = coin(th); ph = np.exp(-1j * V)
    vg = np.sin(k0) * np.cos(th) / np.sin(w0)
    ticks = int((abs(x0) + 6 * sig + 4 * ell) / vg)
    ticks = min(ticks, int(0.45 * Nx / vg) + abs(x0))
    for n in range(ticks):
        R = np.roll(psi[0], 1); L = np.roll(psi[1], -1)
        psi = np.stack([ph * (C[0, 0] * R + C[0, 1] * L), ph * (C[1, 0] * R + C[1, 1] * L)])
    rho = (abs(psi) ** 2).sum(0)
    return dict(transmitted=float(rho[x > 0].sum()), predicted=float(Tpred_num / norm), ticks=ticks)


if __name__ == "__main__":
    t0 = time.time()
    res = dict()
    # ---- R1a / R1b
    scan = {}
    for th in (0.1, 0.05):
        Ec = th ** 2
        for x in (0.25, 0.5, 1, 2, 4, 8, 16, 32, 64):
            if th == 0.05 and x > 8:
                continue
            E = x * Ec
            P0 = crossing(th, E, 0.0)
            Ppi = crossing(th, E, PI)
            PS = np.exp(-PI * th ** 2 / E)
            PD = dykhne(th, E)
            scan[f"th{th}_x{x}"] = dict(th=th, E=E, x=x, P0=P0, Ppi=Ppi, P_schwinger=PS, P_dykhne=PD,
                                        Gamma_lat=E / (2 * PI) * P0, Gamma_S=E / (2 * PI) * PS,
                                        ratio=P0 / PS, ratio_dykhne=P0 / PD, pi_over_0=Ppi / P0 if P0 > 0 else None,
                                        bloch_period=2 * PI / E)
            print(th, x, {k: round(v, 5) if isinstance(v, float) else v for k, v in scan[f"th{th}_x{x}"].items()},
                  flush=True)
    res["scan"] = scan
    passx = (0.25, 0.5, 1, 2, 4)
    res["R1a_PASS"] = bool(all(abs(scan[f"th0.1_x{x}"]["ratio"] - 1) <= 0.10 for x in passx))
    res["R1b_doubler_confirmed"] = bool(all(abs(scan[f"th0.1_x{x}"]["pi_over_0"] - 1) <= 0.01 for x in passx))
    res["valid_window_th0.1"] = [x for x in (0.25, 0.5, 1, 2, 4, 8, 16, 32, 64) if abs(scan[f"th0.1_x{x}"]["ratio"] - 1) <= 0.10]
    res["valid_window_th0.05"] = [x for x in (0.25, 0.5, 1, 2, 4, 8) if abs(scan[f"th0.05_x{x}"]["ratio"] - 1) <= 0.10]
    res["real_space_max_err"] = real_space_check()
    print("R1a", res["R1a_PASS"], "R1b", res["R1b_doubler_confirmed"], "real-space err", res["real_space_max_err"],
          f"{time.time()-t0:.0f}s", flush=True)
    # ---- R1c Sauter
    th, V0 = 0.1, 0.6
    m = th
    st = {}
    for ell in (5, 10, 20, 40):
        ws = np.linspace(m, V0 - m, 121)[1:-1]
        Nl = np.array([lattice_N(w, th, V0, ell) for w in ws])
        Nc = np.array([continuum_N(w, m, V0, ell) for w in ws[::10]])
        Nf = np.array([sauter_formula(w, m, V0, ell) for w in ws])
        Ncf = Nf[::10]
        dw = ws[1] - ws[0]
        Gl = float(np.trapezoid(Nl, ws) / (2 * PI)); Gf = float(np.trapezoid(Nf, ws) / (2 * PI))
        formula_ok = bool(np.max(np.abs(Nc - Ncf)) <= 1e-4)
        Gc = Gf if formula_ok else float(np.trapezoid(Nc, ws[::10]) / (2 * PI))
        E = V0 / (2 * ell)
        st[f"ell{ell}"] = dict(ell=ell, E_over_Ecrit=E / m ** 2, Gamma_lattice=Gl, Gamma_continuum=Gc,
                               ratio=Gl / Gc, formula_vs_ode_maxdiff=float(np.max(np.abs(Nc - Ncf))),
                               formula_ok=formula_ok, N_mid_lattice=float(Nl[len(ws) // 2]),
                               N_mid_formula=float(Nf[len(ws) // 2]),
                               schwinger_local=float(np.exp(-PI * m ** 2 / E)),
                               lattice_vs_formula_maxdiff=float(np.max(np.abs(Nl - Nf))),
                               curve_w=ws[::4].tolist(), curve_lat=Nl[::4].tolist(), curve_formula=Nf[::4].tolist())
        print("ell", ell, {k: v for k, v in st[f"ell{ell}"].items() if not k.startswith("curve")}, f"{time.time()-t0:.0f}s",
              flush=True)
    res["sauter"] = st
    res["R1c_PASS"] = bool(all(abs(st[f"ell{e}"]["ratio"] - 1) <= 0.10 for e in (10, 20, 40)))
    wp = wavepacket(th, V0, 10)
    wp["ok"] = bool(abs(wp["transmitted"] - wp["predicted"]) <= 0.02)
    res["wavepacket"] = wp
    print("wavepacket", wp, flush=True)
    res["R1d_unruh"] = "skipped (time budget; declared optional)"
    res["time_s"] = time.time() - t0
    print("R1c", res["R1c_PASS"], res["time_s"])
    json.dump(res, open(os.path.join(OUT, "redo_r1.json"), "w"), indent=1, default=float)
