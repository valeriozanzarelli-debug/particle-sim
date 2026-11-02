"""emiss6_core.py - study "Emissione vera", phase 4 (docs/EMISSIONE.md §13): absorption and spontaneous emission
measured from the TIME EVOLUTION of the zig-zag rule with the TRUE charge (K = 1), long enough to see the exponential.

The atom is the radial form of the same rule (lab/emiss_radial.py: one Dirac partial wave per walk, no valleys,
any Z, small cells).  The light is the 3-D zig-zag light rule, entered through its exact linear response at the atom
(lab/emiss6_kernel.py, checked against the full box by exam G0).

  QE   one-excitation quantum light (Wigner-Weisskopf on the lattice, docs/EMISSIONE.md §4.2), with the coupling of the
       production runs (centered + transverse, lab/emiss_run.py): a pure 2p3/2 (m = 3/2) decays by itself.
  RABI absorption: the atom in its 1s under a prescribed circular field E0 (cos w t, sin w t, 0) (a coherent beam of
       the light), sigma+ from 1s m=+1/2 reaches only 2p3/2 m=+3/2.  The coupling between the two radial channels is
       the exact exponential of e E.r projected on them (unitary 2x2 per site; the other channels are a declared cut).

Units hbar = c = m = 1.  Every run is cut in pieces of <= `budget` seconds with a checkpoint (exit code 3 = again).
"""
import os, sys, time, json
import numpy as np
import numba as nb
from emiss_core import E_CHARGE, CK, Light3D
from emiss_radial import RadialAtom, VEC, _radial_tick
from zz_core import ALPHA, T_UNIT, MC2_EV

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
HC_EV_NM = 1239.841984


# ------------------------------------------------------------------ the light, through its response at the atom ----
def load_kernel(n=128, sigma=1.0, parity=False):
    """parity=False: the production coupling (one transverse projector for every deposit); True: the amended one
    (the projector of the deposit's own parity, lab/emiss6_kernel.py transverse_odd).  Returns K, (M_even, M_odd)."""
    if parity:
        d = np.load(os.path.join(CK, f"emiss6_kernelP_{n}_{sigma:g}.npz"))
        return d["K"], (d["M"], d["Mo"])
    d = np.load(os.path.join(CK, f"emiss6_kernel_{n}_{sigma:g}.npz"))
    return d["K"], (d["M"], d["M"])


class KernelLight:
    """field at the atom = sum over past deposits of K(age) . J  (two parities), / b^2; self field -M I / b^2."""

    def __init__(self, K, M, b, nmax):
        self.ages = K.shape[1] - 1
        # Kq[q, k] = kernel for a deposit made k steps before a step of parity q (deposit parity (q - k) % 2)
        self.Kq = np.zeros_like(K)
        for q in (0, 1):
            for k in range(K.shape[1]):
                self.Kq[q, k] = K[(q - k) % 2, k]
        self.Kq /= b * b
        self.Ms = M
        self.b = b
        self.J = np.zeros((nmax, 3), complex)
        self.n = 0

    def field(self, s):
        """free field at step s (after the step's tick; deposits made at steps < s)"""
        m = min(s, self.ages)
        if m == 0:
            return np.zeros(3, complex)
        q = (s + 1) % 2
        Jr = self.J[s - 1::-1][:m] if s - 1 >= 0 else None
        return np.einsum("kij,kj->i", self.Kq[q, 1:m + 1], Jr)

    def selff(self, I, s):
        """field at the atom made by the deposit of I at step s (deposit count s + 1)"""
        return -self.Ms[(s + 1) % 2] @ I / self.b**2

    def form_factor2(self, k, sigma=1.0):
        return Light3D(16, self.b, sponge=0, sigma=sigma).form_factor2(k)


# ------------------------------------------------------------------ levels of the radial rule ----------------------
_LEV = {}


def levels(Z, a, rmax):
    key = (Z, a, rmax)
    if key not in _LEV:
        at = RadialAtom(Z, a, rmax)
        oe, ce, re_ = at.eigen("e", 2)
        og, cg, rg = at.eigen("g", 1)
        ce /= np.linalg.norm(ce)
        cg /= np.linalg.norm(cg)
        R = complex(np.sum(at.r * (np.conj(cg[0]) * ce[0] + np.conj(cg[1]) * ce[1])))
        _LEV[key] = dict(at=at, ce=ce, cg=cg, oe=oe, og=og, w=oe - og, R=R, res=(re_, rg))
    return _LEV[key]


def dirac(Z):
    """exam numbers (never inputs): Dirac-Coulomb, point nucleus, infinite mass (lab/emiss_exact.py)"""
    from emiss_exact import rates
    return rates(Z * ALPHA)


# ------------------------------------------------------------------ QE: spontaneous emission ------------------------
def qe_run(name, Z, a, rmax, T_over, cells_per_lambda=30, sigma=1.0, budget=540, every=1, T_abs=None, parity=True):
    """pure 2p3/2 (m=3/2) with one-excitation light; T = T_over / A_rule (or T_abs).  Saves series to CK."""
    ck = os.path.join(CK, f"emiss6_{name}.npz")
    L = levels(Z, a, rmax)
    at, w, R = L["at"], L["w"], L["R"]
    lam = 2 * np.pi / w
    sub = max(1, int(round(lam / cells_per_lambda / a)))
    b = sub * a
    e = E_CHARGE
    K, M = load_kernel(128, sigma, parity)
    W2 = Light3D(16, b, sponge=0, sigma=sigma).form_factor2(w)
    A_rule = 4 / 3 * ALPHA * w**3 * abs(R)**2 / 3
    T = T_abs if T_abs is not None else T_over / (A_rule * W2)
    nsteps = int(np.ceil(T / b))
    KL = KernelLight(K, M, b, nsteps + 2)
    meta = dict(name=name, Z=Z, a=a, rmax=rmax, N=at.r.size, sub=sub, b=b, omega=w, R=abs(R), A_rule=A_rule, W2=W2,
                T=T, nsteps=nsteps, cells_per_lambda=cells_per_lambda, sigma=sigma, ages=KL.ages, parity=parity)
    t0 = time.time()
    if os.path.exists(ck):
        d = np.load(ck)
        psi, phi, s0 = d["psi"].copy(), d["phi"].copy(), int(d["s"])
        KL.J[:s0] = d["J"][:s0]
        ser = {k: list(d["ser_" + k]) for k in ("t", "norm", "D", "p2p")}
        print(f"[{name}] resumed at step {s0}/{nsteps}", flush=True)
    else:
        psi, phi, s0 = L["ce"].copy(), L["cg"].copy(), 0
        ser = dict(t=[0.0], norm=[1.0], D=[complex(np.sum(at.r * (np.conj(phi[0]) * psi[0] + np.conj(phi[1]) * psi[1])))],
                   p2p=[1.0])
    ce = L["ce"]
    for s in range(s0, nsteps):
        D0 = at.dip_matrix(phi, psi)
        h = KL.J[:s]
        if s >= 2 and np.vdot(h[-2], h[-2]).real > 0:
            r = np.vdot(h[-2], h[-1]) / np.vdot(h[-2], h[-2])
            r = r / max(1.0, abs(r))
            jp = h[-1] * r
        elif s == 1:
            jp = h[-1]
        else:
            jp = np.zeros(3, complex)
        field = KL.field(s) + 0.5 * KL.selff(jp, s)
        at.qe_block(psi, [phi], [field], sub, a, e)
        D1 = at.dip_matrix(phi, psi)
        J = -e * (D1 - D0) / b
        dE = 0.5 * KL.selff(J - jp, s)
        at.add_rphi(psi, -1j * b * e * dE, phi)
        KL.J[s] = J
        if (s + 1) % every == 0 or s == nsteps - 1:
            ser["t"].append((s + 1) * b)
            ser["norm"].append(float(np.sum(np.abs(psi)**2)))
            ser["D"].append(complex(np.dot(D1, np.conj(VEC)) / np.dot(VEC, np.conj(VEC)).real))
            # 2p population: the level evolves as e^{-i oe t}; |overlap|^2 is phase free
            ser["p2p"].append(float(abs(np.vdot(ce, psi))**2))
        if time.time() - t0 > budget and s < nsteps - 1:
            np.savez(ck + ".tmp.npz", psi=psi, phi=phi, s=s + 1, J=KL.J[:s + 1], **{"ser_" + k: np.array(v) for k, v in ser.items()})
            os.replace(ck + ".tmp.npz", ck)
            print(f"[{name}] checkpoint at step {s + 1}/{nsteps} (t = {(s + 1) * b:.4g}), norm {ser['norm'][-1]:.6f}", flush=True)
            return None, meta
    res = dict(meta, t=np.array(ser["t"]), norm=np.array(ser["norm"]), D=np.array(ser["D"]), p2p=np.array(ser["p2p"]),
               J=KL.J[:nsteps].copy())
    np.savez(ck + ".tmp.npz", psi=psi, phi=phi, s=nsteps, J=KL.J[:nsteps], **{"ser_" + k: np.array(v) for k, v in ser.items()})
    os.replace(ck + ".tmp.npz", ck)
    print(f"[{name}] done ({time.time() - t0:.0f} s)", flush=True)
    return res, meta


# ------------------------------------------------------------------ RABI: absorption -------------------------------
@nb.njit(cache=True)
def _rabi_block(pe, pg, Me, r0e, r1e, Mg, r0g, r1g, base, w_d, t0, a, nt, buf, rec_every, out_t, out_pe, out_pg,
                ce, cg, irec, tau, tc):
    """nt ticks: both walks tick, then the exact 2x2 exponential of the channel coupling (circular field).
    base[i] = a e E0 r_i / sqrt6; envelope 1 (tau <= 0) or exp(-((t - tc)/tau)^2)."""
    N = pe.shape[1]
    for s in range(nt):
        _radial_tick(pe, Me, r0e, r1e, buf)
        _radial_tick(pg, Mg, r0g, r1g, buf)
        t = t0 + (s + 0.5) * a
        env = 1.0
        if tau > 0:
            x = (t - tc) / tau
            env = np.exp(-x * x)
        # E.conj(VEC) = -E0 env e^{-i w t}/sqrt6:  U = cos(g) + i sin(g) [[0, u],[u*, 0]],  g = base env
        u = np.exp(-1j * w_d * t)
        uc = np.conj(u)
        for i in range(N):
            g = base[i] * env
            C = np.cos(g)
            S = np.sin(g)
            for k in range(2):
                x0 = pe[k, i]
                y0 = pg[k, i]
                pe[k, i] = C * x0 + 1j * S * u * y0
                pg[k, i] = C * y0 + 1j * S * uc * x0
        if (s + 1) % rec_every == 0:
            j = irec[0]
            if j < out_t.size:
                out_t[j] = t0 + (s + 1) * a
                ae = 0j
                ag = 0j
                for i in range(N):
                    for k in range(2):
                        ae += np.conj(ce[k, i]) * pe[k, i]
                        ag += np.conj(cg[k, i]) * pg[k, i]
                out_pe[j] = ae.real * ae.real + ae.imag * ae.imag
                out_pg[j] = ag.real * ag.real + ag.imag * ag.imag
                irec[0] = j + 1


def rabi_run(name, Z, a, rmax, omega_frac=0.025, detune=0.0, cycles=1.3, budget=540, nrec=4000, pulse_tau=None,
             p_target=1e-4):
    """1s m=+1/2 under E0 env(t) (cos, sin, 0) at w_d = w_rule + detune * Omega_rule.
    cw (pulse_tau None): Omega_rule = omega_frac * w, run `cycles` generalized Rabi cycles.
    weak pulse (pulse_tau = tau in optical periods): env = exp(-((t - 3 tau)/tau)^2), T = 6 tau, E0 chosen so that
    first-order theory with the rule's own R gives P_e = p_target; R_dyn = sqrt(6 P_e / pi) / (e E0 tau)."""
    ck = os.path.join(CK, f"emiss6_{name}.npz")
    L = levels(Z, a, rmax)
    at, w, R = L["at"], L["w"], abs(L["R"])
    if pulse_tau is None:
        Om_rule = omega_frac * w                   # Rabi frequency 2V, V = e E0 R / sqrt6 the coupling element
        eE0 = 0.5 * Om_rule * np.sqrt(6) / R       # chosen so that the rule's own Rabi frequency is omega_frac * w
        wd = w + detune * Om_rule
        Om_gen = np.sqrt(Om_rule**2 + (detune * Om_rule)**2)
        T = cycles * 2 * np.pi / Om_gen
        tau, tc = -1.0, 0.0
    else:
        tau = pulse_tau * 2 * np.pi / w
        tc = 3 * tau
        T = 6 * tau
        eE0 = np.sqrt(6 * p_target / np.pi) / (R * tau)
        Om_rule = 2 * eE0 * R / np.sqrt(6)
        wd = w + detune / tau
    nt = int(np.ceil(T / a))
    base = at.r * a * eE0 / np.sqrt(6)             # a e |E.conj(VEC)| r per site
    meta = dict(name=name, Z=Z, a=a, rmax=rmax, N=at.r.size, omega=w, R=R, eE0=eE0, Om_rule=Om_rule, detune=detune,
                omega_drive=wd, T=T, nticks=nt, tau=tau, p_target=p_target if pulse_tau else None)
    rec_every = max(1, nt // nrec)
    nr = nt // rec_every + 2
    if os.path.exists(ck):
        d = np.load(ck)
        pe, pg, s0 = d["pe"].copy(), d["pg"].copy(), int(d["s"])
        out_t, out_pe, out_pg, irec = d["out_t"].copy(), d["out_pe"].copy(), d["out_pg"].copy(), d["irec"].copy()
        print(f"[{name}] resumed at tick {s0}/{nt}", flush=True)
    else:
        pe = np.zeros_like(L["ce"]); pg = L["cg"].copy(); s0 = 0
        out_t, out_pe, out_pg = np.zeros(nr), np.zeros(nr), np.zeros(nr)
        irec = np.zeros(1, np.int64)
    buf = np.empty((2, at.r.size), np.complex128)
    we, wg = at.we, at.wg
    t0w = time.time()
    chunk = max(rec_every, (2_000_000_000 // max(1, at.r.size * 40)) // rec_every * rec_every)
    s = s0
    while s < nt:
        n = min(chunk, nt - s)
        n = max(rec_every, n // rec_every * rec_every) if nt - s >= rec_every else nt - s
        _rabi_block(pe, pg, we.M, we.refl0, we.refl1, wg.M, wg.refl0, wg.refl1, base, wd, s * a, a, n, buf, rec_every,
                    out_t, out_pe, out_pg, L["ce"], L["cg"], irec, tau, tc)
        s += n
        if time.time() - t0w > budget and s < nt:
            np.savez(ck + ".tmp.npz", pe=pe, pg=pg, s=s, out_t=out_t, out_pe=out_pe, out_pg=out_pg, irec=irec)
            os.replace(ck + ".tmp.npz", ck)
            print(f"[{name}] checkpoint at tick {s}/{nt}", flush=True)
            return None, meta
    m = int(irec[0])
    np.savez(ck + ".tmp.npz", pe=pe, pg=pg, s=s, out_t=out_t, out_pe=out_pe, out_pg=out_pg, irec=irec)
    os.replace(ck + ".tmp.npz", ck)
    res = dict(meta, t=out_t[:m], pe=out_pe[:m], pg=out_pg[:m], norm_end=float(np.sum(np.abs(pe)**2) + np.sum(np.abs(pg)**2)))
    print(f"[{name}] done ({time.time() - t0w:.0f} s)", flush=True)
    return res, meta


def fit_rabi(t, pe):
    """P_e = Pmax sin^2(W t / 2): W by a grid search of the least squares (Pmax solved linearly), then refined."""
    from scipy.optimize import curve_fit
    f = lambda tt, P, W: P * np.sin(W * tt / 2)**2
    Ws = np.linspace(0.2, 20, 4000) * 2 * np.pi / t[-1]
    best = None
    for W in Ws:
        s = np.sin(W * t / 2)**2
        P = np.dot(s, pe) / np.dot(s, s)
        r = np.sum((P * s - pe)**2)
        if best is None or r < best[0]:
            best = (r, P, W)
    p, _ = curve_fit(f, t, pe, p0=[best[1], best[2]])
    resid = float(np.max(np.abs(f(t, *p) - pe)))
    return float(p[0]), float(p[1]), resid


# ------------------------------------------------------------------ ABSORPTION OF A LIGHT-WALK PACKET ----------------
def packet_field(hel=+1, k_cells=2 * np.pi / 30, s_cells=15.0, z0=-35, steps=240, n=128, sigma=1.0):
    """A photon packet MADE OF THE LIGHT RULE: F = E + iB = g(z) e^{i k z} (1, i hel, 0)/sqrt2 in the 128^3 box (b = 1),
    g Gaussian of width s_cells centred z0 cells from the atom; the light rule (walk + sponge) carries it.  Returns the
    Gaussian-averaged E = Re F the atom would feel, per light step (step s = after tick s+1), shape (steps, 3); the
    free packet does not depend on Z (the light rule has no scale): with light cell b the time of step s is (s+1) b."""
    f = os.path.join(CK, f"emiss6_packet_{hel:+d}_{k_cells:.5f}_{s_cells:g}_{z0}_{steps}.npy")
    if os.path.exists(f):
        return np.load(f)
    L = Light3D(n, 1.0, sponge=12, sigma=sigma)
    z = np.arange(n) - n // 2
    g = np.exp(-(z - z0)**2 / (2 * s_cells**2)) * np.exp(1j * k_cells * z)
    L.F[0] = g[None, None, :] / np.sqrt(2)
    L.F[1] = 1j * hel * g[None, None, :] / np.sqrt(2)
    L.F[2] = 0
    E = np.zeros((steps, 3))
    for s in range(steps):
        L.tick()
        E[s] = L.field_at_atom()
        L.sponge()
    np.save(f, E)
    return E


@nb.njit(cache=True)
def _drive_block(pe, pg, Me, r0e, r1e, Mg, r0g, r1g, r, a, sub, w_steps, buf):
    """atom driven by a field given per light step: w = e E . conj(VEC) (complex), constant over the sub ticks of a
    step; exact 2x2 exponential of -i a r [[0, w], [w*, 0]] per site after each tick."""
    N = pe.shape[1]
    for st in range(w_steps.size):
        w = w_steps[st]
        aw = abs(w)
        ph = w / aw if aw > 0 else 1.0 + 0j
        for s in range(sub):
            _radial_tick(pe, Me, r0e, r1e, buf)
            _radial_tick(pg, Mg, r0g, r1g, buf)
            if aw == 0:
                continue
            for i in range(N):
                g = a * r[i] * aw
                C = np.cos(g)
                S = np.sin(g)
                for k in range(2):
                    x0 = pe[k, i]
                    y0 = pg[k, i]
                    pe[k, i] = C * x0 - 1j * S * ph * y0
                    pg[k, i] = C * y0 - 1j * S * np.conj(ph) * x0


def packet_absorb(Z, a, rmax, hel=+1, p_target=1e-4, cells_per_lambda=30, steps=240):
    """the 1s meets the light-walk packet; P_e = |<2p|psi>|^2 after it passed.  First order:
    P_e = e^2 R^2 |E_c(w)|^2, E_c(w) = sum_s b E(s).conj(VEC) e^{i w t_s}  ->  R_dyn = sqrt(P_e) / (e |E_c(w)|)."""
    L = levels(Z, a, rmax)
    at, w, R = L["at"], L["w"], abs(L["R"])
    lam = 2 * np.pi / w
    sub = max(1, int(round(lam / cells_per_lambda / a)))
    b = sub * a
    E1 = packet_field(hel, 2 * np.pi * b / lam, 15.0, -35, steps)       # k in cells^-1 matched to this b
    ts = (np.arange(steps) + 0.5) * b                                    # mid-step time (field held over the step)
    Ec = E1 @ np.conj(VEC)
    # spectral amplitude at the line with the same piecewise-constant hold the walk sees
    Ew = np.sum(Ec * np.exp(1j * w * ts)) * b * np.sinc(w * b / (2 * np.pi))
    scale = np.sqrt(p_target) / (E_CHARGE * R * abs(Ew))
    wst = (E_CHARGE * scale * Ec).astype(np.complex128)
    pe = np.zeros_like(L["ce"]); pg = L["cg"].copy()
    buf = np.empty((2, at.r.size), np.complex128)
    t0 = time.time()
    _drive_block(pe, pg, at.we.M, at.we.refl0, at.we.refl1, at.wg.M, at.wg.refl0, at.wg.refl1, at.r, a, sub, wst, buf)
    Pe = float(abs(np.vdot(L["ce"], pe))**2)
    R_dyn = np.sqrt(Pe) / (E_CHARGE * scale * abs(Ew))
    Eneg = abs(np.sum(Ec * np.exp(-1j * w * ts)))
    return dict(Z=Z, a=a, rmax=rmax, N=at.r.size, sub=sub, b=b, omega=w, R_rule=R, R_dyn=float(R_dyn), P_e=Pe,
                hel=hel, packet=dict(lambda_cells=cells_per_lambda, k_cells=2 * np.pi * b / lam, width_cells=15.0,
                                     start_cells=-35, steps=steps, amplitude_F=float(scale)),
                ratio_pos_neg=float(abs(np.sum(Ec * np.exp(1j * w * ts))) / max(Eneg, 1e-300)),
                f32=float(2 / 3 * w * (4 / 6) * R_dyn**2), f_tot_nr=float(2 / 3 * w * R_dyn**2),
                seconds=time.time() - t0)
