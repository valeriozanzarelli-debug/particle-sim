"""IDENTICAL BALLS — a photon made of N identical balls, every ball follows the SAME local rule.

Paraxial 2-D model: x = transverse coordinate, z = distance along the beam (plays the role of time).
The photon's transverse motion is that of a particle of 'mass' k = 2π/λ with ħ = 1:
        i ∂ψ/∂z = −(1/2k) ∂²ψ/∂x²
A ball moves along z at c; the rule decides dx/dz. The balls are ordered along x, so "my neighbours"
= the m balls on my left and the m balls on my right (their rank offsets −m..m are known locally).
Near the tails the window slides inward (the outermost ball uses its 2m inner neighbours).

Rules (all balls identical; a rule may only read the neighbours' positions and what they carry):
  A  HDW  'Many Interacting Worlds' force (Hall-Deckert-Wiseman 2014): U = Σ (1/8k)(1/g_n − 1/g_{n−1})²,
          g = gap to the neighbour; force on a ball involves its 2 neighbours on each side (5 points).
  B  SG   density from a polynomial fit of the neighbours' positions vs their rank (Savitzky-Golay in
          rank space): J = dx/dj  (ρ = 1/(N J)), quantum potential Q = −(1/2k)(√ρ)''/√ρ and its gradient
          from J, J', J'', J''' — all from the same fit.  a = −Q'/k.
  BH       Hamiltonian smoothing of HDW: U = (1/8k) Σ J'²/J⁴ with J, J' from a rank fit over 2m+1 balls.
  P  POLAR each ball also carries amplitude R and phase S: v = S'/k (local least squares of the neighbours'
          phases), dS/dz = k v²/2 − Q, Q = −(1/2k) R''/R (local fit of the neighbours' R), dR/dz = −R v'/2.
  C  PSI  each ball carries a complex amplitude ψ (a Jones vector for polarization). Each tick it reads
          its 2m nearest balls (by distance) and
            1) fits log(ψ_j/ψ_i) ≈ a1 d + a2 d²  (a local Gaussian-beamlet 'carrier' E = e^{a1 d + a2 d²}),
            2) fits the rest ψ_j/E_j − ψ_i ≈ c1 d + c2 d²,
            3) ψ' = a1 ψ + c1,  ψ'' = (2a2 + a1²) ψ + 2 a1 c1 + 2 c2,
               velocity v = Im(ψ*·ψ')/(k|ψ|²)                    (it moves with the local current)
               dψ/dz = (i/2k) ψ'' + v ψ'                           (Schrödinger seen from the moving ball)
          The balls never read a density: the ball density follows |ψ|² because every ball moves with the
          local current (equivariance). Integrator 'cn2': trapezoid (Crank-Nicolson) with the same ball speed
          in the motion and in both operators.  'herm' = unitary nearest-neighbour links + move-and-reread
          (tried; worse).

Everything here is vectorized over the N balls of ONE photon; nothing global enters the rule
(the only global thing is the common z-step of the integrator, see the implicit-step note in the report).
"""
import sys, time, numpy as np
from scipy.special import ndtri, factorial
from scipy.linalg import solve_banded
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from optics import (K, LAM, W0, D, A, L, SLITS, sources, rs_field, flux, screen_reference, coarse_corr,
                    fringe_errors, ideal_baseline, in_slit, XS, EDGES, HALF, BINS)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

W0F = 3.0                       # free Gaussian beam waist (intensity ∝ exp(−2x²/w0²), σ = w0/2)
ZR = np.pi * W0F ** 2 / LAM
ZF = (ZR, 2 * ZR, 5 * ZR)


def w_theory(z):
    return W0F * np.sqrt(1 + (z / ZR) ** 2)


# ───────────────────────────────── neighbour windows ─────────────────────────────────
_WIN = {}


def windows(N, m):
    """Rank window of every ball: ranks lo..lo+2m (slides inward at the tails). Returns idx (N,2m+1), self column."""
    key = (N, m)
    if key not in _WIN:
        lo = np.clip(np.arange(N) - m, 0, max(N - 2 * m - 1, 0))
        idx = lo[:, None] + np.arange(2 * m + 1)[None, :]
        _WIN[key] = (idx, np.arange(N) - lo)
    return _WIN[key]


_SG = {}


def sg_weights(m, P):
    """Savitzky-Golay derivative weights in rank space for every self position t0 = 0..2m:
    W[t0, p, j] -> p-th derivative (d/dj)^p at j = t0 of the degree-P fit through the 2m+1 values."""
    key = (m, P)
    if key not in _SG:
        j = np.arange(2 * m + 1)
        W = np.empty((2 * m + 1, P + 1, 2 * m + 1))
        for t0 in range(2 * m + 1):
            V = (j[:, None] - t0) ** np.arange(P + 1)[None, :]
            W[t0] = np.linalg.pinv(V) * factorial(np.arange(P + 1))[:, None]
        _SG[key] = W
    return _SG[key]


NEIGH = {"mode": "knn"}


def neighbours(x, m):
    """The ball's 2m+1-window (itself included), contiguous in rank.
    'rank': centred on the ball (slides inward at the ends).
    'knn' : the 2m nearest balls by DISTANCE (a ball at the edge of a group listens only to its own group
            until the other group comes closer than its own 2m-th neighbour)."""
    N = x.size
    idx, t0 = windows(N, m)
    if NEIGH["mode"] == "rank":
        return idx, t0
    i = np.arange(N)
    cand = np.clip(i[:, None] - np.arange(2 * m + 1)[None, :], 0, N - 2 * m - 1)     # lo candidates (N, 2m+1)
    reach = np.maximum(x[:, None] - x[cand], x[cand + 2 * m] - x[:, None])
    lo = np.take_along_axis(cand, np.argmin(reach, 1)[:, None], 1)[:, 0]
    return lo[:, None] + np.arange(2 * m + 1)[None, :], i - lo


def xfit(x, m, P):
    """Local weighted least squares in x-space over the rank window (self excluded, polynomial without the
    constant term, i.e. anchored at the ball's own value). Returns d (N,2m+1) and G (N,P,2m+1) such that
    coefficient c_p = Σ_j G[:,p-1,j] (f_j − f_self) of the fit f ≈ f_self + Σ c_p d^p."""
    N = x.size
    idx, t0 = neighbours(x, m)
    d = x[idx] - x[:, None]
    Dm = np.abs(d).max(1) * 1.05 + 1e-300
    t = d / Dm[:, None]
    w = (1 - t ** 2) ** 2
    w[np.arange(N), t0] = 0.0
    Apow = t[..., None] ** np.arange(1, P + 1)                 # (N, 2m+1, P)
    Aw = Apow * w[..., None]
    M = np.einsum("njp,njq->npq", Aw, Apow)
    M += 1e-12 * np.eye(P)[None]
    G = np.linalg.solve(M, Aw.transpose(0, 2, 1))              # (N, P, 2m+1)
    G /= (Dm[:, None] ** np.arange(1, P + 1))[:, :, None]
    xfit.w = w
    return idx, d, G


# ───────────────────────────────── rule A: HDW ───────────────────────────────────────
def acc_hdw(x):
    N = x.size
    a = np.zeros(N + 1)
    a[1:N] = 1.0 / np.diff(x)
    g = np.zeros(N + 1)
    g[1:N] = (2 * a[1:N] - a[0:N - 1] - a[2:N + 1]) / (4 * K)      # dU/da_n, U = (1/8k) Σ (a_{n+1} − a_n)^2
    dU = g[1:] * a[1:] ** 2 - g[:-1] * a[:-1] ** 2
    return -dU / K


# ───────────────────────────────── rule B: rank-space SG ─────────────────────────────
def acc_sg(x, m=4, P=4):
    N = x.size
    idx, t0 = windows(N, m)
    W = sg_weights(m, P)[t0]                                   # (N, P+1, 2m+1)
    dd = np.einsum("npj,nj->np", W, x[idx])
    J, J1, J2, J3 = dd[:, 1], dd[:, 2], dd[:, 3], dd[:, 4]
    J = np.maximum(J, 1e-12)
    dQ = -(1 / (2 * K)) * (4 * J1 * J2 / J ** 4 - 5 * J1 ** 3 / J ** 5 - 0.5 * J3 / J ** 3)
    return -dQ / (K * J)


def acc_sgH(x, m=3, P=2):
    """Rule BH: Hamiltonian smoothing of HDW. Every ball i reads J_i = dx/dj and J'_i = d²x/dj² from a
    degree-P rank fit over its 2m+1 neighbours; the quantum energy is the rank-space Fisher information
        U = (1/8k) Σ_i J'_i² / J_i⁴        (HDW is the special case J = gap, J' = gap difference)
    and the force is −∂U/∂x_l (each ball collects the pushes of the neighbours whose window contains it)."""
    N = x.size
    idx, t0 = windows(N, m)
    W = sg_weights(m, P)[t0]
    xx = x[idx]
    J = np.maximum(np.einsum("nj,nj->n", W[:, 1], xx), 1e-12)
    J1 = np.einsum("nj,nj->n", W[:, 2], xx)
    dUdJ = -(4.0 / (8 * K)) * J1 ** 2 / J ** 5
    dUdJ1 = (2.0 / (8 * K)) * J1 / J ** 4
    contrib = dUdJ[:, None] * W[:, 1] + dUdJ1[:, None] * W[:, 2]      # ∂U_i/∂x_{idx[i,j]}
    grad = np.bincount(idx.ravel(), contrib.ravel(), N)
    return -grad / K


# ───────────────────────────────── explicit Verlet (rules A, B) ─────────────────────
def verlet(x, v, zs, acc, dz_fn, z0=0.0, max_steps=3_000_000, tmax=500, log=None):
    """Velocity Verlet, common step dz = dz_fn(x). Records x at the requested z's.
    Returns snapshots, dict(steps, swaps, time)."""
    x = x.copy(); v = v.copy(); z = z0
    a = acc(x); out = []; steps = swaps = 0; t = time.time()
    for zt in zs:
        while z < zt - 1e-12 and steps < max_steps:
            h = min(dz_fn(x), zt - z)
            v += 0.5 * h * a
            x += h * v
            if np.any(np.diff(x) <= 0):
                o = np.argsort(x, kind="stable"); swaps += int((o != np.arange(x.size)).sum())
                x = x[o]; v = v[o]
            a = acc(x)
            v += 0.5 * h * a
            z += h; steps += 1
            if not np.all(np.isfinite(x)) or time.time() - t > tmax:
                return out + [None] * (len(zs) - len(out)), dict(steps=steps, swaps=swaps, time=time.time() - t, z=z, fail=True)
            if log and steps % log == 0:
                print(f"    z={z:.3f} steps={steps} dz={h:.2e} width={2*x.std():.3f}", flush=True)
        out.append(x.copy())
    return out, dict(steps=steps, swaps=swaps, time=time.time() - t, z=z, fail=steps >= max_steps)


def implicit_midpoint(x, v, zs, acc, band, dz_fn, z0=0.0, newton=10, tol=1e-9, tmax=500):
    """Implicit midpoint (symplectic, A-stable) for x'' = acc(x); Jacobian of acc is banded (half-width
    `band`), obtained with 2·band+1 coloured finite differences; banded Newton solve. A step whose Newton
    iteration does not converge (or that would make two balls cross) is retried with half the step."""
    N = x.size; x = x.copy(); v = v.copy(); z = z0; out = []; steps = swaps = 0; t = time.time()
    nb = 2 * band + 1
    hcur = None

    def try_step(h):
        x1 = x + h * v + 0.5 * h * h * acc(x)
        if np.any(np.diff(x1) <= 0):
            x1 = x + h * v
        for it in range(newton):
            mid = 0.5 * (x + x1)
            a0 = acc(mid)
            Gv = x1 - x - h * v - 0.5 * h * h * a0
            eps = 1e-5 * np.maximum(np.abs(np.diff(mid)).min(), 1e-12)
            ab = np.zeros((nb, N))
            for c in range(nb):
                pert = np.zeros(N); pert[c::nb] = eps
                da = (acc(mid + pert) - a0) / eps                 # column c::nb contributions
                cols = np.arange(c, N, nb)
                for off in range(-band, band + 1):
                    rows = cols + off
                    ok = (rows >= 0) & (rows < N)
                    ab[band + rows[ok] - cols[ok], cols[ok]] = da[rows[ok]]
            ab *= -0.25 * h * h                                      # dG/dx1 = I − (h²/4) J
            ab[band] += 1.0
            try:
                dx = solve_banded((band, band), ab, Gv)
            except (ValueError, np.linalg.LinAlgError):
                return None
            x1 = x1 - dx
            if not np.all(np.isfinite(x1)):
                return None
            if np.abs(dx).max() < tol * (1 + np.abs(x1).max()):
                return x1 if np.all(np.diff(x1) > 0) else None
        return None

    for zt in zs:
        while z < zt - 1e-12:
            hmax = min(dz_fn(x), zt - z)
            h = hmax if hcur is None else min(hmax, 1.5 * hcur)
            x1 = try_step(h)
            while x1 is None and h > 1e-7:
                h *= 0.5
                x1 = try_step(h)
            if x1 is None or time.time() - t > tmax:
                return out + [None] * (len(zs) - len(out)), dict(steps=steps, swaps=swaps, time=time.time() - t, z=z, fail=True)
            hcur = h
            v = 2 * (x1 - x) / h - v
            x = x1; z += h; steps += 1
        out.append(x.copy())
    return out, dict(steps=steps, swaps=swaps, time=time.time() - t, z=z, fail=False)


# ───────────────────────────────── rule P: polar (R, S carried) ─────────────────────
def rhs_polar(x, R, S, m, P):
    idx, d, G = xfit(x, m, P)
    dS = S[idx] - S[:, None]
    dR = R[idx] - R[:, None]
    S1 = np.einsum("nj,nj->n", G[:, 0], dS); S2 = 2 * np.einsum("nj,nj->n", G[:, 1], dS)
    R2 = 2 * np.einsum("nj,nj->n", G[:, 1], dR)
    v = S1 / K
    Q = -(1 / (2 * K)) * R2 / np.where(np.abs(R) > 1e-300, R, 1e-300)
    return v, -0.5 * R * S2 / K, K * v * v / 2 - Q


def run_polar(x, R, S, zs, m, P, dz_fn, z0=0.0, tmax=500):
    x = x.copy(); R = R.copy(); S = S.copy(); z = z0; out = []; steps = swaps = 0; t = time.time()
    for zt in zs:
        while z < zt - 1e-12:
            h = min(dz_fn(x), zt - z)
            k1 = rhs_polar(x, R, S, m, P)
            k2 = rhs_polar(x + .5 * h * k1[0], R + .5 * h * k1[1], S + .5 * h * k1[2], m, P)
            k3 = rhs_polar(x + .5 * h * k2[0], R + .5 * h * k2[1], S + .5 * h * k2[2], m, P)
            k4 = rhs_polar(x + h * k3[0], R + h * k3[1], S + h * k3[2], m, P)
            x = x + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            R = R + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            S = S + h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
            z += h; steps += 1
            if np.any(np.diff(x) <= 0):
                o = np.argsort(x, kind="stable"); swaps += int((o != np.arange(x.size)).sum())
                x, R, S = x[o], R[o], S[o]
            if not np.all(np.isfinite(x)) or time.time() - t > tmax:
                return out + [None] * (len(zs) - len(out)), dict(steps=steps, swaps=swaps, time=time.time() - t, z=z, fail=True)
        out.append(x.copy())
    return out, dict(steps=steps, swaps=swaps, time=time.time() - t, z=z, fail=False)


# ───────────────────────────────── rule C: ψ carried ────────────────────────────────
DEMOD = {"mode": "log", "clip": 3.0}


def psi_fit(x, psi, v, m, P):
    if DEMOD["mode"] == "log":
        return psi_fit_log(x, psi, m, P, DEMOD["clip"])
    return psi_fit_phase(x, psi, v, m, P)


def psi_fit_log(x, psi, m, P, clip=3.0, v_move=None, parts=False):
    """One read of the neighbours, 'log carrier' version.
    1) the ball compares its ψ with each neighbour's: ℓ_j = log(ψ_j/ψ_i) (for a Jones vector: projection
       on its own ψ), fits ℓ ≈ a1 d + a2 d²  -> a local Gaussian-beamlet carrier E(d) = exp(a1 d + a2 d²)
       (amplitude part clipped to |Re| ≤ clip over the window);
    2) fits the residual ψ_j/E_j − ψ_i with a polynomial c1 d + … + cP d^P;
    3) ψ' = a1 ψ + c1,  ψ'' = (2a2 + a1²) ψ + 2 a1 c1 + 2 c2,  v = Im(ψ*·ψ')/(k|ψ|²),
       dψ/dz = (i/2k) ψ'' + u ψ'   (Schrödinger seen from a ball moving with speed u; u = v normally).
    Returns v and the rows of the (frozen-carrier) linear operator dψ/dz = L ψ (for u = v_move or v)."""
    idx, d, G = xfit(x, m, P)
    n2 = (np.abs(psi) ** 2).sum(-1)
    proj = np.einsum("njc,nc->nj", psi[idx], np.conj(psi)) / np.maximum(n2, 1e-300)[:, None]   # ψ_j·ψ_i*/|ψ_i|²
    proj = np.where(np.abs(proj) > 1e-30, proj, 1e-30)
    ph = np.unwrap(np.angle(proj), axis=1)
    ph -= np.take_along_axis(ph, np.argmin(np.abs(d), 1)[:, None], 1)     # own phase = 0
    ell = np.log(np.abs(proj)) + 1j * ph
    G2q = G if P == 2 else xfit(x, m, 2)[2]
    wq = xfit.w
    a1 = np.einsum("nj,nj->n", G2q[:, 0], ell); a2 = np.einsum("nj,nj->n", G2q[:, 1], ell)
    # does the Gaussian-beamlet carrier describe my neighbours? (weighted rms misfit of the log-ratio)
    mis = np.sqrt((wq * np.abs(ell - a1[:, None] * d - a2[:, None] * d * d) ** 2).sum(1) / wq.sum(1))
    good = mis < DEMOD.get("tol", 1e9)
    a1 = np.where(good, a1, 0); a2 = np.where(good, a2, 0)
    psi_fit_log.frac_off = 1 - good.mean()
    Dm = np.abs(d).max(1) + 1e-300
    a1 = np.clip(a1.real, -clip / Dm, clip / Dm) + 1j * a1.imag
    a2 = np.clip(a2.real, -clip / Dm ** 2, clip / Dm ** 2) + 1j * a2.imag
    if DEMOD.get("vmax"):                                            # a phase slope steeper than k is evanescent
        a1 = a1.real + 1j * np.clip(a1.imag, -K * DEMOD["vmax"], K * DEMOD["vmax"])
    if CARRIER == "phase":                                          # keep only the phase part of the carrier
        a1 = 1j * a1.imag; a2 = 1j * a2.imag
    elif CARRIER == "none":
        a1 = 0 * a1; a2 = 0 * a2
    E = np.exp(a1[:, None] * d + a2[:, None] * d * d)
    y = psi[idx] / E[..., None] - psi[:, None, :]
    c1 = np.einsum("nj,njc->nc", G[:, 0], y)
    dpsi = a1[:, None] * psi + c1
    v = np.imag(np.conj(psi) * dpsi).sum(-1) / (K * np.maximum(n2, 1e-300))
    if DEMOD.get("vmax"):                                            # no ball moves sideways faster than light
        v = np.clip(v, -DEMOD["vmax"], DEMOD["vmax"])
    if DEMOD.get("vsmooth"):
        # the ball moves with the straight line that best fits its own and its neighbours' current velocities
        # (a linear flow, e.g. a spreading beam, is kept exactly; ball-to-ball jitter is removed)
        ww = wq.copy(); ww[np.arange(x.size), np.argmin(np.abs(d), 1)] = 1.0
        vj = v[idx]
        S0 = ww.sum(1); S1 = (ww * d).sum(1); S2 = (ww * d * d).sum(1)
        T0 = (ww * vj).sum(1); T1 = (ww * vj * d).sum(1)
        v = (S2 * T0 - S1 * T1) / (S0 * S2 - S1 * S1)
    Acoef = (1j / K) * (a1[:, None] * G[:, 0] + G[:, 1])
    pr = (idx, E, Acoef, G[:, 0], a1, a2)
    psi_fit_log.last = (a1, a2, c1, np.einsum("nj,njc->nc", G[:, 1], y))
    if parts:
        return v, pr
    u = v if v_move is None else v_move
    if not DEMOD.get("move", True):
        u = 0 * v                                                   # diagnostic: balls frozen (Eulerian)
    return (v,) + rows_from_parts(pr, u)


CARRIER = "full"


def rows_from_parts(pr, u):
    idx, E, Acoef, G1, a1, a2 = pr
    coef = Acoef + u[:, None] * G1                                  # weight on (ψ_j/E_j − ψ_i)
    off = coef / E
    diag = (1j / (2 * K)) * (2 * a2 + a1 * a1) + u * a1 - coef.sum(1)
    return idx, off, diag


def psi_fit_phase(x, psi, v, m, P):
    """One read of the neighbours. Returns Bohm velocity, the operator rows (idx, off-diagonal, diagonal)
    of dψ/dz = i k v²/2 ψ + (i/k) c2."""
    idx, d, G = xfit(x, m, P)
    ph = np.exp(-1j * K * v[:, None] * d)                       # remove own phase slope
    y = psi[idx] * ph[..., None] - psi[:, None, :]              # (N, 2m+1, nc)
    c1 = np.einsum("nj,njc->nc", G[:, 0], y)
    n2 = (np.abs(psi) ** 2).sum(-1)
    vB = v + (np.imag(np.conj(psi) * c1).sum(-1)) / (K * np.maximum(n2, 1e-300))
    # rows with the NEW velocity
    ph = np.exp(-1j * K * vB[:, None] * d)
    off = (1j / K) * G[:, 1] * ph                               # coefficient on ψ_j
    diag = 1j * K * vB ** 2 / 2 - (1j / K) * G[:, 1].sum(1)
    return vB, idx, off, diag


def apply_rows(idx, off, diag, psi):
    return diag[:, None] * psi + np.einsum("nj,njc->nc", off, psi[idx])


def banded(idx, off, diag, m, scale):
    """Banded matrix I − scale·L (bandwidth 2m)."""
    N = idx.shape[0]; b = 2 * m
    ab = np.zeros((2 * b + 1, N), complex)
    rows = np.repeat(np.arange(N), idx.shape[1])
    cols = idx.ravel()
    ab[b + rows - cols, cols] = -scale * off.ravel()
    ab[b] += 1.0 - scale * diag                                 # self column of off is 0
    return ab


def herm_step(x, psi, h):
    """Eulerian half of the rule 'H': the balls stand still and exchange wave with their two nearest
    neighbours through symmetric links (a discrete Schrödinger equation on the balls):
        h_i dψ_i/dz = (i/2k) [ c_+ (ψ_{i+1} − ψ_i)/g_+  −  c_− (ψ_i − ψ_{i−1})/g_− ],   h_i = (g_+ + g_−)/2
    c = (θ/2)²/sin²(θ/2) (θ = phase difference across the link) makes it exact for a local plane wave.
    The links are symmetric, so Σ h_i |ψ_i|² is conserved exactly (Crank-Nicolson): it cannot blow up,
    whatever the arrangement of the balls."""
    N = x.size
    g = np.diff(x)
    th = np.abs(np.angle(np.einsum("nc,nc->n", psi[1:], np.conj(psi[:-1]))))
    th = np.minimum(th, 2.5)
    c = np.where(th > 1e-6, (th / 2) ** 2 / np.maximum(np.sin(th / 2) ** 2, 1e-300), 1.0)
    kap = c / g                                                  # link strength
    hh = np.empty(N); hh[1:-1] = 0.5 * (g[1:] + g[:-1]); hh[0] = 0.5 * g[0]; hh[-1] = 0.5 * g[-1]
    lo = np.concatenate([[0.0], kap]); up = np.concatenate([kap, [0.0]])
    diag = -(lo + up) / hh
    sub = kap / hh[1:]                                           # A[i, i-1]
    sup = kap / hh[:-1]                                          # A[i, i+1]
    a = 1j / (2 * K)
    def Aop(q):
        r = diag[:, None] * q
        r[1:] += sub[:, None] * q[:-1]; r[:-1] += sup[:, None] * q[1:]
        return a * r
    rhs = psi + 0.5 * h * Aop(psi)
    ab = np.zeros((3, N), complex)
    ab[0, 1:] = -0.5 * h * a * sup
    ab[1] = 1 - 0.5 * h * a * diag
    ab[2, :-1] = -0.5 * h * a * sub
    return solve_banded((1, 1), ab, rhs)


def run_psi(x, psi, zs, m=5, P=4, dz_fn=None, z0=0.0, scheme="cn", tmax=500, v0=None, log=None, grow=None, cfl_x=None,
            dzmin=1e-6):
    """Rule C. scheme 'cn': trapezoidal (Crank-Nicolson) predictor-corrector, 'rk4': explicit RK4."""
    x = x.copy(); psi = np.array(psi, complex).reshape(x.size, -1).copy()
    v = np.zeros(x.size) if v0 is None else v0.copy()
    z = z0; out = []; outpsi = []; steps = swaps = 0; t = time.time(); b = 2 * m
    lab = np.arange(x.size)                                     # who is who (to follow carried weights)
    for _ in range(2):                                          # settle the velocity state on the initial ψ
        v, *_ = psi_fit(x, psi, v, m, P)
    for zt in zs:
        while z < zt - 1e-12:
          try:
            h = min(dz_fn(x, z), zt - z)
            if cfl_x:                                           # neighbours must not overtake each other in one step
                closing = np.maximum(-(np.diff(v)), 1e-12)
                h = max(min(h, cfl_x * np.min(np.diff(x) / closing)), min(dzmin, zt - z))
            if scheme == "herm":
                # 1) wave exchange between standing balls (unitary), 2) every ball moves with the local
                #    current and reads the wave at its new place from its neighbours (local quadratic fit)
                psi = herm_step(x, psi, h)
                v, _ = psi_fit_log(x, psi, m, P, DEMOD["clip"], parts=True)
                a1, a2, c1, c2 = psi_fit_log.last
                dl = h * v
                E = np.exp(a1 * dl + a2 * dl * dl)
                psi = E[:, None] * (psi + c1 * dl[:, None] + c2 * (dl * dl)[:, None])
                x = x + dl
            elif scheme == "cn2":
                # consistent trapezoid: the SAME ball speed ū is used for the motion and in both operators
                v0_, pr0 = psi_fit_log(x, psi, m, P, DEMOD["clip"], parts=True)
                i0, o0, d0 = rows_from_parts(pr0, v0_)
                ps_ = solve_banded((b, b), banded(i0, o0, d0, m, 0.5 * h), psi + 0.5 * h * apply_rows(i0, o0, d0, psi))
                xs_ = x + h * v0_
                o = np.argsort(xs_, kind="stable")
                vs_ = np.empty_like(v0_); vs_[o] = psi_fit_log(xs_[o], ps_[o], m, P, DEMOD["clip"], parts=True)[0]
                vbar = 0.5 * (v0_ + vs_)
                i0, o0, d0 = rows_from_parts(pr0, vbar)
                th_ = DEMOD.get("theta", 0.5)                   # 0.5 = Crank-Nicolson, > 0.5 damps ball-scale noise
                rhs = psi + (1 - th_) * h * apply_rows(i0, o0, d0, psi)
                x1 = x + h * vbar
                o = np.argsort(x1, kind="stable")
                if (o != np.arange(o.size)).any():
                    swaps += int((o != np.arange(o.size)).sum())
                x1, rhs, vbar, lab = x1[o], rhs[o], vbar[o], lab[o]
                _, i1, o1, d1 = psi_fit_log(x1, psi[o], m, P, DEMOD["clip"], v_move=vbar)   # carrier from old ψ (frozen)
                psi = solve_banded((b, b), banded(i1, o1, d1, m, th_ * h), rhs)
                x = x1; v = vbar
            elif scheme == "rk4":
                def f(xx, pp, vv):
                    vB, idx, off, diag = psi_fit(xx, pp, vv, m, P)
                    return vB, apply_rows(idx, off, diag, pp)
                v1, p1 = f(x, psi, v)
                v2, p2 = f(x + .5 * h * v1, psi + .5 * h * p1, v1)
                v3, p3 = f(x + .5 * h * v2, psi + .5 * h * p2, v2)
                v4, p4 = f(x + h * v3, psi + h * p3, v3)
                x = x + h / 6 * (v1 + 2 * v2 + 2 * v3 + v4)
                psi = psi + h / 6 * (p1 + 2 * p2 + 2 * p3 + p4)
                v = v4
            else:
                v0_, idx, off, diag = psi_fit(x, psi, v, m, P)
                rhs = psi + 0.5 * h * apply_rows(idx, off, diag, psi)
                xs_ = x + h * v0_
                ps_ = solve_banded((b, b), banded(idx, off, diag, m, 0.5 * h), rhs)
                if np.any(np.diff(xs_) <= 0):                   # predictor crossed: keep ranks, re-sort
                    o = np.argsort(xs_, kind="stable")
                else:
                    o = None
                # corrector: operator at the predicted end point (in the predictor's rank order)
                xp, pp, vp = (xs_, ps_, v0_) if o is None else (xs_[o], ps_[o], v0_[o])
                v1_, idx1, off1, diag1 = psi_fit(xp, pp, vp, m, P)
                if o is not None:
                    inv = np.empty_like(o); inv[o] = np.arange(o.size)
                    # express the corrector system in the original order
                    Lp = lambda q: apply_rows(idx1, off1, diag1, q[o])[inv]
                    # fall back to fixed-point iterations (rare)
                    p_new = ps_.copy()
                    for _ in range(20):
                        p_new = rhs + 0.5 * h * Lp(p_new)
                    v1o = v1_[inv]
                else:
                    p_new = solve_banded((b, b), banded(idx1, off1, diag1, m, 0.5 * h), rhs)
                    v1o = v1_
                x = x + 0.5 * h * (v0_ + v1o)
                psi = p_new; v = v1o
          except (ValueError, np.linalg.LinAlgError):              # NaN/inf reached the linear algebra
            return (out + [None] * (len(zs) - len(out)), outpsi,
                    dict(steps=steps, swaps=swaps, time=time.time() - t, z=z, fail=True))
          z += h; steps += 1
          if True:
            if np.any(np.diff(x) <= 0):
                o = np.argsort(x, kind="stable"); swaps += int((o != np.arange(x.size)).sum())
                x, psi, v, lab = x[o], psi[o], v[o], lab[o]
            if not (np.all(np.isfinite(x)) and np.all(np.isfinite(psi))) or time.time() - t > tmax:
                return (out + [None] * (len(zs) - len(out)), outpsi,
                        dict(steps=steps, swaps=swaps, time=time.time() - t, z=z, fail=True))
            if log and steps % log == 0:
                print(f"    z={z:.3f} steps={steps} dz={h:.2e} width={2*x.std():.3f}", flush=True)
        out.append(x.copy()); outpsi.append(psi.copy())
    return out, outpsi, dict(steps=steps, swaps=swaps, time=time.time() - t, z=z, fail=False, lab=lab)


# ───────────────────────────────── initial conditions ───────────────────────────────
def free_start(N, seed=None):
    """Balls emitted at the waist with density |ψ|² (regular quantiles if seed is None, else random)."""
    sig = W0F / 2
    if seed is None:
        x = sig * ndtri((np.arange(N) + 0.5) / N)
    else:
        x = np.sort(np.random.default_rng(seed).normal(0, sig, N))
    return x, np.exp(-x ** 2 / W0F ** 2).astype(complex)


def slit_start(N, z0=1.0, seed=None, xr=12.0, ng=24001, m_src=400):
    """Flux-based start behind the slits (previous study): positions distributed as the exact RS flux
    j_z(x, z0); each ball also receives the local exact field value U(x_i, z0) (what the slit emits)."""
    s, w = sources(m_src)
    xg = np.linspace(-xr, xr, ng)
    U, Ux, Uz = rs_field(xg, z0, s, w)
    jz = np.maximum(flux(U, Ux, Uz)[1], 0)
    cdf = np.concatenate([[0], np.cumsum(0.5 * (jz[1:] + jz[:-1]) * np.diff(xg))]); cdf /= cdf[-1]
    q = (np.arange(N) + 0.5) / N if seed is None else np.sort(np.random.default_rng(seed).random(N))
    x = np.interp(q, cdf, xg)
    Ui, Uxi, _ = rs_field(x, z0, s, w)
    return x, Ui[:, 0], Uxi[:, 0]


_REF = {}


def slit_refs(z0=1.0):
    """exact RS flux on the screen (the target) and the paraxial evolution of the same start field
    (FFT, 'what a perfect paraxial rule gives'), both per fine screen bin, normalized."""
    if z0 in _REF:
        return _REF[z0]
    ref = screen_reference("flux")
    s, w = sources(400)
    n = 2 ** 16; X = 800.0
    xg = (np.arange(n) - n // 2) * (X / n)
    near = np.abs(xg) < 40
    U0 = np.zeros(n, complex)
    U0[near] = rs_field(xg[near], z0, s, w)[0][:, 0]
    kx = 2 * np.pi * np.fft.fftfreq(n, X / n)
    UL = np.fft.ifft(np.fft.fft(U0) * np.exp(-1j * kx ** 2 * (L - z0) / (2 * K)))
    I = np.abs(UL) ** 2
    par = np.array([I[(xg >= a) & (xg < b)].sum() for a, b in zip(EDGES[:-1], EDGES[1:])])
    par /= par.sum()
    _REF[z0] = (ref, par)
    return ref, par


def quantile_sampler(ref, N):
    """Balls placed at the regular quantiles of a reference pattern (the ideal deterministic photon)."""
    cdf = np.concatenate([[0], np.cumsum(ref)]); cdf /= cdf[-1]
    return np.interp((np.arange(N) + 0.5) / N, cdf, EDGES)


def score(xl, ref):
    c = coarse_corr(xl, ref)
    _, rms, mx = fringe_errors(xl, ref)
    return c, rms, mx
