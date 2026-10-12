"""knot2_core.py - Skyrme lattice engine with a SELECTABLE stencil, for the knot2 study (S1-S4).
Reuses knot_core.py (read-only) for the radial ODE, hedgehog, qmul, baryon density, sponge, isosurface.

Schemes (set with set_scheme):
  "old" : forward/backward 2nd-order differences in the static energy; kinetic matrix
          M = (1 + 1/2 sum_s |g^s|^2) I - sum_i gbar_i gbar_i^T  with gbar = (g+ + g-)/2  (knot_core.py, K1-K5).
          Static and kinetic Skyrme pieces are built with DIFFERENT stencils.
  "c2"  : same static energy; kinetic matrix built like the static one (time = one more direction):
          M = (1 + 1/2 sum_{i,s} |g_i^s|^2) I - 1/2 sum_{i,s} Gt_i^s Gt_i^s^T     (s = forward, backward)
  "c4"  : like c2 but with the 4th-order STAGGERED forward/backward differences
          g+ = (-phi[x+2] + 27 phi[x+1] - 27 phi[x] + phi[x-1]) / (24 h),  g-(x) = g+(x-1)
          (no zero 'checkerboard' mode: 27 sin(q) - sin(3q) = 0 only at q = 0), same stencil in kinetic and static.
PERIODIC (module flag) replaces the Dirichlet vacuum padding by periodic padding (used only by the S1 wave test).
"""
import math, numpy as np, torch
import knot_core as kc

torch.set_num_threads(2)
SCHEME = "c4"
PERIODIC = False
TWO_PI2 = kc.TWO_PI2
grid, hedgehog, qmul, baryon_density, normalize, tangent = kc.grid, kc.hedgehog, kc.qmul, kc.baryon_density, kc.normalize, kc.tangent
sponge_profile, isosurface, radial_bvp, derivs = kc.sponge_profile, kc.isosurface, kc.radial_bvp, kc.derivs


def set_scheme(s):
    global SCHEME
    assert s in ("old", "c2", "c4")
    SCHEME = s


def pad(phi):
    if PERIODIC:
        return torch.nn.functional.pad(phi[None], (2, 2, 2, 2, 2, 2), mode="circular")[0]
    return kc.pad(phi)


def derivs_fb(phi, h):
    P = pad(phi)
    n = phi.shape[1:]
    gp, gm = [], []
    for ax in range(3):
        sl = [slice(None), slice(2, 2 + n[0]), slice(2, 2 + n[1]), slice(2, 2 + n[2])]

        def sh(o):
            s = list(sl); s[ax + 1] = slice(2 + o, 2 + o + n[ax]); return P[tuple(s)]
        c = sh(0)
        if SCHEME == "c4":
            p1, m1 = sh(1), sh(-1)
            gp.append((-sh(2) + 27 * p1 - 27 * c + m1) / (24 * h))
            gm.append((-p1 + 27 * c - 27 * m1 + sh(-2)) / (24 * h))
        else:
            gp.append((sh(1) - c) / h); gm.append((c - sh(-1)) / h)
    return gp, gm


def static_density(phi, h, m=0.0, g=None):
    return kc.static_density(phi, h, m, derivs_fb(phi, h) if g is None else g)


def energy(phi, h, m=0.0):
    return static_density(phi, h, m).sum() * h**3


def _U(phi, g):
    gp, gm = g
    a = 1 + 0.5 * sum((gp[i] * gp[i]).sum(0) + (gm[i] * gm[i]).sum(0) for i in range(3))
    if SCHEME == "old":
        G = [0.5 * (gp[i] + gm[i]) for i in range(3)]
    else:
        s = math.sqrt(0.5)
        G = [s * x for i in range(3) for x in (gp[i], gm[i])]
    U = [gi - phi * (phi * gi).sum(0, keepdim=True) for gi in G]
    return a, U


def _minv(phi, g, pi):
    """M^-1 pi, M = a I - sum_k U_k U_k^T on the tangent space (Woodbury, local k x k solve)."""
    a, U = _U(phi, g)
    if SCHEME == "old":
        return kc._minv(phi, g, pi)
    # B = sum_k U_k U_k^T (4x4, B phi = 0); A = a - B on the 3-dim tangent space; Cayley-Hamilton:
    # A^-1 = (A^2 - e1 A + e2) / e3 with the invariants of A written through the principal minors of B.
    idx = [(i, j) for i in range(4) for j in range(i, 4)]
    Bm = {}
    for i, j in idx:
        Bm[i, j] = sum(u[i] * u[j] for u in U); Bm[j, i] = Bm[i, j]
    mv = lambda x: torch.stack([sum(Bm[i, j] * x[j] for j in range(4)) for i in range(4)])
    e1 = Bm[0, 0] + Bm[1, 1] + Bm[2, 2] + Bm[3, 3]
    e2 = sum(Bm[i, i] * Bm[j, j] - Bm[i, j] ** 2 for i in range(4) for j in range(i + 1, 4))

    def m3(i, j, k):
        return (Bm[i, i] * (Bm[j, j] * Bm[k, k] - Bm[j, k] ** 2) - Bm[i, j] * (Bm[i, j] * Bm[k, k] - Bm[j, k] * Bm[i, k])
                + Bm[i, k] * (Bm[i, j] * Bm[j, k] - Bm[j, j] * Bm[i, k]))
    e3 = m3(1, 2, 3) + m3(0, 2, 3) + m3(0, 1, 3) + m3(0, 1, 2)
    E1 = 3 * a - e1; E2 = 3 * a * a - 2 * a * e1 + e2; E3 = a**3 - a * a * e1 + a * e2 - e3
    Bp = mv(pi); BBp = mv(Bp)
    Ap = a * pi - Bp
    A2p = a * a * pi - 2 * a * Bp + BBp
    return (A2p - E1 * Ap + E2 * pi) / E3


def mass_apply(phi, v, h):
    """M v (tangent)."""
    g = derivs_fb(phi, h)
    a, U = _U(phi, g)
    v = tangent(phi, v)
    return tangent(phi, a * v - sum(u * (u * v).sum(0, keepdim=True) for u in U))


def momentum_from_velocity(phi, v, h):
    return 2 * mass_apply(phi, v, h)


def ham_terms(phi, pi, h, m=0.0, need_grad=True):
    p = phi.detach().requires_grad_(need_grad)
    g = derivs_fb(p, h)
    mi = _minv(p, g, pi)
    kin = 0.25 * (pi * mi).sum(0)
    st = static_density(p, h, m, g)
    H = (kin + st).sum() * h**3
    gr = torch.autograd.grad(H, p)[0] / h**3 if need_grad else None
    return float(H), kin.detach(), st.detach(), tangent(phi, 0.5 * mi.detach()), gr


def rhs(phi, pi, h, m=0.0, sponge=None):
    _, _, _, phid, gr = ham_terms(phi, pi, h, m)
    pid = -tangent(phi, gr) - phi * (pi * phid).sum(0, keepdim=True)
    if sponge is not None:
        pid = pid - sponge * pi
    return phid, pid


def rk4(phi, pi, h, dt, m=0.0, sponge=None):
    k1 = rhs(phi, pi, h, m, sponge)
    k2 = rhs(phi + 0.5 * dt * k1[0], pi + 0.5 * dt * k1[1], h, m, sponge)
    k3 = rhs(phi + 0.5 * dt * k2[0], pi + 0.5 * dt * k2[1], h, m, sponge)
    k4 = rhs(phi + dt * k3[0], pi + dt * k3[1], h, m, sponge)
    phi = normalize(phi + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]))
    pi = tangent(phi, pi + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]))
    return phi, pi


def relax_lbfgs(phi, h, m=0.0, maxit=3000, tol=1e-11, log=None, every=5, wall=None, hist=None):
    import time as _t
    t0 = _t.time()
    psi = phi.detach().clone().requires_grad_(True)
    opt = torch.optim.LBFGS([psi], lr=1, max_iter=20, history_size=30, line_search_fn="strong_wolfe",
                            tolerance_grad=1e-14, tolerance_change=1e-16)
    hist = list(hist or [])

    def clo():
        opt.zero_grad(); E = energy(normalize(psi), h, m); E.backward(); return E
    done = False
    while len(hist) < maxit:
        opt.step(clo)
        with torch.no_grad():
            psi.data = normalize(psi.data)
        hist.append(float(energy(normalize(psi.detach()), h, m)))
        if log and len(hist) % every == 0:
            log(f"    lbfgs {len(hist):5d} (x20 inner)  E/12pi^2 = {hist[-1]/(12*math.pi**2):.8f}")
        if len(hist) > 10 and abs(hist[-6] - hist[-1]) / hist[-1] < 5 * tol:
            done = True; break
        if wall and _t.time() - t0 > wall:
            break
    return normalize(psi.detach()), hist, done


# ------------------------------------------------------------------ analytic lattice dispersion (vacuum, massless)
def omega_1d(kh, scheme):
    """omega * h for a lattice wave along one axis (kh in (0, pi])."""
    q = np.asarray(kh) / 2
    if scheme == "c4":
        return (27 * np.sin(q) - np.sin(3 * q)) / 12
    return 2 * np.sin(q)
