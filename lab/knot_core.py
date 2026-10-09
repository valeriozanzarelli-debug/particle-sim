"""knot_core.py - Skyrme model on a 3D lattice of 'balls' (each ball = unit 4-vector phi = (sigma, pi1, pi2, pi3)).

Skyrme units (energy f_pi/(4e) with ANW F_pi convention, length 2/(e F_pi)):
  static energy density  e = |d phi|^2 + 1/2 [ (tr D)^2 - tr(D^2) ] + 2 m^2 (1 - sigma),   D_ij = d_i phi . d_j phi
  kinetic energy density T = |v|^2 (1 + |d phi|^2) - sum_i (v . d_i phi)^2                (v = d phi / dt)
  -> Lorentz invariant, waves travel at c = 1, Faddeev-Bogomolny bound E >= 12 pi^2 |B|.
Derivatives: 4th-order central differences (Battye-Sutcliffe style), Dirichlet boundary phi = (1,0,0,0).
Dynamics: Hamiltonian form, pi = 2 M v, H = sum h^3 [ pi M^-1 pi / 4 + e ], Dirac-bracket projection onto the
sphere, RK4 + renormalisation (constraint projection) every step. M^-1 by Woodbury (3x3 adjugate, all local).
"""
import math, numpy as np, torch
from scipy.integrate import solve_ivp

torch.set_num_threads(3)
TWO_PI2 = 2 * math.pi ** 2


# ---------------------------------------------------------------------------------------- radial hedgehog (exact ref)
def _rhs(x, y):        # ANW units x = e F_pi r = 2 r_skyrme
    F, Fp = y
    s, s2F = np.sin(F), np.sin(2 * F)
    return [Fp, -(x * Fp / 2 + s2F * Fp**2 - s2F / 4 - s * s * s2F / x**2) / (x**2 / 4 + 2 * s * s)]


def radial_profile():
    """Shooting for the massless B=1 hedgehog. Returns r (Skyrme units), f(r), and dict of radial integrals."""
    def shoot(a, xmax=30.0):
        x0 = 1e-3
        ez = lambda x, y: y[0]; ez.terminal = True
        eu = lambda x, y: y[1]; eu.terminal = True
        return solve_ivp(_rhs, [x0, xmax], [np.pi - a * x0, -a], rtol=1e-12, atol=1e-14,
                         events=[ez, eu], dense_output=True)
    lo, hi = 1.0, 3.0
    for _ in range(60):
        am = 0.5 * (lo + hi)
        if shoot(am).t_events[0].size: hi = am
        else: lo = am
    a = 0.5 * (lo + hi)
    s_ = shoot(a)
    xx = np.linspace(1e-3, s_.t[-1], 200001)
    FF = s_.sol(xx)[0]
    im = int(np.argmax(FF < 0.05))
    A = FF[im] * xx[im] ** 2
    xt = np.geomspace(xx[im], 5e3, 20001)[1:]
    x = np.concatenate([xx[:im + 1], xt]); F = np.concatenate([FF[:im + 1], A / xt**2])
    r = x / 2.0
    f = F
    fp = np.gradient(f, r)
    s2 = np.sin(f) ** 2
    I = lambda g: np.trapezoid(g, r)
    E2 = 4 * np.pi * I(r**2 * fp**2 + 2 * s2)
    E4 = 4 * np.pi * I(2 * s2 * fp**2 + s2**2 / r**2)
    B = -(2 / np.pi) * I(s2 * fp)
    Lam = (16 * np.pi / 3) * I(r**2 * s2 * (1 + fp**2 + s2 / r**2))     # T = Lam w^2 / 2
    eden = (fp**2 + 2 * s2 / r**2) + (2 * s2 * fp**2 / r**2 + s2**2 / r**4)
    bden = -(1 / (2 * np.pi**2)) * s2 * fp / r**2
    r_half = float(np.interp(np.pi / 2, f[::-1], r[::-1]))
    rE2 = 4 * np.pi * I(r**4 * eden) / (E2 + E4)
    rB2 = 4 * np.pi * I(r**4 * bden)
    return dict(r=r, f=f, fp=fp, a=a * 2, A_tail=A / 4, E=E2 + E4, E2=E2, E4=E4, B=B, Lam=Lam,
                r_half=r_half, rE2=rE2, rB2=rB2)


def ball_energy(R, n=4000):
    """Hedgehog minimised inside a ball of radius R with f(R)=0 (Dirichlet), 1D grid. Returns E."""
    r = torch.linspace(0, R, n + 1, dtype=torch.float64)
    h = R / n
    f0 = np.pi * np.exp(-(r.numpy() / 0.9) ** 1.5) * (1 - (r.numpy() / R) ** 4)
    fi = torch.tensor(f0[1:-1], requires_grad=True)
    opt = torch.optim.LBFGS([fi], lr=1, max_iter=4000, tolerance_grad=1e-12, tolerance_change=1e-15,
                            history_size=50, line_search_fn="strong_wolfe")
    rm = 0.5 * (r[1:] + r[:-1])

    def E():
        f = torch.cat([torch.tensor([np.pi], dtype=torch.float64), fi, torch.tensor([0.0], dtype=torch.float64)])
        fp = (f[1:] - f[:-1]) / h
        s2 = torch.sin(0.5 * (f[1:] + f[:-1])) ** 2
        return 4 * np.pi * h * torch.sum(rm**2 * fp**2 + 2 * s2 + 2 * s2 * fp**2 + s2**2 / rm**2)

    def clo():
        opt.zero_grad(); e = E(); e.backward(); return e
    for _ in range(3):
        opt.step(clo)
    with torch.no_grad():
        return float(E())


# ---------------------------------------------------------------------------------------- lattice
VAC = torch.tensor([1.0, 0, 0, 0])


def grid(shape, h, dtype=torch.float64):
    ax = [(torch.arange(n, dtype=dtype) - (n - 1) / 2) * h for n in shape]
    return torch.meshgrid(*ax, indexing="ij")


def hedgehog(X, Y, Z, prof, center=(0, 0, 0), gx=1.0, sign=1, iso=None, taper=None):
    """phi = (cos f, sign sin f xhat) with x contracted by gx (boost) about center. iso: 3x3 isorotation."""
    x = (X - center[0]) * gx; y = Y - center[1]; z = Z - center[2]
    r = torch.sqrt(x * x + y * y + z * z)
    rn = prof["r"]; fn = prof["f"]
    f = torch.tensor(np.interp(r.cpu().numpy().ravel(), rn, fn, left=np.pi, right=0.0).reshape(r.shape), dtype=X.dtype)
    if taper is not None:      # smooth cut of the 1/r^2 tail so that phi -> vacuum before the Dirichlet wall
        R1, R2 = taper
        w = 0.5 * (1 + torch.cos(math.pi * ((r - R1) / (R2 - R1)).clamp(0, 1)))
        f = f * w
    rs = torch.where(r > 0, r, torch.ones_like(r))
    s = torch.sin(f) / rs * sign
    pi = torch.stack([s * x, s * y, s * z])
    if iso is not None:
        pi = torch.einsum("ab,b...->a...", torch.tensor(iso, dtype=X.dtype), pi)
    return torch.cat([torch.cos(f)[None], pi])


def qmul(a, b):
    """quaternion product of fields (4,...) : U1 U2 for U = sigma + i pi.tau"""
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return torch.stack([a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
                        a0 * b1 + a1 * b0 - (a2 * b3 - a3 * b2),
                        a0 * b2 + a2 * b0 - (a3 * b1 - a1 * b3),
                        a0 * b3 + a3 * b0 - (a1 * b2 - a2 * b1)])


def pad(phi):
    vac = VAC.to(phi.dtype).view(4, 1, 1, 1)
    return torch.nn.functional.pad((phi - vac)[None], (2, 2, 2, 2, 2, 2))[0] + vac


def derivs(phi, h):
    P = pad(phi)
    n = phi.shape[1:]
    out = []
    for ax in range(3):
        sl = [slice(None), slice(2, 2 + n[0]), slice(2, 2 + n[1]), slice(2, 2 + n[2])]
        def sh(o):
            s = list(sl); s[ax + 1] = slice(2 + o, 2 + o + n[ax]); return P[tuple(s)]
        out.append((-sh(2) + 8 * sh(1) - 8 * sh(-1) + sh(-2)) / (12 * h))
    return out


def derivs_fb(phi, h):
    """forward and backward first differences (no zero 'checkerboard' mode)."""
    P = pad(phi)
    n = phi.shape[1:]
    gp, gm = [], []
    for ax in range(3):
        sl = [slice(None), slice(2, 2 + n[0]), slice(2, 2 + n[1]), slice(2, 2 + n[2])]
        def sh(o):
            s = list(sl); s[ax + 1] = slice(2 + o, 2 + o + n[ax]); return P[tuple(s)]
        c = sh(0)
        gp.append((sh(1) - c) / h); gm.append((c - sh(-1)) / h)
    return gp, gm


def static_density(phi, h, m=0.0, g=None):
    """sigma term: mean of forward/backward |d phi|^2; Skyrme term: sum_{i<j} |g_i ^ g_j|^2 averaged over the
    4 forward/backward combinations (2nd order, symmetric, positive, no zero modes)."""
    gp, gm = derivs_fb(phi, h) if g is None else g
    n2p = [(x * x).sum(0) for x in gp]; n2m = [(x * x).sum(0) for x in gm]
    e = 0.5 * sum(n2p[i] + n2m[i] for i in range(3))
    for i in range(3):
        for j in range(i + 1, 3):
            for gi, ni in ((gp[i], n2p[i]), (gm[i], n2m[i])):
                for gj, nj in ((gp[j], n2p[j]), (gm[j], n2m[j])):
                    e = e + 0.25 * (ni * nj - (gi * gj).sum(0) ** 2)
    if m: e = e + 2 * m * m * (1 - phi[0])
    return e


def baryon_density(phi, h):
    """B0 = -(1/2pi^2) det[phi, d1 phi, d2 phi, d3 phi] (4th-order central derivatives; measurement only)."""
    g1, g2, g3 = derivs(phi, h)
    def m3(a, b, c, i, j, k):
        return (a[i] * (b[j] * c[k] - b[k] * c[j]) - a[j] * (b[i] * c[k] - b[k] * c[i]) + a[k] * (b[i] * c[j] - b[j] * c[i]))
    det = (phi[0] * m3(g1, g2, g3, 1, 2, 3) - phi[1] * m3(g1, g2, g3, 0, 2, 3)
           + phi[2] * m3(g1, g2, g3, 0, 1, 3) - phi[3] * m3(g1, g2, g3, 0, 1, 2))
    return -det / TWO_PI2


def energy(phi, h, m=0.0):
    return static_density(phi, h, m).sum() * h**3


def normalize(phi):
    return phi / phi.norm(dim=0, keepdim=True)


def tangent(phi, w):
    return w - phi * (phi * w).sum(0, keepdim=True)


def relax(phi, h, m=0.0, dt=None, maxit=4000, tol=1e-10, log=None, every=200):
    """Arrested Newton flow. Returns phi, history of E."""
    dt = dt or 0.25 * h
    v = torch.zeros_like(phi)
    Eold = None; hist = []
    for it in range(maxit):
        p = phi.detach().requires_grad_(True)
        E = energy(p, h, m)
        (gr,) = torch.autograd.grad(E, p)
        E = float(E.detach())
        if Eold is not None and E > Eold:
            v.zero_()
        F = -tangent(phi, gr) / h**3
        v = tangent(phi, v + dt * F)
        phi = normalize(phi + dt * v)
        hist.append(E)
        if log and it % every == 0:
            log(f"    it {it:5d}  E/12pi^2 = {E/(12*math.pi**2):.6f}  |F|max = {float(F.abs().max()):.2e}")
        if it > 300 and abs(hist[-100] - E) / E < tol * 100:
            break
        Eold = E
    return phi.detach(), hist


# ---------------------------------------------------------------------------------------- dynamics
def _minv(phi, g, pi):
    """M^-1 pi on the tangent space (Woodbury, all local). g = (gp, gm)."""
    gp, gm = g
    a = 1 + 0.5 * sum((gp[i] * gp[i]).sum(0) + (gm[i] * gm[i]).sum(0) for i in range(3))
    g = [0.5 * (gp[i] + gm[i]) for i in range(3)]
    Gt = [gi - phi * (phi * gi).sum(0, keepdim=True) for gi in g]
    D = [[(Gt[i] * Gt[j]).sum(0) for j in range(3)] for i in range(3)]
    S = [[(a if i == j else 0) - D[i][j] for j in range(3)] for i in range(3)]
    c00 = S[1][1] * S[2][2] - S[1][2] * S[2][1]
    c01 = S[0][2] * S[2][1] - S[0][1] * S[2][2]
    c02 = S[0][1] * S[1][2] - S[0][2] * S[1][1]
    c11 = S[0][0] * S[2][2] - S[0][2] * S[2][0]
    c12 = S[0][2] * S[1][0] - S[0][0] * S[1][2]
    c22 = S[0][0] * S[1][1] - S[0][1] * S[1][0]
    det = S[0][0] * c00 + S[0][1] * (S[1][2] * S[2][0] - S[1][0] * S[2][2]) + S[0][2] * (S[1][0] * S[2][1] - S[1][1] * S[2][0])
    Inv = [[c00, c01, c02], [c01, c11, c12], [c02, c12, c22]]
    w = [(Gt[i] * pi).sum(0) for i in range(3)]
    y = [(Inv[i][0] * w[0] + Inv[i][1] * w[1] + Inv[i][2] * w[2]) / det for i in range(3)]
    return (pi + Gt[0] * y[0] + Gt[1] * y[1] + Gt[2] * y[2]) / a


def momentum_from_velocity(phi, v, h):
    """pi = 2 M v (tangent)."""
    gp, gm = derivs_fb(phi, h)
    a = 1 + 0.5 * sum((gp[i] * gp[i]).sum(0) + (gm[i] * gm[i]).sum(0) for i in range(3))
    g = [0.5 * (gp[i] + gm[i]) for i in range(3)]
    Gt = [gi - phi * (phi * gi).sum(0, keepdim=True) for gi in g]
    v = tangent(phi, v)
    return 2 * tangent(phi, a * v - sum(Gt[i] * (Gt[i] * v).sum(0, keepdim=True) for i in range(3)))


def ham_terms(phi, pi, h, m=0.0, need_grad=True):
    """returns H (float), kinetic density, static density, phidot, dH/dphi (per volume)"""
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


def sponge_profile(shape, h, width, strength, dtype):
    X, Y, Z = grid(shape, h, dtype)
    s = torch.zeros_like(X)
    for A, n in zip((X, Y, Z), shape):
        L = (n - 1) / 2 * h
        d = (A.abs() - (L - width)).clamp(min=0) / width
        s = torch.maximum(s, strength * d**2)
    return s[None]


# ---------------------------------------------------------------------------------------- marching tetrahedra
_TETS = [(0, 1, 3, 7), (0, 3, 2, 7), (0, 2, 6, 7), (0, 6, 4, 7), (0, 4, 5, 7), (0, 5, 1, 7)]
_CORN = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]])


def isosurface(F, level, h, origin):
    """Marching tetrahedra on a 3D numpy array. Returns triangles (T,3,3) in physical coords."""
    nx, ny, nz = F.shape
    I, J, K = np.meshgrid(np.arange(nx - 1), np.arange(ny - 1), np.arange(nz - 1), indexing="ij")
    base = np.stack([I.ravel(), J.ravel(), K.ravel()], 1)
    vals = np.stack([F[I + c[0], J + c[1], K + c[2]].ravel() for c in _CORN], 1) - level
    keep = (vals.min(1) < 0) & (vals.max(1) > 0)
    base, vals = base[keep], vals[keep]
    tris = []
    for t in _TETS:
        v = vals[:, t]; P = (base[:, None, :] + _CORN[list(t)][None]).astype(float)
        inside = v > 0
        cnt = inside.sum(1)
        for c1 in (1, 3):     # one vertex separated
            sel = np.where(cnt == c1)[0]
            if not len(sel): continue
            odd = np.where(inside[sel] if c1 == 1 else ~inside[sel])[1]
            oth = np.array([[j for j in range(4) if j != o] for o in odd])
            tri = []
            for q in range(3):
                b = oth[:, q]
                va = v[sel, odd]; vb = v[sel, b]
                w = (va / (va - vb))[:, None]
                tri.append(P[sel, odd] + w * (P[sel, b] - P[sel, odd]))
            tris.append(np.stack(tri, 1))
        sel = np.where(cnt == 2)[0]
        if len(sel):
            ins = np.argsort(~inside[sel], axis=1, kind="stable")    # first two: inside
            a0, a1, b0, b1 = ins[:, 0], ins[:, 1], ins[:, 2], ins[:, 3]

            def cp(a, b):
                va = v[sel, a]; vb = v[sel, b]
                w = (va / (va - vb))[:, None]
                return P[sel, a] + w * (P[sel, b] - P[sel, a])
            p1, p2, p3, p4 = cp(a0, b0), cp(a0, b1), cp(a1, b1), cp(a1, b0)
            tris.append(np.stack([p1, p2, p3], 1)); tris.append(np.stack([p1, p3, p4], 1))
    T = np.concatenate(tris, 0) if tris else np.zeros((0, 3, 3))
    return T * h + np.asarray(origin)[None, None, :]
