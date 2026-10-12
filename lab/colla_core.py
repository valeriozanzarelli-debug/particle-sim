"""colla_core.py - the 'glue' study: what is missing to get nuclear binding energies right from Skyrme knots.
Reuses knot_core.py / knot2_core.py / nuclei_core.py READ-ONLY (c4 stencil, hedgehog, rational maps, baryon density).

Generalised static energy (Skyrme units of knot_core; phi = (sigma, pi) unit 4-vector):
    e = ks |d phi|^2 + Skyrme + 2 m^2 ks (1 - sigma) + alpha (1 - sigma)^4 + c6 B0^2
    ks = 1 - alpha.  alpha = 0, c6 = 0 : standard (massive) Skyrme model.
    alpha > 0 : Gillard-Harland-Speight 2015 'lightly bound' model (Nucl. Phys. B 895, 272; eq. 3.1),
                there with m = 1 (their m = 2 m_pi sqrt(1-alpha)/(F_pi g)).
    c6 > 0  : sextic (omega-meson-like) term  c6 * B0^2  (B0 = baryon density, 4th-order central differences).
Kinetic metric (for inertia tensors and normal modes): T = sum h^3 v.M v,
    M v = a v - sum_k U_k (U_k . v),  a = ks + 1/2 sum_{i,s} |g_i^s|^2  (c4 forward/backward differences, as knot2_core 'c4').
Lagrangian L = T - E  ->  normal modes  K xi = w^2 (2M) xi,  K = Hessian of E on the sphere bundle.
Physical units: energy unit E_u = F_pi/(4e), length 2/(e F_pi), hbar = 2 e^2 (Skyrme units, hbar c = 197.327 MeV fm).
"""
import math, os, sys, time, numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot_core as kc, knot2_core as k2, nuclei_core as nc
torch.set_num_threads(2)
k2.set_scheme("c4")
PI2 = 12 * math.pi**2
HBARC = 197.327
M_N, M_D = 938.919, 1232.0       # nucleon (p,n average) and Delta masses (ANW / AN fit inputs)
M_PI = 138.0                      # pion mass used by Adkins-Nappi 1984
EXP_BIND = {2: 2.2246, 3: 7.718, 4: 28.296}   # deuteron, 3He, 4He binding energies (MeV)


# ------------------------------------------------------------------------------------------ 1D hedgehog
def hedgehog_1d(m=0.0, alpha=0.0, c6=0.0, R=14.0, n=7000):
    """B=1 hedgehog of the generalised model by direct minimisation (f(0)=pi, f(R)=0).
    Returns r, f, E, Lam (Lam without the sextic contribution)."""
    ks = 1 - alpha
    r = torch.linspace(0, R, n + 1, dtype=torch.float64); h = R / n
    rn = r.numpy()
    f0 = np.pi * np.exp(-(rn / 1.0)**1.5) * (1 - (rn / R)**4)
    fi = torch.tensor(f0[1:-1], requires_grad=True)
    rm = 0.5 * (r[1:] + r[:-1])
    pi_ = torch.tensor([np.pi], dtype=torch.float64); z_ = torch.tensor([0.0], dtype=torch.float64)

    def parts(fi):
        f = torch.cat([pi_, fi, z_]); fp = (f[1:] - f[:-1]) / h
        fm = 0.5 * (f[1:] + f[:-1]); s2 = torch.sin(fm)**2
        sig = rm**2 * fp**2 + 2 * s2
        sky = 2 * s2 * fp**2 + s2**2 / rm**2
        pot = rm**2 * (2 * m * m * ks * (1 - torch.cos(fm)) + alpha * (1 - torch.cos(fm))**4)
        sex = c6 * (s2 * fp)**2 / (4 * math.pi**4 * rm**2)
        return f, fp, fm, s2, 4 * math.pi * h * torch.sum(ks * sig + sky + pot + sex)
    opt = torch.optim.LBFGS([fi], lr=1, max_iter=400, tolerance_grad=1e-13, tolerance_change=1e-16,
                            history_size=60, line_search_fn="strong_wolfe")

    def clo():
        opt.zero_grad(); e = parts(fi)[-1]; e.backward(); return e
    opt.step(clo)
    # damped Newton with the exact tridiagonal Hessian (3 colour HVPs), banded solve
    from scipy.linalg import solve_banded
    x = fi.detach().clone()
    N = x.numel(); idx = torch.arange(N)
    for it in range(60):
        xr = x.clone().requires_grad_(True)
        E0 = parts(xr)[-1]
        (gr,) = torch.autograd.grad(E0, xr, create_graph=True)
        ab = np.zeros((3, N))
        for k in range(3):
            v = (idx % 3 == k).double()
            (hv,) = torch.autograd.grad(gr, xr, grad_outputs=v, retain_graph=True)
            hv = hv.detach().numpy()
            cols = np.arange(k, N, 3)
            ab[1, cols] = hv[cols]                                   # diagonal
            c1 = cols[cols + 1 < N]; ab[0, c1 + 1] = hv[c1 + 1]  # H[j, j+1] stored at ab[0, j+1]
            c2 = cols[cols - 1 >= 0]; ab[2, c2 - 1] = hv[c2 - 1]                   # H[j, j-1] stored at ab[2, j-1]
        gnp = gr.detach().numpy()
        try:
            dx = solve_banded((1, 1), ab, -gnp)
        except Exception:
            dx = -gnp
        step = 1.0; E0v = float(E0)
        while step > 1e-4:
            xn = x + step * torch.tensor(dx)
            if float(parts(xn)[-1]) <= E0v: break
            step *= 0.5
        x = xn
        if np.abs(dx).max() * step < 1e-11: break
    fi = x
    with torch.no_grad():
        f, fp, fm, s2, E = parts(fi)
        Lam = (16 * math.pi / 3) * h * torch.sum(rm**2 * s2 * (ks + fp**2 + s2 / rm**2))
        return dict(r=rn, f=f.numpy(), E=float(E), Lam=float(Lam))


def units_from_e_F(e, F):
    return dict(e=e, F_pi=F, E_unit=F / (4 * e), L_unit_fm=2 * HBARC / (e * F), hbar=2 * e * e)


def fit_units(m_pi=0.0, MN=M_N, MD=M_D):
    """Fit (e, F_pi) to M_N, M_Delta with M = E_u (E1 + J(J+1) hbar^2/(2 Lam)), J = 1/2, 3/2; m = 2 m_pi/(e F_pi)
    solved self-consistently when m_pi > 0. Massless: E1, Lam from knot_core.radial_bvp (exact tail)."""
    m = 0.0
    for it in range(40):
        if m_pi == 0:
            p = kc.radial_bvp(); E1, Lam = p["E"], p["Lam"]
        else:
            hh = hedgehog_1d(m=m); E1, Lam = hh["E"], hh["Lam"]
        # M_D - M_N = (F/4e)(3/2)(4e^4)/Lam = 1.5 F e^3/Lam ; M_N/(M_D-M_N) = (E1 + 1.5 e^4/Lam)/(6 e^4/Lam)
        rr = MN / (MD - MN)
        e = (E1 * Lam / (6 * rr - 1.5))**0.25
        F = (MD - MN) * Lam / (1.5 * e**3)
        m_new = 2 * m_pi / (e * F) if m_pi > 0 else 0.0
        conv = abs(m_new - m) < 1e-9
        m = m_new
        if m_pi == 0 or conv:
            break
    u = units_from_e_F(e, F); u.update(m=m, E1=E1, Lam=Lam, m_pi=m_pi)
    return u


# ------------------------------------------------------------------------------------------ lattice energy
def density(phi, h, m=0.0, alpha=0.0, c6=0.0, g=None):
    g = k2.derivs_fb(phi, h) if g is None else g
    gp, gm = g
    e = kc.static_density(phi, h, 0.0, g)           # sigma + Skyrme (c4 stencils)
    if alpha:
        sig = 0.5 * sum((gp[i] * gp[i]).sum(0) + (gm[i] * gm[i]).sum(0) for i in range(3))
        e = e - alpha * sig + alpha * (1 - phi[0])**4
    if m:
        e = e + 2 * m * m * (1 - alpha) * (1 - phi[0])
    if c6:
        e = e + c6 * kc.baryon_density(phi, h)**2
    return e


def energy(phi, h, **kw):
    return density(phi, h, **kw).sum() * h**3


def relax(phi, h, maxit=600, tol=1e-9, log=None, every=10, wall=None, hist=None, **kw):
    t0 = time.time()
    psi = phi.detach().clone().requires_grad_(True)
    opt = torch.optim.LBFGS([psi], lr=1, max_iter=20, history_size=30, line_search_fn="strong_wolfe",
                            tolerance_grad=1e-14, tolerance_change=1e-16)
    hist = list(hist or [])

    def clo():
        opt.zero_grad(); E = energy(kc.normalize(psi), h, **kw); E.backward(); return E
    done = False
    while len(hist) < maxit:
        opt.step(clo)
        with torch.no_grad():
            psi.data = kc.normalize(psi.data)
        hist.append(float(energy(kc.normalize(psi.detach()), h, **kw)))
        if log and len(hist) % every == 0:
            log(f"    lbfgs {len(hist):4d}  E = {hist[-1]:.8f}")
        if len(hist) > 10 and abs(hist[-6] - hist[-1]) / abs(hist[-1]) < 5 * tol:
            done = True; break
        if wall and time.time() - t0 > wall:
            break
    return kc.normalize(psi.detach()), hist, done


# ------------------------------------------------------------------------------------------ kinetic metric, inertia
def mass_apply(phi, v, h, ks=1.0, g=None):
    g = k2.derivs_fb(phi, h) if g is None else g
    gp, gm = g
    a = ks + 0.5 * sum((gp[i] * gp[i]).sum(0) + (gm[i] * gm[i]).sum(0) for i in range(3))
    s = math.sqrt(0.5)
    G = [s * x for i in range(3) for x in (gp[i], gm[i])]
    U = [gi - phi * (phi * gi).sum(0, keepdim=True) for gi in G]
    v = kc.tangent(phi, v)
    return kc.tangent(phi, a * v - sum(u * (u * v).sum(0, keepdim=True) for u in U))


def generators(phi, h):
    """iso generators A_a phi = (0, e_a x pi) and spatial S_a phi = -(e_a x x).grad phi (4th-order central), about the
    baryon centroid. With these the hedgehog has S_a = -A_a, i.e. W = -Lambda (symmetry = equal angles)."""
    X, Y, Z = kc.grid(phi.shape[1:], h)
    b = kc.baryon_density(phi, h); tot = b.sum()
    c = [float((b * A).sum() / tot) for A in (X, Y, Z)]
    x = [X - c[0], Y - c[1], Z - c[2]]
    d = kc.derivs(phi, h)
    pi = phi[1:]
    A, S = [], []
    for a in range(3):
        ea = [0.0, 0.0, 0.0]; ea[a] = 1.0
        cr = torch.stack([ea[1] * pi[2] - ea[2] * pi[1], ea[2] * pi[0] - ea[0] * pi[2], ea[0] * pi[1] - ea[1] * pi[0]])
        A.append(torch.cat([torch.zeros_like(phi[:1]), cr]))
        u = [ea[1] * x[2] - ea[2] * x[1], ea[2] * x[0] - ea[0] * x[2], ea[0] * x[1] - ea[1] * x[0]]
        S.append(-(u[0] * d[0] + u[1] * d[1] + u[2] * d[2]))
    return A, S, c


def inertia(phi, h, ks=1.0):
    """U, V, W with T = 1/2 (w.U.w + 2 w.W.W' + W'.V.W') (w iso, W' spatial angular velocity)."""
    g = k2.derivs_fb(phi, h)
    A, S, c = generators(phi, h)
    MA = [mass_apply(phi, x, h, ks, g) for x in A]
    MS = [mass_apply(phi, x, h, ks, g) for x in S]
    dot = lambda p, q: float(2 * (p * q).sum() * h**3)
    U = np.array([[dot(A[i], MA[j]) for j in range(3)] for i in range(3)])
    V = np.array([[dot(S[i], MS[j]) for j in range(3)] for i in range(3)])
    W = np.array([[dot(A[i], MS[j]) for j in range(3)] for i in range(3)])
    return U, V, W, c


# ------------------------------------------------------------------------------------------ Hessian
def hvp_factory(phi, h, **kw):
    """K(xi): second variation of E on the sphere bundle, per site (divided by h^3); tangent in/out."""
    p0 = phi.detach().clone().requires_grad_(True)
    E = energy(p0, h, **kw)
    (g0,) = torch.autograd.grad(E, p0, create_graph=True)
    lam = (phi * g0.detach()).sum(0, keepdim=True)

    def K(xi):
        xi = kc.tangent(phi, xi)
        (hv,) = torch.autograd.grad(g0, p0, grad_outputs=xi, retain_graph=True)
        return kc.tangent(phi, hv - lam * xi) / h**3
    return K
