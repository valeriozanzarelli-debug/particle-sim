"""nucleo1s_core.py - the zig-zag radial rule on a GRADED radial mesh (cells fine near the nucleus, coarse far away),
and the tools of the study "1s pulito vicino al nucleo" (docs/NUCLEO_1S.md).

THE RULE IS THE SAME as lab/zz_core.py RadialWalk (frozen 2026-10-09): two radial walkers per site, chi+ (moves out at c)
and chi- (moves in at c); at every site the rotation exp(-i dt (m sigma_z + kappa/r sigma_x + V(r))) on (g, f); at r = 0
the inward walker turns into the outward one (g(0) = 0). Only one thing changes: the cells need not all be equal.
Cell i is [e_i, e_{i+1}], width h_i, site at its centre r_i; the rotation of site i uses dt = h_i (the time a ball at c
spends crossing the cell); the walker going from site i to site i+1 travels the distance l = r_{i+1} - r_i at c.

Stationary states ("listening for an infinite time", as in ZIGZAG.md): a state with frequency omega picks up the phase
e^{i omega l} along every link. With x = post-rotation amplitudes (chi+_i, chi-_i):
        x_i = M_i ( e^{i w l_{i-1,i}} chi+_{i-1} , e^{i w l_{i,i+1}} chi-_{i+1} )
i.e. x = U(w) x with U(w) = M S Phi(w) unitary. For equal cells U(w) = e^{i w a} W and this is exactly RadialWalk.eigen.
For unequal cells the same equation is the stationary limit of a unitary dynamical rule on the finest tick in which the
coarse zone has "delay slots" (a ball crosses a coarse cell in several fine ticks and is rotated once, by the whole
coarse angle): the agglomerator's "blocca" applied to the rotations. (Live dynamics on nested grids is NOT done here.)

Units hbar = c = m_e = 1 (length lambda_C, energy m c^2). Inputs: m_e, alpha, Z.  Dirac formula = CHECK ONLY.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA, MC2_EV, MP_ME, U_CHI, dirac_formula, RadialWalk  # noqa: F401


def rot_chi(dt, kappa, r, Za, mass=1.0, vmode="point", edges=None):
    """per-site rotation in the chi basis, exactly as zz_core.RadialWalk (vmode='point': V = -Za/r at the site)."""
    hx = kappa / r
    hz = mass * np.ones_like(r)
    V = -Za / r
    nrm = np.sqrt(hx**2 + hz**2)
    phi = dt * nrm
    c, s = np.cos(phi), np.sin(phi)
    M = np.zeros((len(r), 2, 2), complex)
    M[:, 0, 0] = c - 1j * s * hz / nrm
    M[:, 1, 1] = c + 1j * s * hz / nrm
    M[:, 0, 1] = -1j * s * hx / nrm
    M[:, 1, 0] = -1j * s * hx / nrm
    M *= np.exp(-1j * V * dt)[:, None, None]
    return np.einsum("ab,nbc,cd->nad", U_CHI.conj().T, M, U_CHI)


class Mesh:
    """a radial mesh given by its cell edges e_0 = 0 < e_1 < ... < e_N"""

    def __init__(self, edges):
        self.e = np.asarray(edges, float)
        self.h = np.diff(self.e)
        self.r = 0.5 * (self.e[1:] + self.e[:-1])
        self.N = len(self.r)
        # link lengths: inner reflection (site 0 -> origin -> site 0) = h_0; i -> i+1 = r_{i+1} - r_i; outer = h_{N-1}
        self.l_in = 2 * self.r[0]
        self.l = np.diff(self.r)
        self.l_out = self.h[-1]

    @staticmethod
    def uniform(a, rmax):
        n = int(np.ceil(rmax / a))
        return Mesh(np.arange(n + 1) * a)

    @staticmethod
    def graded(a_coarse, rmax, zones):
        """zones = [(R_k, b_k), ...]: inside radius R_k the cell is a_coarse / b_k (R_k multiples of a_coarse,
        decreasing R with increasing b). Edges are nested: every coarse cell is split into b equal fine cells."""
        n = int(np.ceil(rmax / a_coarse))
        edges = [0.0]
        for i in range(n):
            lo = i * a_coarse
            b = 1
            for (R, bk) in zones:
                if lo < R - 1e-12:
                    b = max(b, bk)
            for q in range(1, b + 1):
                edges.append(lo + q * a_coarse / b)
        return Mesh(np.array(edges))


class GradedWalk:
    def __init__(self, mesh, kappa, Za, mass=1.0, refl0=-1.0, refl1=-1.0, M=None):
        self.m, self.kappa, self.Za, self.mass = mesh, kappa, Za, mass
        self.refl0, self.refl1 = refl0, refl1
        self.M = rot_chi(mesh.h, kappa, mesh.r, Za, mass) if M is None else M
        N = mesh.N
        # S Phi(w): incoming chi+ at i comes from x[2(i-1)] over link l[i-1] (i>=1) or from refl0 * x[1] over l_in;
        #           incoming chi- at i comes from x[2(i+1)+1] over l[i] (i<=N-2) or from refl1 * x[2(N-1)] over l_out
        i = np.arange(N)
        self.src_p = np.where(i >= 1, 2 * (i - 1), 1)
        self.fac_p = np.where(i >= 1, 1.0, refl0)
        self.len_p = np.where(i >= 1, np.concatenate([[0.0], mesh.l]), mesh.l_in)
        self.src_q = np.where(i <= N - 2, 2 * (i + 1) + 1, 2 * (N - 1))
        self.fac_q = np.where(i <= N - 2, 1.0, refl1)
        self.len_q = np.where(i <= N - 2, np.concatenate([mesh.l, [0.0]]), mesh.l_out)

    def U(self, w):
        N = self.m.N
        i = np.arange(N)
        ph_p = self.fac_p * np.exp(1j * w * self.len_p)
        ph_q = self.fac_q * np.exp(1j * w * self.len_q)
        rows, cols, vals = [], [], []
        for s in (0, 1):
            rows += [2 * i + s, 2 * i + s]
            cols += [self.src_p, self.src_q]
            vals += [self.M[:, s, 0] * ph_p, self.M[:, s, 1] * ph_q]
        return sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(2 * N, 2 * N))

    def link_weight(self, x):
        """d arg(mu)/d w for the eigenvector x of U(w): sum over incoming links of |incoming amplitude|^2 * length"""
        p = np.abs(x[self.src_p])**2 * self.len_p
        q = np.abs(x[self.src_q])**2 * self.len_q
        return (np.sum(p) + np.sum(q)) / np.sum(np.abs(x)**2)

    def eigen(self, w0, tol=1e-14, maxit=40, inner=6, x0=None):
        """frequency w with x = U(w) x, by shift-invert on U(w) at 1 plus Newton on arg mu(w)."""
        w = w0
        x = np.ones(2 * self.m.N, complex) if x0 is None else x0.copy()
        I = sp.identity(2 * self.m.N, format="csc")
        for it in range(maxit):
            Uw = self.U(w)
            lu = spl.splu((Uw - I).tocsc())
            for _ in range(inner):
                y = lu.solve(x)
                x = y / np.linalg.norm(y)
            mu = np.vdot(x, Uw @ x)
            dw = -np.angle(mu) / self.link_weight(x)
            w += dw
            if abs(dw) < tol * max(1.0, abs(w)):
                break
        Uw = self.U(w)
        mu = np.vdot(x, Uw @ x)
        resid = np.linalg.norm(Uw @ x - x)
        return w, x, dict(iters=it + 1, resid=resid, mu_arg=float(np.angle(mu)))

    def gf(self, x):
        """(g, f) at the sites (palindromic point not applied), normalised with the cell widths"""
        chi = x.reshape(self.m.N, 2).T
        g = (chi[0] + chi[1]) / np.sqrt(2)
        f = 1j * (chi[0] - chi[1]) / np.sqrt(2)
        k = np.argmax(np.abs(g))
        ph = np.exp(-1j * np.angle(g[k]))
        g, f = g * ph, f * ph
        nrm = np.sqrt(np.sum((np.abs(g)**2 + np.abs(f)**2) * self.m.h))
        return g / nrm, f / nrm


def dirac_binding(n, kappa, Za):
    return 1.0 - dirac_formula(n, kappa, Za)


def level(mesh, n, kappa, Za, guess_rel=1.05, **kw):
    """binding of level (n, kappa) on a mesh: listen near guess_rel x the Dirac binding (the Dirac number is only the
    place where we start listening; the rule decides where the level is). Returns (binding, info)."""
    ref = dirac_binding(n, kappa, Za)
    gw = GradedWalk(mesh, kappa, Za, **kw)
    w, x, info = gw.eigen(1.0 - guess_rel * ref)
    return 1.0 - w, dict(info, N=mesh.N, gw=gw, x=x)


def gamma_of(kappa, Za):
    return np.sqrt(kappa**2 - Za**2)


def mesh_zoom(a_out, rmax, levels, n_per_level, n_inner=None):
    """the zoom mesh: cells a_out far away; going inward, `levels` nested levels with cells a_out/2^k (k = 1..levels),
    each level `n_per_level` cells thick; the innermost level (k = levels) reaches r = 0 (n_inner cells if given).
    Nested: every level boundary is an edge of the coarser level (2 fine cells per coarse cell)."""
    # build from the inside out
    edges = [0.0]
    h = a_out / 2**levels
    n_in = n_inner if n_inner is not None else n_per_level
    for _ in range(n_in):
        edges.append(edges[-1] + h)
    for k in range(levels - 1, 0, -1):
        h = a_out / 2**k
        # align to the coarser grid: pad with fine cells until edges[-1] is a multiple of h
        while abs(edges[-1] / h - round(edges[-1] / h)) > 1e-9:
            edges.append(edges[-1] + h / 2)
        for _ in range(n_per_level):
            edges.append(edges[-1] + h)
    while abs(edges[-1] / a_out - round(edges[-1] / a_out)) > 1e-9:
        edges.append(edges[-1] + a_out / 2)
    while edges[-1] < rmax - 1e-12:
        edges.append(edges[-1] + a_out)
    return Mesh(np.array(edges))


# ----------------------------------------------------------------------------------------------------------------
# automatic zoom ("controlla e ripeti"): split the cells whose splitting changes the observed level the most
# ----------------------------------------------------------------------------------------------------------------
def split_cells(edges, which):
    """split the cells with indices in `which` into two equal halves"""
    e = list(edges)
    out = [e[0]]
    w = set(int(i) for i in which)
    for i in range(len(e) - 1):
        if i in w:
            out.append(0.5 * (e[i] + e[i + 1]))
        out.append(e[i + 1])
    return np.array(out)


def balance(edges, ratio=2.0):
    """2:1 balance: neighbouring cells differ by at most `ratio` (split the bigger one until true)"""
    e = np.asarray(edges, float)
    while True:
        h = np.diff(e)
        bad = set()
        for i in range(len(h) - 1):
            if h[i] > ratio * h[i + 1] * (1 + 1e-9):
                bad.add(i)
            elif h[i + 1] > ratio * h[i] * (1 + 1e-9):
                bad.add(i + 1)
        if not bad:
            return e
        e = split_cells(e, sorted(bad))


class AutoZoom:
    """start from equal cells a0 (the lab lattice); repeat: for every cell, split it in two and measure how much the
    observed level(s) move (eta); estimate the remaining error of every cell with the convergence law
    (eps_c = eta_c / (1 - 2^-p), p = 2 in the bulk, p = 2 gamma for the cell touching the nucleus); stop when the sum
    is below the budget; otherwise split the cells that move the level the most (all with eta >= frac * max)."""

    def __init__(self, a0, rmax, states, Za, budget, frac=0.3, max_iter=40, ratio=2.0, log=print):
        self.a0, self.rmax, self.states, self.Za = a0, rmax, states, Za
        self.budget, self.frac, self.max_iter, self.ratio, self.log = budget, frac, max_iter, ratio, log
        self.edges = Mesh.uniform(a0, rmax).e
        self.history = []

    def levels(self, edges, guess=None):
        out = []
        for j, (n, kappa) in enumerate(self.states):
            ref = dirac_binding(n, kappa, self.Za)
            g = guess[j] / ref if guess is not None else (1.3 if kappa == -1 else 1.05)
            b, _ = level(Mesh(edges), n, kappa, self.Za, guess_rel=g)
            out.append(b)
        return np.array(out)

    def run(self):
        E0 = self.levels(self.edges)
        for it in range(self.max_iter):
            N = len(self.edges) - 1
            eta = np.zeros((N, len(self.states)))
            for c in range(N):
                e2 = balance(split_cells(self.edges, [c]), self.ratio)
                eta[c] = (self.levels(e2, guess=E0) - E0) / E0          # relative change of each observed level
            g = gamma_of(-1, self.Za)
            p = np.full(N, 2.0)
            p[0] = min(2.0, 2 * g)
            eps = eta / (1 - 2.0**(-p))[:, None]
            est = np.abs(eps.sum(axis=0))                               # estimated remaining relative error per level
            score = np.abs(eta).max(axis=1)
            rec = dict(iter=it, cells=N, levels=E0.tolist(), est_err=est.tolist(), max_eta=float(score.max()),
                       hmin=float(np.diff(self.edges).min()))
            self.history.append(rec)
            self.log(f"  zoom it {it}: cells {N}, h_min {rec['hmin']:.3e}, levels {E0}, estimated error {est}")
            if est.max() < self.budget:
                break
            pick = np.where(score >= self.frac * score.max())[0]
            self.edges = balance(split_cells(self.edges, pick), self.ratio)
            E0 = self.levels(self.edges, guess=E0)
        self.E = E0
        return self


# ----------------------------------------------------------------------------------------------------------------
# the nucleus "agglomerated": the fine cells inside r < R replaced by what they do to a ball that comes in at R and
# goes out again (a reflection with a phase that depends on the frequency).  Learned from the fine rule itself.
# ----------------------------------------------------------------------------------------------------------------
def inner_reflection(edges_in, kappa, Za, w, a_out, refl0=-1.0):
    """rho(w) = (chi+ arriving at the first coarse site) / (chi- leaving it), for the fine zig-zag cells edges_in
    (edges_in[0] = 0, edges_in[-1] = R) and a first coarse cell of width a_out just outside R."""
    m = Mesh(edges_in)
    gw = GradedWalk(m, kappa, Za, refl0=refl0, refl1=0.0)      # outer end: no reflection, a source instead
    N = m.N
    U = gw.U(w)
    lb = 0.5 * (m.h[-1] + a_out)
    src = np.zeros(2 * N, complex)
    # chi- arriving at the last inner site from outside: x_{N-1} += M_{N-1}[:, 1] * e^{i w lb} * 1
    src[2 * (N - 1)] = gw.M[N - 1, 0, 1] * np.exp(1j * w * lb)
    src[2 * (N - 1) + 1] = gw.M[N - 1, 1, 1] * np.exp(1j * w * lb)
    x = spl.spsolve((sp.identity(2 * N, format="csc") - U).tocsc(), src)
    return np.exp(1j * w * lb) * x[2 * (N - 1)]


class OuterWalk(GradedWalk):
    """the coarse walk outside R with the nucleus agglomerated into a reflection rho(w) (callable)"""

    def __init__(self, mesh, kappa, Za, rho, **kw):
        super().__init__(mesh, kappa, Za, **kw)
        self.rho = rho
        self.len_p = self.len_p.copy()
        self.len_p[0] = 0.0                       # the phase of the inner trip is inside rho

    def U(self, w):
        self.fac_p = self.fac_p.astype(complex)
        self.fac_p[0] = self.rho(w)
        return super().U(w)

    def link_weight(self, x):
        # d arg mu / d w: links + the delay of the inner trip (numerical derivative of arg rho)
        dw = 1e-6
        tau = np.angle(self.rho(self._w + dw) / self.rho(self._w - dw)) / (2 * dw)
        p = np.abs(x[self.src_p])**2 * self.len_p
        p[0] = np.abs(x[self.src_p[0]])**2 * tau
        q = np.abs(x[self.src_q])**2 * self.len_q
        return (np.sum(p) + np.sum(q)) / np.sum(np.abs(x)**2)

    def eigen(self, w0, **kw):
        self._w = w0
        orig_U = self.U

        def U_track(w):
            self._w = w
            return orig_U(w)
        self.U = U_track
        try:
            return super().eigen(w0, **kw)
        finally:
            self.U = orig_U
