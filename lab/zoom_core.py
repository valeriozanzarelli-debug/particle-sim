"""zoom_core: ONE agglomerator for every system and every zoom level (learned coarse-graining / RG).

The same five operations are used, without code changes, for a spring chain, the 2-D Ising model and
bundles of photon balls (docs/ZOOM.md):

  block(state, k)   group neighbouring balls (k per side) into super-balls.  A super-ball is described by
                    a small FIXED set of descriptors chosen by the DIMENSION d of space, never by system:
                      n      number of elementary balls inside            (scalar)
                      q      total charges (mass, spin, photon count ...)  (c values)
                      X      centre (count-weighted, periodic-safe)        (d values)
                      V      centre velocity (momentum / mass)             (d values)
                      s2     internal 2nd moment  <|x - X|^2>              (scalar)
                      eint   internal kinetic energy (heat of the block)   (scalar)
                    If the charges live in the alphabet {-1,+1} (spins), the coarse charge is projected back
                    to the same alphabet by majority (sign of the total, random tie-break), so that every
                    level has the SAME state space: that is what makes the rule repeatable (a fractal).
  learn(data)       fit the super-balls' rule with ONE generic model (LocalModel below):
                      'force'    force matching   (dynamics: mean force on super-balls)
                      'transfer' regression of descriptor changes over one coarse step
                      'stat'     maximum pseudo-likelihood (equilibrium statistics; the exponential-family
                                 relative-entropy / moment-matching fit, local version)
  step(state)       evolve super-balls with the learned rule (Verlet / Langevin / heat bath / transfer).
  check(...)        spot-check: re-simulate a (random, small) region at the finer level, return the error,
                    flag the region if the error exceeds a tolerance (refinement hook).
  repeat            the caller simply calls block + learn again on the super-balls (level 2, 3, ...).

The generic model ("a flexible structure that lets balls take geometries") is a local cluster expansion:
    W = sum_{pairs i<j, r<rcut} q_i q_j phi(r_ij / h)  +  H sum_i q_i  +  C3 sum_triangles q q q
        + C4 sum_quads q q q q
with phi a piecewise-linear function (hat basis, knots every DR spacings, out to RCUT spacings), triangles
and quads = compact cliques of balls (all pair distances <= 1.5 h).  h = the level's ball spacing, so
the same hyper-parameters mean the same thing at every level.  The model is linear in its parameters:
the force is -grad of a pair potential, the descriptor change is a sum of pair terms, the
log-probability is W.  Hyper-parameters (RCUT, DR) are the same for every system and level.
"""
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.optimize import minimize

RCUT = 6.5      # interaction range of the generic model, in units of the level spacing h
DR = 0.25       # hat-basis knot spacing, in units of h
RCLIQUE = 1.5   # max pair distance inside a compact cluster (triangle / quad), units of h


# ----------------------------------------------------------------------------------------------- state
class State:
    """A set of balls.  pos (N,d); q (N,c) charges (column 0 = primary charge); n (N,) elementary counts;
    ref (N,d) reference (label) positions used to decide membership in blocks; vel (N,d) or None;
    mass (N,) or None (dynamics only); s2, eint (N,) internal descriptors; box (d,) periodic or None;
    h = spacing of this level; alphabet = (-1, 1) for spin-like charges."""

    def __init__(self, pos, q, h, box=None, n=None, ref=None, vel=None, mass=None, s2=None, eint=None,
                 alphabet=None):
        self.pos = np.asarray(pos, float).reshape(len(pos), -1)
        N, self.d = self.pos.shape
        q = np.asarray(q, float)
        self.q = q.reshape(N, -1)
        self.h = float(h)
        self.box = None if box is None else np.asarray(box, float).reshape(self.d)
        self.n = np.ones(N) if n is None else np.asarray(n, float)
        self.ref = self.pos.copy() if ref is None else np.asarray(ref, float).reshape(N, self.d)
        self.vel = None if vel is None else np.asarray(vel, float).reshape(N, self.d)
        self.mass = None if mass is None else np.asarray(mass, float)
        self.s2 = np.zeros(N) if s2 is None else np.asarray(s2, float)
        self.eint = np.zeros(N) if eint is None else np.asarray(eint, float)
        self.alphabet = alphabet

    @property
    def N(self):
        return self.pos.shape[0]

    def copy(self, **kw):
        a = dict(pos=self.pos.copy(), q=self.q.copy(), h=self.h, box=self.box, n=self.n.copy(),
                 ref=self.ref, vel=None if self.vel is None else self.vel.copy(), mass=self.mass,
                 s2=self.s2.copy(), eint=self.eint.copy(), alphabet=self.alphabet)
        a.update(kw)
        return State(**a)


def min_image(dx, box):
    if box is None:
        return dx
    return dx - box * np.round(dx / box)


# ----------------------------------------------------------------------------------------------- block
class Blocking:
    """Membership map fine -> coarse (decided once from the reference positions), reusable on any
    configuration with the same labels.  Generic in d."""

    def __init__(self, fine, k):
        self.k = k
        H = k * fine.h
        org = fine.ref.min(0) - 0.5 * fine.h
        cell = np.floor((fine.ref - org) / H + 1e-9).astype(np.int64)
        if fine.box is not None:
            nc = np.maximum(np.round(fine.box / H).astype(np.int64), 1)
            cell %= nc
        ids, inv = np.unique(cell, axis=0, return_inverse=True)
        self.inv = inv.ravel()
        first = np.full(ids.shape[0], -1)
        first[self.inv[::-1]] = np.arange(self.inv.size)[::-1]
        self.first = first
        self.M = ids.shape[0]
        self.H = H
        self.cell = ids

    def sum(self, arr):
        arr = np.asarray(arr, float)
        if arr.ndim == 1:
            return np.bincount(self.inv, arr, self.M)
        flat = arr.reshape(arr.shape[0], -1)
        out = np.stack([np.bincount(self.inv, flat[:, c], self.M) for c in range(flat.shape[1])], 1)
        return out.reshape((self.M,) + arr.shape[1:])


def block(state, k, rng=None, blocking=None):
    """Group neighbouring balls (k per side) into super-balls and compute the generic descriptors.
    Returns (coarse_state, blocking)."""
    B = blocking if blocking is not None else Blocking(state, k)
    w = state.n
    nI = B.sum(w)
    # centre: unwrap members relative to the first member of their block (periodic-safe)
    first = B.first
    anchor = state.pos[first]
    dx = min_image(state.pos - anchor[B.inv], state.box)
    X = anchor + B.sum(w[:, None] * dx) / nI[:, None]
    if state.box is not None:
        X = X % state.box
    dref = state.ref - state.ref[first][B.inv]
    R = state.ref[first] + B.sum(w[:, None] * dref) / nI[:, None]
    rel = min_image(state.pos - X[B.inv], state.box)
    s2 = B.sum(w * (state.s2 + (rel ** 2).sum(1))) / nI
    qI = B.sum(state.q)
    vel = mass = None
    eint = B.sum(state.eint)
    if state.vel is not None:
        m = state.mass if state.mass is not None else w
        mass = B.sum(m)
        vel = B.sum(m[:, None] * state.vel) / mass[:, None]
        eint = eint + B.sum(0.5 * m * ((state.vel - vel[B.inv]) ** 2).sum(1))
    if state.alphabet is not None:      # project back to the same alphabet: majority, random tie-break
        rng = rng or np.random.default_rng(0)
        tot = qI[:, 0]
        tie = rng.choice(np.array(state.alphabet, float), size=tot.size)
        qI = np.where(tot > 0, 1.0, np.where(tot < 0, -1.0, tie))[:, None] * np.ones_like(qI)
    return State(X, qI, B.H, box=state.box, n=nI, ref=R, vel=vel, mass=mass, s2=s2, eint=eint,
                 alphabet=state.alphabet), B


# ----------------------------------------------------------------------------------------------- model
class Geometry:
    """Neighbour pairs (with hat-basis weights) and compact cliques of one configuration of positions."""

    def __init__(self, st, rcut=RCUT, dr=DR, want_cliques=True):
        self.N, self.d = st.N, st.d
        h = st.h
        if st.box is not None:
            tree = cKDTree(np.mod(st.pos, st.box), boxsize=st.box)
        else:
            tree = cKDTree(st.pos)
        P = tree.query_pairs(rcut * h, output_type="ndarray")
        if P.size == 0:
            P = np.zeros((0, 2), np.int64)
        i, j = P[:, 0], P[:, 1]
        rv = min_image(st.pos[i] - st.pos[j], st.box)
        r = np.sqrt((rv ** 2).sum(1))
        self.i, self.j, self.rv, self.r = i, j, rv, r
        self.u = r / h                                    # distance in level units
        self.knots = np.arange(dr, rcut + 1e-9, dr)
        self.dr = dr
        self.nb = self.knots.size
        # hat basis: each pair touches two knots (lower / upper) with linear weights
        t = (self.u - self.knots[0]) / dr
        lo = np.floor(t).astype(np.int64)
        f = t - lo
        self.lo, self.f = lo, f
        self.h = h
        self.cl3 = self.cl4 = None
        if want_cliques:
            self._cliques()

    def hats(self):
        """(P, nb) dense hat values (pairs below the first knot are clipped to it)."""
        G = np.zeros((self.r.size, self.nb))
        lo, f = self.lo, self.f
        c = lo < 0                                        # closer than first knot: flat continuation
        a = (~c) & (lo < self.nb)
        G[np.where(a)[0], lo[a]] += 1 - f[a]
        b = (~c) & (lo + 1 < self.nb)
        G[np.where(b)[0], lo[b] + 1] += f[b]
        G[np.where(c)[0], 0] = 1.0
        return G

    def hat_int(self):
        """(P, nb) integral of each hat from u to +inf (for the potential of a pair force)."""
        u = self.u[:, None]
        k = self.knots[None, :]
        dr = self.dr
        x = (u - k) / dr                                  # position relative to the knot, in knot units
        full = dr                                         # area of a whole hat
        I = np.where(x <= -1, full,
            np.where(x <= 0, full - 0.5 * dr * (1 + x) ** 2,
            np.where(x <= 1, 0.5 * dr * (1 - x) ** 2, 0.0)))
        return I

    def _cliques(self):
        close = self.u <= RCLIQUE + 1e-9
        ii, jj = self.i[close], self.j[close]
        nbr = [set() for _ in range(self.N)]
        for a, b in zip(ii.tolist(), jj.tolist()):
            nbr[a].add(b); nbr[b].add(a)
        c3, c4 = [], []
        for a in range(self.N):
            hi = sorted(x for x in nbr[a] if x > a)
            for x1 in range(len(hi)):
                b = hi[x1]
                for x2 in range(x1 + 1, len(hi)):
                    c = hi[x2]
                    if c in nbr[b]:
                        c3.append((a, b, c))
                        for x3 in range(x2 + 1, len(hi)):
                            e = hi[x3]
                            if e in nbr[b] and e in nbr[c]:
                                c4.append((a, b, c, e))
        self.cl3 = np.array(c3, np.int64).reshape(-1, 3)
        self.cl4 = np.array(c4, np.int64).reshape(-1, 4)


class LocalModel:
    """The ONE generic model.  theta = [phi hat weights (nb), H, C3, C4]."""

    def __init__(self, rcut=RCUT, dr=DR):
        self.rcut, self.dr = rcut, dr
        self.theta = None
        self.mode = None
        self.qscale = 1.0

    # ---- feature builders (all linear in theta)
    def _w(self, st, g):
        qt = st.q[:, 0] / self.qscale
        return qt[g.i] * qt[g.j]

    def vec_features(self, st, g):
        """(N*d, nb): pair-force features, F_i = sum_b theta_b feat."""
        G = g.hats() * self._w(st, g)[:, None]
        e = g.rv / np.maximum(g.r, 1e-12)[:, None]        # unit vector j -> i
        out = np.zeros((st.N, st.d, g.nb))
        for a in range(st.d):
            Ga = G * e[:, a:a + 1]
            np.add.at(out[:, a, :], g.i, Ga)
            np.add.at(out[:, a, :], g.j, -Ga)
        return out.reshape(st.N * st.d, g.nb)

    def sca_features(self, st, g):
        """(N, nb+1): [1, sum_j w g_b(r_ij)]."""
        G = g.hats() * self._w(st, g)[:, None]
        out = np.zeros((st.N, g.nb))
        np.add.at(out, g.i, G); np.add.at(out, g.j, G)
        return np.hstack([np.ones((st.N, 1)), out])

    def field_features(self, st, g):
        """(N, nb+3): dW/dq_i per parameter, for the +-1 alphabet (pseudo-likelihood / heat bath)."""
        q = st.q[:, 0]
        N, nb = st.N, g.nb
        out = np.zeros((N, nb + 3))
        lo = np.clip(g.lo, -1, nb - 1)
        f = np.where(g.lo < 0, 0.0, g.f)
        lo0 = np.maximum(lo, 0)
        hi = lo + 1
        okh = (hi < nb) & (g.lo >= 0)
        acc = np.zeros(N * nb)
        for a, b in ((g.i, g.j), (g.j, g.i)):
            acc += np.bincount(a * nb + lo0, (1 - f) * q[b], N * nb)
            acc += np.bincount((a * nb + hi)[okh], (f * q[b])[okh], N * nb)
        out[:, :nb] = acc.reshape(N, nb)
        out[:, nb] = 1.0
        for col, c in ((nb + 1, g.cl3), (nb + 2, g.cl4)):
            if c is None or len(c) == 0:
                continue
            P = np.prod(q[c], axis=1)
            for a in range(c.shape[1]):
                out[:, col] += np.bincount(c[:, a], P * q[c[:, a]], N)
        return out

    # ---- learn
    def learn_force(self, states, forces, alpha=1e-8):
        self.mode = "force"
        self.qscale = np.mean(np.abs(states[0].q[:, 0]))
        A, b = [], []
        for st, F in zip(states, forces):
            g = Geometry(st, self.rcut, self.dr, want_cliques=False)
            A.append(self.vec_features(st, g)); b.append(np.asarray(F).ravel())
        A = np.vstack(A); b = np.concatenate(b)
        keep = np.abs(A).sum(0) > 0
        th = np.zeros(A.shape[1])
        AtA = A[:, keep].T @ A[:, keep]
        lam = alpha * np.trace(AtA) / max(keep.sum(), 1)
        th[keep] = np.linalg.solve(AtA + lam * np.eye(keep.sum()), A[:, keep].T @ b)
        self.theta = th
        res = b - A @ th
        return dict(rel_residual=float(np.sqrt((res ** 2).mean() / (b ** 2).mean())), nparam=int(keep.sum()))

    def learn_transfer(self, states, dX, ds2, alpha=1e-6):
        """Descriptor change over one coarse step: dX = sum of pair vector terms (in units of h);
        the positive internal moment changes by a factor: log(s2_out/s2_in) = bias + pair terms
        (log, because positive descriptors change by factors and their tails would dominate a
        squared loss)."""
        self.mode = "transfer"
        self.qscale = np.mean(np.abs(states[0].q[:, 0]))
        Av, bv, As, bs = [], [], [], []
        for st, a, b in zip(states, dX, ds2):
            g = Geometry(st, self.rcut, self.dr, want_cliques=False)
            Av.append(self.vec_features(st, g)); bv.append(np.asarray(a).ravel() / st.h)
            tiny = 1e-12 * st.h ** 2
            As.append(self.sca_features(st, g)); bs.append(np.log((st.s2 + np.asarray(b) + tiny) / (st.s2 + tiny)))
        fits = []
        for A, b in ((np.vstack(Av), np.concatenate(bv)), (np.vstack(As), np.concatenate(bs))):
            keep = np.abs(A).sum(0) > 0
            th = np.zeros(A.shape[1])
            if keep.any():
                AtA = A[:, keep].T @ A[:, keep]
                lam = alpha * np.trace(AtA) / keep.sum() + 1e-12
                th[keep] = np.linalg.solve(AtA + lam * np.eye(keep.sum()), A[:, keep].T @ b)
            fits.append(th)
        self.theta_v, self.theta_s = fits
        return dict(nparam=int((self.theta_v != 0).sum() + (self.theta_s != 0).sum()))

    def learn_stat(self, states, geom=None, alpha=1e-4, sites=None, rng=None):
        """Maximum pseudo-likelihood for +-1 charges: P(q_i | rest) = sigmoid(2 q_i f_i . theta)."""
        self.mode = "stat"
        g = geom or Geometry(states[0], self.rcut, self.dr)
        rng = rng or np.random.default_rng(0)
        D, Q = [], []
        for st in states:
            F = self.field_features(st, g)
            idx = np.arange(st.N) if sites is None or sites >= st.N else rng.choice(st.N, sites, replace=False)
            D.append(F[idx]); Q.append(st.q[idx, 0])
        D = np.vstack(D); Q = np.concatenate(Q)
        keep = np.abs(D).sum(0) > 0
        Dk = D[:, keep] * Q[:, None]
        n = Dk.shape[0]

        def f(th):
            z = 2 * Dk @ th
            ll = -np.logaddexp(0, -z).sum() / n
            gr = (2 * Dk * (1 / (1 + np.exp(z)))[:, None]).sum(0) / n
            return -ll + alpha * th @ th, -gr + 2 * alpha * th

        r = minimize(f, np.zeros(keep.sum()), jac=True, method="L-BFGS-B", options=dict(maxiter=2000))
        th = np.zeros(D.shape[1]); th[keep] = r.x
        self.theta = th
        self.geom = g
        return dict(nll=float(r.fun), nparam=int(keep.sum()), ok=bool(r.success))

    # ---- evaluate
    def phi(self, u):
        """Learned pair function at distances u (level units)."""
        k = np.arange(self.dr, self.rcut + 1e-9, self.dr)
        th = self.theta if self.mode != "transfer" else self.theta_v
        return np.interp(u, k, th[:k.size], left=th[0], right=0.0)

    def dphi(self, u, eps=1e-4):
        return (self.phi(u + eps) - self.phi(u - eps)) / (2 * eps)

    def force(self, st, g=None):
        """Learned force (same numbers as vec_features @ theta, evaluated pair by pair)."""
        g = g or Geometry(st, self.rcut, self.dr, want_cliques=False)
        ph = self.phi(g.u) * self._w(st, g)
        fv = ph[:, None] * g.rv / np.maximum(g.r, 1e-12)[:, None]
        F = np.empty((st.N, st.d))
        for a in range(st.d):
            F[:, a] = np.bincount(g.i, fv[:, a], st.N) - np.bincount(g.j, fv[:, a], st.N)
        return F

    def energy(self, st, g=None):
        """Potential of the learned pair force: U = sum_pairs w * h * sum_b theta_b int_u^inf g_b."""
        g = g or Geometry(st, self.rcut, self.dr, want_cliques=False)
        return float((g.hat_int() @ self.theta * self._w(st, g)).sum() * st.h)

    def transfer(self, st):
        g = Geometry(st, self.rcut, self.dr, want_cliques=False)
        dX = (self.vec_features(st, g) @ self.theta_v).reshape(st.N, st.d) * st.h
        ds2 = st.s2 * (np.exp(self.sca_features(st, g) @ self.theta_s) - 1)
        return dX, ds2

    def couplings(self):
        """Readable statistical couplings: phi at the lattice shells is reported by the caller."""
        nb = self.geom.nb
        return dict(H=self.theta[nb], C3=self.theta[nb + 1], C4=self.theta[nb + 2])


# ----------------------------------------------------------------------------------------------- step
def verlet(model, st, dt, nsteps, record=None):
    """NVE velocity Verlet with the learned force (positions/velocities of super-balls)."""
    st = st.copy()
    m = (st.mass if st.mass is not None else st.n)[:, None]
    F = model.force(st)
    for t in range(nsteps):
        st.vel += 0.5 * dt * F / m
        st.pos += dt * st.vel
        F = model.force(st)
        st.vel += 0.5 * dt * F / m
        if record is not None:
            record(t, st, F)
    return st


def langevin(model_or_force, st, dt, nsteps, T, gamma, rng, every=0, out=None):
    """BAOAB Langevin; model_or_force = LocalModel or a callable force(state)."""
    st = st.copy()
    m = (st.mass if st.mass is not None else st.n)[:, None]
    force = model_or_force.force if hasattr(model_or_force, "force") else model_or_force
    c1 = np.exp(-gamma * dt); c2 = np.sqrt((1 - c1 ** 2) * T / m)
    F = force(st)
    for t in range(nsteps):
        st.vel += 0.5 * dt * F / m
        st.pos += 0.5 * dt * st.vel
        st.vel = c1 * st.vel + c2 * rng.standard_normal(st.vel.shape)
        st.pos += 0.5 * dt * st.vel
        F = force(st)
        st.vel += 0.5 * dt * F / m
        if every and (t + 1) % every == 0 and out is not None:
            out.append((st.copy(), F.copy()))
    return st


class HeatBath:
    """Heat-bath sweeps of the learned statistical model on a fixed lattice (+-1 charges).
    Colouring by reference cell index mod m, with m*h > interaction range, so same-colour balls do not
    interact (generic in d)."""

    def __init__(self, model, st):
        g = model.geom if getattr(model, "geom", None) is not None and model.geom.N == st.N else Geometry(st, model.rcut, model.dr)
        self.g = g
        nb = g.nb
        th = model.theta
        w = model.phi(g.u)
        J = sp.coo_matrix((np.r_[w, w], (np.r_[g.i, g.j], np.r_[g.j, g.i])), shape=(st.N, st.N)).tocsr()
        self.J, self.H, self.C3, self.C4 = J, th[nb], th[nb + 1], th[nb + 2]
        cell = np.round((st.ref - st.ref.min(0)) / st.h).astype(np.int64)
        ncell = cell.max(0) + 1
        need = int(np.ceil(model.rcut)) + 1
        col = np.zeros(st.N, np.int64)
        for a in range(st.d):
            m = next((x for x in range(need, ncell[a] + 1) if ncell[a] % x == 0), ncell[a])
            col = col * m + cell[:, a] % m
        self.colours = [np.where(col == c)[0] for c in np.unique(col)]
        self.Jc = [J[c] for c in self.colours]
        self.local = np.zeros(st.N, np.int64)
        colid = np.zeros(st.N, np.int64)
        for ci, c in enumerate(self.colours):
            self.local[c] = np.arange(c.size); colid[c] = ci
        self.cmap = []
        for ci in range(len(self.colours)):
            ent = []
            for cl in (g.cl3, g.cl4):
                if cl is None or len(cl) == 0:
                    ent.append((np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros((0, 1), np.int64)))
                    continue
                r, p = np.nonzero(colid[cl] == ci)
                ent.append((r, p, cl))
            self.cmap.append(ent)

    def field(self, q, ci):
        idx = self.colours[ci]
        f = self.Jc[ci] @ q + self.H
        for C, (rows, pos, cl) in zip((self.C3, self.C4), self.cmap[ci]):
            if C == 0 or rows.size == 0:
                continue
            P = np.prod(q[cl[rows]], axis=1) * q[cl[rows, pos]]
            f = f + C * np.bincount(self.local[cl[rows, pos]], P, idx.size)
        return f

    def sweep(self, q, rng, n=1):
        q = q.copy()
        for _ in range(n):
            for ci, idx in enumerate(self.colours):
                f = self.field(q, ci)
                p = 1 / (1 + np.exp(-2 * f))
                q[idx] = np.where(rng.random(idx.size) < p, 1.0, -1.0)
        return q


# ----------------------------------------------------------------------------------------------- check
def check(region, fine_fn, coarse_fn, err_fn, tol):
    """Spot-check / refinement hook: re-simulate `region` at the finer level (fine_fn), simulate it with
    the coarse rule (coarse_fn), compare with err_fn.  Returns dict(err, flagged): a flagged region must
    be refined (kept at the finer level)."""
    a = fine_fn(region)
    b = coarse_fn(region)
    e = float(err_fn(a, b))
    return dict(region=region, err=e, flagged=bool(e > tol))
