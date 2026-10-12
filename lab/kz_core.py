"""kz_core.py - engine for the 'sincronia' (Kibble-Zurek) tests.

Classical XY model (planar rotors, 'metronomes'): H = -J sum_<ij> cos(theta_i - theta_j), J = 1, simple cubic / square
lattice, periodic boundaries.  Dynamics:
  * metropolis(): checkerboard Metropolis, proposal theta -> theta + U(-D, D), D = 1 rad (model A, no conservation).
    One 'sweep' = both colours = 1 time unit.  Exact equilibrium at temperature T (no time-step bias).
  * sw_update(): Swendsen-Wang embedded-Ising cluster update (Wolff 1989 reflection trick), for equilibrium only (K2).
  * inertial(): second-order (Hamiltonian) dynamics  theta'' = -gamma theta' + sum_j sin(theta_j - theta_i)  (K3).
  * vortex counting by plaquette winding; gradient-flow 'cooling' (theta' = sum_j sin(theta_j - theta_i), T = 0).
torch on CPU, <= 2 threads.
"""
import math, numpy as np, torch

torch.set_num_threads(2)
TC3 = 2.20184     # 3D XY simple cubic, J = 1 (Campostrini et al. 2006, beta_c = 0.454165)
NU_XY = 0.6717    # 3D XY nu (Campostrini 2006; helium 0.6709)
TWO_PI = 2 * math.pi


def wrap(a):
    return torch.remainder(a + math.pi, TWO_PI) - math.pi


def dims(th):
    return tuple(range(th.dim()))


def local_field(th):
    c, s = torch.cos(th), torch.sin(th)
    hx = torch.zeros_like(th); hy = torch.zeros_like(th)
    for d in dims(th):
        for sh in (1, -1):
            hx += torch.roll(c, sh, d); hy += torch.roll(s, sh, d)
    return hx, hy


def parity_mask(shape):
    idx = torch.zeros(shape, dtype=torch.int64)
    for d, n in enumerate(shape):
        v = torch.arange(n).view([-1 if k == d else 1 for k in range(len(shape))])
        idx = idx + v
    return (idx % 2 == 0)


def metropolis(th, T, gen, masks, D=1.0):
    """one sweep (two colours) at temperature T; returns acceptance."""
    acc = 0.0
    for m in masks:
        hx, hy = local_field(th)
        new = th + (torch.rand(th.shape, generator=gen, dtype=th.dtype) * 2 - 1) * D
        dE = -(hx * (torch.cos(new) - torch.cos(th)) + hy * (torch.sin(new) - torch.sin(th)))
        if T <= 0:
            ok = dE <= 0
        else:
            ok = torch.rand(th.shape, generator=gen, dtype=th.dtype) < torch.exp(torch.clamp(-dE / T, max=0.0))
        ok = ok & m
        th = torch.where(ok, wrap(new), th)
        acc += ok.sum().item()
    return th, acc / th.numel()


def energy(th):
    e = 0.0
    for d in dims(th):
        e -= torch.cos(torch.roll(th, -1, d) - th).sum().item()
    return e


def vortex_count(th):
    """(sum |w|, number of w>0) over all plaquettes, all orientations.  3D: vortex-line length in lattice units;
    2D: number of point vortices (+ and -)."""
    D = th.dim()
    diffs = [wrap(torch.roll(th, -1, d) - th) for d in range(D)]
    tot = 0
    pos = 0
    for a in range(D):
        for b in range(a + 1, D):
            circ = diffs[a] + torch.roll(diffs[b], -1, a) - torch.roll(diffs[a], -1, b) - diffs[b]
            w = torch.round(circ / TWO_PI)
            tot += w.abs().sum().item()
            pos += (w > 0).sum().item()
    return tot, pos


def vortex_map_2d(th):
    dx = wrap(torch.roll(th, -1, 0) - th); dy = wrap(torch.roll(th, -1, 1) - th)
    circ = dx + torch.roll(dy, -1, 0) - torch.roll(dx, -1, 1) - dy
    return torch.round(circ / TWO_PI)


def force(th):
    c, s = torch.cos(th), torch.sin(th)
    hx, hy = local_field(th)
    return hy * c - hx * s     # sum_j sin(th_j - th_i)


def cool(th, tflow, dt=0.1):
    """gradient flow at T = 0 for time tflow (removes small thermal loops)."""
    n = int(round(tflow / dt))
    for _ in range(n):
        th = th + dt * force(th)
    return wrap(th)


def inertial(th, p, steps, dt=0.05, gamma=0.0):
    """velocity-Verlet (damping as half-step exponential) for theta'' = -gamma theta' + F."""
    f = force(th)
    damp = math.exp(-gamma * dt / 2)
    for _ in range(steps):
        p = p * damp + 0.5 * dt * f
        th = th + dt * p
        f = force(th)
        p = (p + 0.5 * dt * f) * damp
    return th, p


def kin_energy(p):
    return 0.5 * (p * p).sum().item()


# ---------------- Swendsen-Wang embedded cluster (equilibrium only) ----------------
def sw_update(th_np, beta, rng):
    """one Swendsen-Wang sweep with the Wolff reflection embedding. th_np: numpy array (any dim)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    phi = rng.uniform(0, TWO_PI)
    s = np.cos(th_np - phi)
    N = th_np.size
    idx = np.arange(N).reshape(th_np.shape)
    rows, cols = [], []
    for d in range(th_np.ndim):
        sn = np.roll(s, -1, d)
        p = 1.0 - np.exp(np.minimum(0.0, -2.0 * beta * s * sn))
        bond = rng.random(th_np.shape) < p
        rows.append(idx[bond]); cols.append(np.roll(idx, -1, d)[bond])
    r = np.concatenate(rows); c = np.concatenate(cols)
    g = coo_matrix((np.ones(len(r), dtype=np.int8), (r, c)), shape=(N, N))
    nc, lab = connected_components(g, directed=False)
    flip = rng.random(nc) < 0.5
    f = flip[lab].reshape(th_np.shape)
    # reflection about the line perpendicular to phi: theta -> pi + 2 phi - theta
    th_new = np.where(f, np.pi + 2 * phi - th_np, th_np)
    return np.mod(th_new + np.pi, TWO_PI) - np.pi


def overrelax_np(th_np, par):
    """microcanonical over-relaxation (checkerboard): reflect each spin about its local field (energy preserving)."""
    for col in (0, 1):
        c, s = np.cos(th_np), np.sin(th_np)
        hx = sum(np.roll(c, sh, d) for d in range(th_np.ndim) for sh in (1, -1))
        hy = sum(np.roll(s, sh, d) for d in range(th_np.ndim) for sh in (1, -1))
        ang = np.arctan2(hy, hx)
        th_np = np.where(par == col, 2 * ang - th_np, th_np)
    return np.mod(th_np + np.pi, TWO_PI) - np.pi


# ---------------- linear quench ----------------
def quench(L, D, tauQ, count_times, seed, Tc=TC3, Ti_rel=1.2, n_eq=200, cool_times=(0.0,), dtype=torch.float32,
           keep_final=False, cool_dt=0.15):
    """T(t) = Tc (1 - t/tauQ) (clamped at 0), start at T = Ti_rel*Tc after n_eq sweeps there; t in sweeps.
    count_times: times t after the Tc crossing at which vortices are counted.  At each count time a COPY of the field
    is cooled by gradient flow (T = 0), incrementally, and counted after each cooling time in cool_times (a list, or a
    dict {t_count: list}); the quench continues on the uncooled field.
    returns {(t_count, t_cool): (sum|w|, n_pos)} and optionally the final (uncooled) configuration."""
    gen = torch.Generator().manual_seed(seed)
    shape = (L,) * D
    masks = [parity_mask(shape), ~parity_mask(shape)]
    th = (torch.rand(shape, generator=gen, dtype=dtype) * 2 - 1) * math.pi
    Ti = Ti_rel * Tc
    for _ in range(n_eq):
        th, _ = metropolis(th, Ti, gen, masks)
    t = -float(int(math.ceil((Ti_rel - 1.0) * tauQ)))
    res = {}
    targets = sorted(count_times)
    k = 0
    while True:
        while k < len(targets) and t >= targets[k] - 1e-9:
            tk = targets[k]
            cts = sorted(cool_times[tk] if isinstance(cool_times, dict) else cool_times)
            thc = th.clone(); done = 0.0
            for tc in cts:
                if tc > done:
                    n = int(round((tc - done) / cool_dt))
                    for _ in range(n):
                        thc = thc + cool_dt * force(thc)
                    thc = wrap(thc); done = tc
                res[(tk, tc)] = vortex_count(thc)
            k += 1
        if k >= len(targets):
            break
        T = max(Tc * (1.0 - (t + 0.5) / tauQ), 0.0)
        th, _ = metropolis(th, T, gen, masks)
        t += 1.0
    return (res, th) if keep_final else res
