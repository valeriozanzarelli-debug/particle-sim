"""big_core.py - shared tools for the "big step" (self-built geometry) and its two repairs R1, R2.

Imports pre_core / geo_core (unchanged). Contents:
  1. REPAIRED SHAPE READER (R1).  Blind, from the mutual-information matrix I only:
       largest informational piece (geo_core.components) -> d = -log(I/I0) -> Isomap geodesic G (geo_core)
       -> diameter D = max G -> max-min landmarks (<= 220, geo_core.landmarks) -> Vietoris-Rips H1 with ripser,
       FULL filtration (no threshold, every H1 bar has a finite death) -> relative persistence
       rel = (death - birth) / D.  A HOLE is a bar with rel > 0.4.  holes = number of such bars.
     b2 (cavities) is NOT part of the reader (declared consequence: it cannot see a sphere-like cavity; a torus
     is recognised through its H1 bars only). Dimension (geo_core.ball_dimension) and MDS are reported, not graded.
  2. Alternative reader (pre-declared, used ONLY if the first fails R1): bootstrap stability. A hole counts only
     if rel > 0.4 on the full data AND on >= 8 of 10 random 70 % site subsamples (d restricted, geodesic recomputed).
  3. Surfer runs (lam = 0) and the reinforcement rule with link weights (pre_core.LatticeWalk).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import numpy as np
import scipy.sparse.linalg as spl
import geo_core as G
from pre_core import LatticeWalk, modes_state, surfer_tick, sample_sites, hbar, windings, fit_tau

THR = 0.4
NLAND = 220


# ------------------------------------------------------------------ R1 reader
def _h1_rel(Gm, nland=NLAND):
    from ripser import ripser
    diam = float(np.max(Gm[np.isfinite(Gm)]))
    L = G.landmarks(Gm, nland)
    D1 = ripser(Gm[np.ix_(L, L)], distance_matrix=True, maxdim=1)["dgms"][1]
    if len(D1) == 0:
        return np.zeros(0), diam, D1
    assert np.all(np.isfinite(D1[:, 1])), "H1 bar with infinite death"
    pers = D1[:, 1] - D1[:, 0]
    o = np.argsort(pers)[::-1]
    return pers[o] / diam, diam, D1[o]


def read_shape(I, nland=NLAND, with_mds=True):
    nc, lab = G.components(I)
    sizes = np.bincount(lab)
    c = int(np.argmax(sizes))
    idx = np.where(lab == c)[0]
    d = G.mi_distance(I[np.ix_(idx, idx)])
    Gm, tau = G.geodesic(d)
    rel, diam, bars = _h1_rel(Gm, nland)
    dsm, dlg, h, _ = G.ball_dimension(Gm)
    out = dict(pieces=int(nc), size=int(len(idx)), sizes=sorted(sizes.tolist(), reverse=True)[:6], diam=diam,
               rel=rel[:5].tolist(), holes=int((rel > THR).sum()), D_small=float(dsm), D_large=float(dlg),
               bars=bars[:5].tolist())
    if with_mds:
        X, ev = G.mds(Gm, 3)
        out["X"] = X; out["ev"] = ev[:8]; out["sites"] = idx; out["d"] = d
    return out


def read_shape_boot(I, nboot=10, frac=0.7, need=8, seed=0):
    """pre-declared alternative: hole must survive >= need of nboot 70 % subsamples"""
    base = read_shape(I, with_mds=True)
    d = base["d"]
    n = len(d)
    rng = np.random.default_rng(seed)
    cnt = 0
    hb = []
    for b in range(nboot):
        keep = np.sort(rng.choice(n, int(frac * n), replace=False))
        Gm, _ = G.geodesic(d[np.ix_(keep, keep)])
        rel, _, _ = _h1_rel(Gm)
        k = int((rel > THR).sum())
        hb.append(k)
        cnt += k >= 1
    base["boot_holes"] = hb
    base["holes_boot"] = int(base["holes"] >= 1 and cnt >= need)
    return base


# ------------------------------------------------------------------ vacuum MI of a walk
def vacuum_mi(walk):
    C, gap, nocc = G.vacuum(walk.matrix())
    return G.mutual_information(C, walk.n * walk.n), gap


def low_spectrum(walk, k=12):
    """lowest |omega| quasi-energies of the one-tick U (near omega = 0) and near omega = pi"""
    U = walk.matrix().tocsc()
    D = U.shape[0]
    v0 = np.random.default_rng(3).normal(size=D) + 0j
    out = {}
    for name, s in (("zero", 1.0 + 0j), ("pi", -1.0 + 0j)):
        vals, vecs = spl.eigs(U, k=k, sigma=s * np.exp(-1j * 1e-3), which="LM", tol=1e-10, v0=v0)
        om = -np.angle(vals) if name == "zero" else np.angle(-vals)
        o = np.argsort(np.abs(om))
        out[name] = np.abs(om[o])
        if name == "zero":
            V = vecs[:, o]
            rho = (np.abs(V)**2).reshape(4, walk.n * walk.n, -1).sum(0)
            rho /= rho.sum(0, keepdims=True)
            out["ipr_frac"] = 1.0 / (walk.n * walk.n * (rho**2).sum(0))   # fraction of sites occupied
            out["rho_low"] = rho[:, 0].reshape(walk.n, walk.n)
    return out


# ------------------------------------------------------------------ surfer runs
def gini(x):
    x = np.sort(np.asarray(x, float)); k = len(x)
    return float((2 * np.arange(1, k + 1) - k - 1).dot(x) / (k * x.sum()))


def top10(x):
    x = np.sort(np.asarray(x, float))[::-1]
    return float(x[:max(1, len(x) // 10)].sum() / x.sum())


def interior(W, n):
    return np.concatenate([W[0][:-1, :].ravel(), W[1][:, :-1].ravel()])


def set_interior(vals, n):
    W0 = np.zeros((n, n)); W1 = np.zeros((n, n))
    m = (n - 1) * n
    W0[:-1, :] = vals[:m].reshape(n - 1, n)
    W1[:, :-1] = vals[m:].reshape(n, n - 1)
    return [W0, W1]


def run(n, theta, M, start, T, N=20000, cell=4, psi_seed=101, surf_seed=201, lam=0.0, omega0=None,
        every=10, reinf=None, w0=1.0, vort_every=None, snap=()):
    """box n x n (mirror edges). reinf = dict(eta, mu, wbar, wmin, K) or None. w0 = starting link weight."""
    if omega0 is None:
        omega0 = 1.8 * theta
    wk = LatticeWalk(n, theta, periodic=False, w=[np.full((n, n), w0)] * 2)
    psi, om = modes_state(wk, M, omega0=omega0, seed=psi_seed)
    rng = np.random.default_rng(surf_seed)
    rho = np.sum(np.abs(psi)**2, 0)
    xx, yy = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    prob = rho if start == "psi2" else np.ones((n, n))
    X = sample_sites(rng, N, prob)
    counts = [np.zeros((n, n)), np.zeros((n, n))] if reinf else None
    ts, H, V, Gi = [], [], [], []
    snaps = {}
    ncell = 1
    for t in range(T + 1):
        if t % every == 0:
            rho = np.sum(np.abs(psi)**2, 0)
            h, ncell = hbar(X, rho, cell)
            ts.append(t); H.append(h)
            if reinf:
                Gi.append(gini(interior(wk.W, n)))
        if vort_every and t % vort_every == 0:
            V.append(windings(psi))
        if t in snap:
            snaps[t] = [w.copy() for w in wk.W]
        if t == T:
            break
        psi, X = surfer_tick(wk, psi, X, rng, lam, counts)
        if reinf and (t + 1) % reinf["K"] == 0:
            J = interior(counts, n)
            w = interior(wk.W, n)
            Jm = J**reinf["mu"]
            w = (1 - reinf["eta"]) * w + reinf["eta"] * reinf["wbar"] * Jm / max(Jm.mean(), 1e-12)
            for _ in range(3):
                w = np.clip(w, reinf["wmin"], 1.0)
                w *= reinf["wbar"] / w.mean()
            w = np.clip(w, reinf["wmin"], 1.0)
            wk.set_w(set_interior(w, n))
            counts = [np.zeros((n, n)), np.zeros((n, n))]
    F = (ncell - 1) / (2 * N)
    H = np.array(H)
    return dict(t=np.array(ts), H=H, F=F, tau=fit_tau(ts, H, F), ratio=float(H[-1] / H[0]), H0=float(H[0]),
                V=np.array(V), Gini=np.array(Gi), walk=wk, psi=psi, snaps=snaps, omega=om)
