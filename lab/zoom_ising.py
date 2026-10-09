"""ZOOM test b: 2-D Ising model, block-spin renormalisation with the SAME agglomerator (zoom_core).

Exact answers (Onsager 1944 and the 2-D Ising CFT): T_c = 2 / ln(1 + sqrt 2) = 2.269185 (K_c = 0.440687),
thermal exponent y_t = 1/nu = 1 (nu = 1), magnetic exponent y_h = 15/8 = 1.875 (eta = 1/4).

Pipeline:
  fine level   Swendsen-Wang cluster Monte Carlo of the nearest-neighbour model, L = 256, periodic
               (brute force; the only system-specific code here).
  block        zoom_core.block(k=2): super-ball charge = total spin, projected back to +-1 by majority
               (random tie-break) -> levels 128, 64, 32, 16.  Same code at every level.
  learn        zoom_core.LocalModel.learn_stat: maximum pseudo-likelihood of the generic cluster model
               (pair function phi(r) out to 6.5 spacings, field, 3- and 4-ball compact clusters).
               Same hyper-parameters at every level.  Couplings at T_c must stay the same from level to
               level (fixed point = the "fractal"); at T != T_c they must flow away.
  exponents    Swendsen MCRG: linearised map dK^(n)/dK^(n-1) = A^-1 B from correlations of the model's
               own operators (pair shells r = 1, sqrt2, 2 and the 4-ball cluster; odd: field and
               3-ball cluster); eigenvalues lambda = 2^y.  Errors by jackknife (20 blocks).
  check        heat-bath simulation of the LEARNED level-n rule started from blocked fine
               configurations; local statistics must not drift (spot-check vs the fine level).
This file: the Monte Carlo (writes lab/out/zoom_ising_mc.npz).  Learning, MCRG, checks and figures:
lab/zoom_ising_rg.py.
Run: timeout 590 python lab/zoom_ising.py && timeout 590 python lab/zoom_ising_rg.py
"""
import os, sys, time, json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zoom_core as zc

T0 = time.time()
OUT = os.path.join(HERE, "out")
TC = 2 / np.log(1 + np.sqrt(2))
KC = 1 / TC
L0 = 256
NLEV = 4
RUNS = [  # (name, T/Tc, L, sweeps, store_every)
    ("Tc", 1.00, L0, 3000, 10),
    ("below", 0.96, L0, 1000, 5),
    ("above", 1.04, L0, 1000, 5),
    ("Tc_small", 1.00, 64, 1500, 5),
]
NTHERM = 150


def log(*a):
    print(f"[{time.time() - T0:6.1f}s]", *a, flush=True)


def lattice_state(L, s):
    x, y = np.meshgrid(np.arange(L) + 0.5, np.arange(L) + 0.5, indexing="ij")
    return zc.State(np.c_[x.ravel(), y.ravel()], s.astype(float), 1.0, box=[L, L], alphabet=(-1, 1))


def sw_sweep(s, p, rng, idx, rr, dn):
    N = s.size
    b1 = (s == s[rr]) & (rng.random(N) < p)
    b2 = (s == s[dn]) & (rng.random(N) < p)
    i = np.r_[idx[b1], idx[b2]]; j = np.r_[rr[b1], dn[b2]]
    G = sp.coo_matrix((np.ones(i.size, np.int8), (i, j)), shape=(N, N)).tocsr()
    nc, lab = connected_components(G, directed=False)
    flip = rng.random(nc) < 0.5
    s[flip[lab]] *= -1


class Ops:
    """The model's own operators at one level (pair shells u^2 = 1, 2, 4; 4- and 3-ball clusters; field)."""
    NAMES_EVEN = ["pair r=1", "pair r=sqrt2", "pair r=2", "4-ball"]
    NAMES_ODD = ["field", "3-ball"]

    def __init__(self, st):
        g = zc.Geometry(st, rcut=2.01, dr=zc.DR)
        u2 = np.round(g.u ** 2, 3)
        self.sh = [(g.i[u2 == v], g.j[u2 == v]) for v in (1.0, 2.0, 4.0)]
        self.c3, self.c4 = g.cl3, g.cl4

    def __call__(self, q):
        e = [np.dot(q[i], q[j]) for i, j in self.sh] + [np.prod(q[self.c4], 1).sum()]
        o = [q.sum(), np.prod(q[self.c3], 1).sum()]
        return np.array(e + o)


def run(args):
    name, tt, L, nsw, every, seed = args
    rng = np.random.default_rng(seed)
    T = tt * TC
    p = 1 - np.exp(-2 / T)
    idx = np.arange(L * L)
    I = idx.reshape(L, L)
    rr = np.roll(I, -1, 1).ravel(); dn = np.roll(I, -1, 0).ravel()
    s = rng.choice(np.array([-1, 1], np.int8), L * L) if tt > 1 else np.ones(L * L, np.int8)
    nlev = int(np.log2(L // 16)) + 1 if L > 64 else 2
    st = lattice_state(L, s)
    blks, ops = [], [Ops(st)]
    c = st
    for n in range(nlev - 1):
        c, B = zc.block(c, 2, rng)
        blks.append(B); ops.append(Ops(c))
    for _ in range(NTHERM):
        sw_sweep(s, p, rng, idx, rr, dn)
    rec, store = [], [[] for _ in range(nlev)]
    t0 = time.time()
    for t in range(nsw):
        sw_sweep(s, p, rng, idx, rr, dn)
        st = lattice_state(L, s)
        levs = [st]
        for B in blks:
            levs.append(zc.block(levs[-1], 2, rng, blocking=B)[0])
        rec.append(np.array([o(v.q[:, 0]) for o, v in zip(ops, levs)]))
        if (t + 1) % every == 0:
            for n, v in enumerate(levs):
                store[n].append(v.q[:, 0].astype(np.int8))
    return dict(name=name, tt=tt, L=L, ops=np.array(rec), store=[np.array(x) for x in store],
                time=time.time() - t0, sweeps=nsw)


def main():
    log(f"T_c = {TC:.6f} (K_c = {KC:.6f});  runs: " + ", ".join(f"{r[0]} T/Tc={r[1]} L={r[2]} sweeps={r[3]}" for r in RUNS))
    with Pool(3) as pool:
        results = pool.map(run, [r + (1000 + i,) for i, r in enumerate(RUNS)])
    save = {}
    for r in results:
        log(f"MC {r['name']}: L={r['L']}, {r['sweeps']} SW sweeps in {r['time']:.0f}s ({1000 * r['time'] / r['sweeps']:.1f} ms/sweep incl. blocking+operators)")
        nm = r["name"]
        save[nm + "_ops"] = r["ops"]
        save[nm + "_meta"] = np.array([r["tt"], r["L"], r["time"], r["sweeps"]])
        for n, st in enumerate(r["store"]):
            save[f"{nm}_store{n}"] = st if (n > 0 or r["L"] <= 64) else st[-3:]
    # cost of the bare cluster update (no blocking), for the efficiency table
    rng = np.random.default_rng(5)
    L = L0; idx = np.arange(L * L); I = idx.reshape(L, L)
    rr = np.roll(I, -1, 1).ravel(); dn = np.roll(I, -1, 0).ravel()
    s = rng.choice(np.array([-1, 1], np.int8), L * L)
    for _ in range(20):
        sw_sweep(s, 1 - np.exp(-2 / TC), rng, idx, rr, dn)
    t0 = time.time()
    for _ in range(50):
        sw_sweep(s, 1 - np.exp(-2 / TC), rng, idx, rr, dn)
    save["sw_ms"] = np.array([1000 * (time.time() - t0) / 50])
    log(f"bare Swendsen-Wang sweep L={L}: {save['sw_ms'][0]:.1f} ms")
    np.savez_compressed(os.path.join(OUT, "zoom_ising_mc.npz"), **save)
    log("saved lab/out/zoom_ising_mc.npz; next: python lab/zoom_ising_rg.py")
    log("done")


if __name__ == "__main__":
    main()
