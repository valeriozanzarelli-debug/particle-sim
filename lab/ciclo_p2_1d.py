"""ciclo_p2_1d.py - phase 2 (1-D) production: python lab/ciclo_p2_1d.py <free|coulomb> <2|3|fib>   or   anchor
Writes lab/out/ciclo_p2_1d_<mode>_<scheme>.json (one run per call, each well under 10 min)."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ciclo_p2_core as c

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
TH, B0, WB = 0.05, 32, 8
G_C, S_C = 0.015, 10.0                       # soft Coulomb strength (dimensionless) and core (lambda_C)
SCHEMES = {"2": (2, 2048, None), "3": (3, 2187, None), "fib": ("fib", 2584, 16)}
KOBS = 2 * np.pi * 4 / 2048                  # base wave number of the observables (same physical k for all)
WMAX = 8


def soft_V(x):
    return -G_C / np.sqrt(x * x + S_C * S_C)


def base_potential(n0):
    """block mean of the fine per-tick potential V(x) theta over the B0 fine cells of each base cell"""
    xf = (np.arange(n0 * B0) - n0 * B0 / 2 + 0.5) * TH
    vf = soft_V(xf) * TH
    return vf.reshape(n0, B0).mean(1), vf


def run(mode, sch):
    kind, n0, nsub = SCHEMES[sch]
    g, berr = c.make_base(TH, B0, WB)
    H0 = c.circulant(g, n0)
    kk = 2 * np.pi * np.fft.fftfreq(n0)
    nobs = int(np.sum(np.abs(kk) <= KOBS + 1e-12))
    V0 = base_potential(n0)[0] if mode == "coulomb" else None
    t0 = time.time()
    out = c.run_chain(kind, H0, nlev=12, eps_level=1e-3, wmax=WMAX, nobs=nobs, V0=V0, nb=3, cmax=2, eps_c=2e-2, nsub=nsub,
                      log=lambda s: print(s, flush=True))
    out.update(mode=mode, base_err=berr, base_stencil=g.tolist(), theta=TH, B0=B0, nobs=nobs, wall=time.time() - t0,
               params=dict(FRAC_CHECK=c.FRAC_CHECK, FRAC_SLOW=c.FRAC_SLOW, MIN_CHECK=c.MIN_CHECK, W_FAST=c.W_FAST,
                           REL_FLOOR=c.REL_FLOOR, NMIN_CELLS=c.NMIN_CELLS, WMAX=WMAX))
    json.dump(out, open(os.path.join(OUT, f"ciclo_p2_1d_{mode}_{sch}.json"), "w"), indent=1, default=float)
    print("done", mode, sch, round(time.time() - t0, 1), "s")


def anchor():
    """fine zig-zag bound states (N = 32 x 2048 fine cells) vs the base level (2048 cells)"""
    n0 = 2048
    g, berr = c.make_base(TH, B0, WB)
    V0, vf = base_potential(n0)
    eb = np.sort(np.linalg.eigvalsh(c.circulant(g, n0).toarray() + np.diag(V0)))[:3]
    N = n0 * B0
    U = c.fine_U(N, TH, vf)
    sig = np.exp(-1j * (TH + 1.2 * eb[0]))
    lam = spla.eigs(U, k=12, sigma=sig, return_eigenvectors=False, tol=1e-12, maxiter=50000)
    om = np.sort(-np.angle(lam)) - TH
    ob = om[om < 0][:3]
    rel = (np.abs(eb - ob) / np.abs(ob)).tolist()
    res = dict(base_bound=eb.tolist(), fine_bound=ob.tolist(), rel=rel, base_err=berr)
    json.dump(res, open(os.path.join(OUT, "ciclo_p2_1d_anchor.json"), "w"), indent=1)
    print(res)


if __name__ == "__main__":
    if sys.argv[1] == "anchor":
        anchor()
    else:
        run(sys.argv[1], sys.argv[2])
