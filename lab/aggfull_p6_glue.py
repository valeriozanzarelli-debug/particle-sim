"""aggfull_p6_glue.py - phase 6 glue test: the zoom chain of NUCLEO_1S phase 8 (lab/nucleo1s_f8_core.py, imported READ
ONLY, not modified) glued to the coarse level that now carries the field-motion coupling LEARNED in the large-component
language (lab/aggfull_p6.py -> lab/out/aggfull_p6_learn.pkl).

Outer operator (s wave, u = r psi, cells of the T level):
    H = K + g0 V + (a+b+c) V P2 - (2a+c) grad V.grad,    P2 = 2 m_kin K
from the learned small-k coefficients Gamma_LC = g0 + a k'^2 + b k^2 + c k'.k, rewritten with q = k' - k:
g0 + (a+b+c) k^2 + (2a+c) q.k + a q^2. The a q^2 part is the contact term (-a) lap V = Darwin's delta at the nucleus,
which lives inside the zoom region (r < R_f) and is NOT added outside (the chain carries it). (A first version used the
lattice products P2 V and V P2: the odd-reflected V u jumps at r = 0 and the lattice delta leaks 1/n^2 tails to every
row, wrong sign far out: rejected in development, declared.) Debt: the operator basis (V p^2, grad V.grad, lap V) and the
8th-order radial derivative are mine; the three coefficients are learned. The imaginary part of a, b (i theta/4 (k'^2 - k^2): where in the fine tick the
field acts; a unitary relabelling of time, zero on eigenvalues) is dropped: the chain's stationary solution is read at the
palindromic point (real), the same convention. K = the T level's own s-wave kinetic matrix (nucleo1s_f8_core.CoarseT).

Two-Z difference as in nucleo1s phase 8: (Z = 1, B = 512) and (Z = 2, B = 256), same geometry in a0/Z.
Usage: timeout 590 python lab/aggfull_p6_glue.py run    (repeat until 'all done'; results lab/out/aggfull_p6_glue.json)
"""
import os, sys, json, pickle, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nucleo1s_f8_core as n8                                     # read only
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
DEV = os.environ.get("AGGFULL_DEV") == "1"
PFX = "dev_" if DEV else ""


def lc_coeffs(B, src=None):
    import aggfull_p6 as p6
    st = pickle.load(open(src or os.path.join(OUT, PFX + "aggfull_p6_learn.pkl"), "rb"))
    s = [x for x in st if x["B"] == B][0]
    ph = p6.lc_physics(s["C"], s["offs"], s["h"])
    a, b, c, g0 = (complex(ph[x]).real for x in ("a", "b", "c", "g0"))
    return dict(g0=g0, vp2=a + b + c, gg=2 * a + c, contact=-a, a=a, b=b, c=c)


def deriv_matrix(N, h, order=8):
    """central finite-difference d/dr on cells r_j = (j + 1/2) h, odd reflection u_{-j-1} = -u_j, u = 0 beyond N"""
    m = order // 2
    offs = np.arange(-m, m + 1)
    A = np.vander(offs, increasing=True).T.astype(float)
    rhs = np.zeros(len(offs)); rhs[1] = 1.0
    w = np.linalg.solve(A, rhs)
    D = np.zeros((N, N))
    for j in range(N):
        for o, wo in zip(offs, w):
            i = j + o
            if i >= N or wo == 0:
                continue
            if i < 0:
                D[j, -i - 1] -= wo
            else:
                D[j, i] += wo
    return D / h


def lc_outer(c, N, Za, coef, K):
    """the learned field-motion term outside the nucleus, delta-free:  Gamma = g0 + a k'^2 + b k^2 + c k'.k
    = g0 + (a+b+c) k^2 + (2a+c) q.k + a q^2   ->   (g0-1) V + (a+b+c) V p^2 - (2a+c) grad V . grad  [+ (-a) lap V: contact,
    inside the chain, not added]. s wave, u = r psi: grad V . grad psi -> V'(r) (u' - u/r)."""
    r = (np.arange(N) + 0.5) * c.h
    V = c.V(N, Za)
    dV = Za / r ** 2
    P2 = 2 * c.m_kin * K
    D = deriv_matrix(N, c.h)
    return ((coef["g0"] - 1) * np.diag(V) + coef["vp2"] * (V[:, None] * P2)
            - coef["gg"] * (dV[:, None] * (D - np.diag(1 / r))))


class GlueLC(n8.Glue8):
    def __init__(self, c, Z, J0, rmax, coef, zoom_levels=12, npl=16):
        super().__init__(c, Z, J0, rmax, zoom_levels, npl)
        if coef is not None:
            self.Hm = self.Hm + lc_outer(c, self.N, self.Za, coef, self.K)


def run(cfg, resfile):
    res = json.load(open(resfile)) if os.path.exists(resfile) else {}
    T0 = time.time()
    coarse, coefs = {}, {}
    for mode in cfg["modes"]:
        for J0 in cfg["J0"]:
            for pz in cfg["pair"]:
                Z, B = pz["Z"], pz["B"]
                for n in cfg["states"]:
                    key = f"{mode}_J{J0}_Z{Z}_n{n}"
                    if key in res:
                        continue
                    if time.time() - T0 > cfg.get("budget", 480):
                        print("time budget: run again", flush=True)
                        return res, False
                    if B not in coarse:
                        coarse[B] = n8.CoarseT(B, cfg.get("band", "learned"), ndir=cfg["ndir"])
                        coefs[B] = lc_coeffs(B)
                    t = time.time()
                    g = GlueLC(coarse[B], Z, J0, cfg["rmax_a0_over_Z"] * n8.A0 / Z, coefs[B] if mode == "lc" else None,
                               cfg["zoom_levels"], cfg["n_per_level"])
                    r = g.solve(n)
                    if r is None:
                        r = dict(B=None)
                    r.update(sec=time.time() - t)
                    res[key] = r
                    json.dump(res, open(resfile, "w"), indent=1, default=float)
                    print(key, f"{r.get('rel', float('nan')):+.4e}", f"{r['sec']:.1f}s", flush=True)
    return res, True


def errors(res, cfg):
    out = {}
    for mode in cfg["modes"]:
        for J0 in cfg["J0"]:
            for n in cfg["states"]:
                b1, b2 = res.get(f"{mode}_J{J0}_Z1_n{n}", {}).get("B"), res.get(f"{mode}_J{J0}_Z2_n{n}", {}).get("B")
                if b1 is None or b2 is None:
                    continue
                out[f"{mode}_J{J0}_{n}s"] = n8.delta_pair(b1, b2) / n8.delta_dirac(n) - 1
    return out


if __name__ == "__main__":
    if sys.argv[1] == "run":
        P = json.load(open(os.path.join(OUT, "aggfull_p6_prereg.json")))["config"]["glue"]
        cfg = dict(modes=P["modes"], J0=P["J0"], pair=[dict(Z=1, B=512), dict(Z=2, B=256)], states=P["states"],
                   rmax_a0_over_Z=P["rmax_a0_over_Z"], zoom_levels=P["zoom_levels"], n_per_level=P["n_per_level"], ndir=P["ndir"])
        res, done = run(cfg, os.path.join(OUT, "aggfull_p6_glue.json"))
        if done:
            json.dump(errors(res, cfg), open(os.path.join(OUT, "aggfull_p6_glue_err.json"), "w"), indent=1)
            print("all done")
    if sys.argv[1] == "dev":
        cfg = dict(modes=["base", "lc"], J0=[int(x) for x in sys.argv[2].split(",")], pair=[dict(Z=1, B=512), dict(Z=2, B=256)],
                   states=[1, 2], rmax_a0_over_Z=90, zoom_levels=12, n_per_level=16, ndir=120)
        res, done = run(cfg, os.path.join(OUT, "dev_aggfull_p6_glue.json"))
        for k, v in errors(res, cfg).items():
            print(k, f"{v:+.4f}")
