"""zstar_run.py - production walk runs of ZSTAR (after lab/out/zstar_prereg.json was committed).
The frozen radial rule with the variants at the origin; Z_cr(a) by the giro-1 tracker (supercrit_run_c1).
Checkpoint: lab/out/zstar_runs.json. Each call stops starting new work after ~400 s.

    python zstar_run.py VARIANT [a ...]       # VARIANT in zstar_prereg.VARIANTS or 'split' / 'order' / 'refl1'
    python zstar_run.py K+2 [a ...]           # contact-state crossing scan (d3/2)
"""
import json, os, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supercrit_core import SWalk, coins, level, gap_levels
from supercrit_run_c1 import track_to_critical
from zstar_core import coin_variant
from zz_core import ALPHA

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "zstar_runs.json")
PRE = json.load(open(os.path.join(HERE, "out", "zstar_prereg.json")))


class ZWalk(SWalk):
    """SWalk with sites r_i = (i + s) a, complex refl0, optional coin variant"""

    def __init__(self, kappa, a, rmax, Za, s=0.5, refl0=-1.0, refl1=-1.0, variant="frozen", mass=1.0):
        self.kappa, self.a, self.m = kappa, a, mass
        self.Za = Za
        self.N = int(np.ceil(rmax / a))
        self.r = (np.arange(self.N) + s) * a
        self.theta = mass * a
        self.refl0, self.refl1 = refl0, refl1
        self.V = -Za / self.r
        self.M = coin_variant(kappa, a, self.r, self.V, mass, variant)
        self.Mh = coins(kappa, a, self.r, self.V, mass, frac=0.5)   # only used for plotting states


def spec(name):
    from zstar_prereg import VARIANTS
    if name in VARIANTS:
        kappa, s, refl0, _ = VARIANTS[name]
        return dict(kappa=kappa, s=s, refl0=refl0, refl1=-1.0, variant="frozen")
    if name == "split":
        return dict(kappa=-1, s=0.5, refl0=-1.0, refl1=-1.0, variant="split")
    if name == "refl1":
        return dict(kappa=-1, s=0.5, refl0=-1.0, refl1=1.0, variant="frozen")
    raise KeyError(name)


def z0_for(name, sp_):
    if name in PRE["predictions"]:
        roots = PRE["predictions"][name]["Z_roots"]
        if sp_["kappa"] == -2:
            return 180.0
        if roots:
            return min(120.0, roots[-1] - 20.0) if abs(sp_["kappa"]) == 1 else 120.0
    return 120.0


def run_dive(name, alist, res, t0):
    sp_ = spec(name)
    rmax = PRE["configs"]["rmax"]
    for a in alist:
        key = f"{name}|{a:.6g}"
        if key in res or time.time() - t0 > 400:
            continue
        mk = lambda Z: ZWalk(sp_["kappa"], a, rmax, Z * ALPHA, sp_["s"], sp_["refl0"], sp_["refl1"], sp_["variant"])
        Z0 = z0_for(name, sp_)
        print(name, "a =", a, "Z0 =", Z0, flush=True)
        t1 = time.time()
        Zcr, track = track_to_critical(mk, Z0, log=lambda *x: None, Zmax=300.0)
        res[key] = dict(variant=name, a=a, Zcr=Zcr, Z0=Z0, n_track=len(track), track_tail=track[-4:],
                        seconds=time.time() - t1)
        print(f"  {name} a={a}: Z_cr = {Zcr:.5f} ({time.time() - t1:.0f}s)", flush=True)
        json.dump(res, open(OUT, "w"), indent=1)


def contact_levels(w, sigma=0.0, k=8, rin=50):
    """gap eigenstates of W near omega = sigma; returns (omega, weight inside r < rin a)"""
    W = w.matrix()
    vals, vecs = spl.eigs(W, k=k, sigma=np.exp(-1j * sigma * w.a), which="LM", tol=1e-12)
    out = []
    for j in range(len(vals)):
        om = -np.angle(vals[j]) / w.a
        x = vecs[:, j].reshape(w.N, 2).T
        rho = np.abs(x[0])**2 + np.abs(x[1])**2
        out.append((om, rho[:rin].sum() / rho.sum()))
    return out


def run_contact(alist, res, t0):
    """K+2: find Z where the localized contact state crosses omega = 0"""
    rmax = 10.0
    for a in alist:
        key = f"K+2|{a:.6g}"
        if key in res or time.time() - t0 > 400:
            continue
        def om_contact(Z):
            w = ZWalk(2, a, rmax, Z * ALPHA)
            L = [x for x in contact_levels(w) if x[1] > 0.5]
            return L
        scan = []
        Zs = np.linspace(40.0, 44.0, 41)
        for Z in Zs:
            L = om_contact(Z)
            scan.append((float(Z), [float(o) for o, f in L]))
        # Z where a contact state has omega = 0: sign change between consecutive scan points (+ bisection)
        Zx = None
        for i in range(len(scan) - 1):
            if scan[i][1] and scan[i + 1][1]:
                o1 = min(scan[i][1], key=abs); o2 = min(scan[i + 1][1], key=abs)
                if np.sign(o1) != np.sign(o2):
                    lo, hi, flo = scan[i][0], scan[i + 1][0], o1
                    for _ in range(30):
                        mid = 0.5 * (lo + hi)
                        L = om_contact(mid)
                        if not L:
                            break
                        om = min([o for o, f in L], key=abs)
                        if np.sign(om) == np.sign(flo):
                            lo, flo = mid, om
                        else:
                            hi = mid
                    Zx = 0.5 * (lo + hi)
                    break
        res[key] = dict(variant="K+2", a=a, Zx=Zx, scan=scan)
        print(f"  K+2 a={a}: Z_x = {Zx}", flush=True)
        json.dump(res, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    name = sys.argv[1]
    alist = [float(x) for x in sys.argv[2:]] or PRE["configs"]["a"]
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    t0 = time.time()
    if name == "K+2":
        run_contact(alist, res, t0)
    elif name == "order":
        # W = M S versus W' = S M at a = 1e-3, Z = 132 (identical spectrum: similar matrices)
        a, Z = 1e-3, 132.0
        w = ZWalk(-1, a, 40.0, Z * ALPHA)
        W = w.matrix()
        M = sp.block_diag([sp.csc_matrix(w.M[i]) for i in range(w.N)], format="csc")
        Minv = sp.block_diag([sp.csc_matrix(np.conj(w.M[i].T)) for i in range(w.N)], format="csc")
        W2 = ((Minv @ W) @ M).tocsc()                     # S M  with  S = M^-1 W
        o1, *_ = level(w, -0.5)
        lam0 = np.exp(-1j * o1 * a)
        vals1 = spl.eigs(W.tocsc(), k=4, sigma=lam0, which="LM", return_eigenvectors=False)
        vals2 = spl.eigs(W2, k=4, sigma=lam0, which="LM", return_eigenvectors=False)
        d = np.max(np.abs(np.sort_complex(vals1) - np.sort_complex(vals2)))
        res["order"] = dict(a=a, Z=Z, omega=o1, max_dlambda=float(d))
        print("order: max |d lambda| =", d)
        json.dump(res, open(OUT, "w"), indent=1)
    else:
        run_dive(name, alist, res, t0)
