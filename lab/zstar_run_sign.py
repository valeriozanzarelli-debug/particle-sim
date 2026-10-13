"""zstar_run_sign.py - AMENDMENT A3 (lab/zstar_amend.py): Z_cr(a) of the walk WITHOUT eigen-tracking.
At omega = -1 (lambda = e^{+i a}) the walk's own regular solution from r = 0 (outward recursion of W chi = lambda chi,
exact, same rule) is followed to r_end = 20; it is e^{i phi} x real with phi ~ +-pi/2, and Im g(r_end) changes sign
exactly when a level of W crosses omega = -1. Z is scanned upward with step h from Z_lo, first sign change, then
bisection to 1e-7. Validated on every variant where the tracker worked. Checkpoint lab/out/zstar_runs_sign.json.
    python zstar_run_sign.py VARIANT Z_lo Z_hi h [a ...]"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zstar_run import ZWalk, spec
from zstar_core import _fwd
from zz_core import ALPHA

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "zstar_runs_sign.json")


OMEGA = float(os.environ.get("ZSTAR_OMEGA", "-1"))     # crossing energy (K+2 contact state: also omega = 0)


def h_of(Z, a, sp_, rend=20.0):
    w = ZWalk(sp_["kappa"], a, rend, Z * ALPHA, sp_["s"], sp_["refl0"], sp_["refl1"], sp_["variant"])
    cp, cm, _ = _fwd(w.M, complex(sp_["refl0"]), complex(np.exp(-1j * OMEGA * a)))
    g = cp[-1] + cm[-1]
    return float(g.imag / abs(g)), float(abs(np.cos(np.angle(g))))


def first_crossing(a, sp_, Zlo, Zhi, h):
    Z = Zlo
    f0, _ = h_of(Z, a, sp_)
    re_max = 0.0
    while Z < Zhi:
        Zn = min(Z + h, Zhi)
        f1, re = h_of(Zn, a, sp_)
        re_max = max(re_max, re)
        if np.sign(f1) != np.sign(f0):
            lo, hi, flo = Z, Zn, f0
            while hi - lo > 1e-7:
                mid = 0.5 * (lo + hi)
                fm, _ = h_of(mid, a, sp_)
                if np.sign(fm) == np.sign(flo):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return 0.5 * (lo + hi), re_max
        Z, f0 = Zn, f1
    return None, re_max


if __name__ == "__main__":
    name, Zlo, Zhi, h = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    alist = [float(x) for x in sys.argv[5:]] or [3e-3, 1e-3, 3e-4, 1e-4]
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    sp_ = spec(name)
    t0 = time.time()
    for a in alist:
        key = f"{name}|{a:.6g}" + ("" if OMEGA == -1 else f"|omega={OMEGA:g}")
        if key in res or time.time() - t0 > 450:
            continue
        t1 = time.time()
        Zc, re = first_crossing(a, sp_, Zlo, Zhi, h)
        res[key] = dict(variant=name, a=a, Zcr=Zc, Zlo=Zlo, Zhi=Zhi, h=h, max_real_part=re, seconds=time.time() - t1)
        print(f"{name} a={a}: Z_cr = {Zc}  (|cos arg g| max {re:.1e}, {time.time() - t1:.0f}s)", flush=True)
        json.dump(res, open(OUT, "w"), indent=1)
