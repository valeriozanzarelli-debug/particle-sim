"""supercrit_run_c1.py - exam C1: point nucleus, the lowest s level of the rule followed in Z until it touches
-mc^2 (Z_cr(a)). Tracking: shift-invert eigenvector of W from the previous omega, steps in Z, Newton steps on
omega(Z) + 1 = 0 with the local slope, final value by quadratic extrapolation of the last 3 points (omega > -0.999).
Checkpointed in lab/out/supercrit_c1.json. Usage: python supercrit_run_c1.py [a ...]"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supercrit_core import SWalk, V_point, level, gap_levels, ALPHA

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "supercrit_c1.json")


DZMAX = 0.5          # amendment A1 (after the first pass): max step 0.5 (4.0 skipped over a dive at a = 1e-4, 3e-5)


def track_to_critical(make_walk, Z0, log=print, om_stop=-0.9995, Zmax=260.0, sigma0=0.3):
    """follow the lowest localized gap level from Z0 upward until omega -> -1. returns Z_cr and the track."""
    w = make_walk(Z0)
    L = [x for x in gap_levels(w, sigma0, k=12) if x[2] > 0.99]
    om, nodes = L[0][0], L[0][1]
    track = [(Z0, om, nodes)]
    Z = Z0
    dZ = 1.0
    fresh = True
    slope = None
    while True:
        # predict next Z (Newton step only right after a successful step; after a failure keep the halved dZ)
        if len(track) >= 2:
            (Za_, oa), (Zb_, ob) = track[-2][:2], track[-1][:2]
            slope = (ob - oa) / (Zb_ - Za_)
            if slope < 0 and fresh:
                dZ_newton = (-1 - ob) / slope
                dZ = min(max(dZ_newton * 0.5, 1e-5), DZMAX)
            om_pred = ob + slope * dZ
        else:
            om_pred = om - 0.02
        fresh = False
        Zn = Z + dZ
        if Zn > Zmax:
            raise RuntimeError("no critical charge below Zmax")
        w = make_walk(Zn)
        omn, rho, resid, smooth = level(w, max(om_pred, -0.9999))
        inside = rho[: w.N // 2].sum() / rho.sum()
        jump = len(track) >= 2 and abs(omn - om_pred) > 0.5 * abs(om_pred - track[-1][1]) + 0.01
        if omn > track[-1][1] + 1e-9 or omn < -1 or inside < 0.99 or jump:
            # lost the level (overshoot past -1, or jumped to another state): halve the step
            dZ *= 0.4
            if dZ < 1e-6:
                break
            continue
        Z = Zn
        fresh = True
        track.append((Z, omn, inside))
        log(f"   Z={Z:.5f} omega={omn:.6f} inside={inside:.5f} resid={resid:.1e}")
        if omn < om_stop:
            break
    Zs = np.array([t[0] for t in track[-3:]])
    Os = np.array([t[1] for t in track[-3:]])
    c = np.polyfit(Os, Zs, 2)                    # Z as a function of omega, evaluated at -1
    return float(np.polyval(c, -1.0)), track


if __name__ == "__main__":
    from supercrit_prereg import P
    alist = [float(x) for x in sys.argv[1:]] or P["configs"]["C1"]["a"]
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    t0 = time.time()
    for a in alist:
        key = f"{a:.6g}"
        if key in res:
            continue
        if time.time() - t0 > 400:
            break
        rmax = P["configs"]["C1"]["rmax"]
        mk = lambda Z: SWalk(-1, a, rmax, lambda r: V_point(r, Z * ALPHA))
        print("a =", a, flush=True)
        Zcr, track = track_to_critical(mk, 120.0)
        res[key] = dict(a=a, Zcr=Zcr, zeta_cr=Zcr * ALPHA, track=track, N=int(np.ceil(rmax / a)))
        print(f"a={a}: Z_cr = {Zcr:.5f}  ({time.time() - t0:.0f}s)", flush=True)
        json.dump(res, open(OUT, "w"), indent=1)
