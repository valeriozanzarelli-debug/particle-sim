"""eye_calib.py - L1 calibration: 4 null classes x 10 (seed 20261009+31) -> THR2 (pre-registered rule).
Writes lab/out/eye_calib.json."""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import eye_core as E

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
if __name__ == "__main__":
    t0 = time.time()
    p = os.path.join(OUT, "eye_nullbase.npy")
    Ib = np.load(p) if os.path.exists(p) else E.null_base()
    rng = np.random.default_rng(20261009 + 31)
    res = {c: [] for c in E.CLASSES}
    for k in range(10):
        for c in E.CLASSES:
            r = E.strip(E.eye(E.make_null(c, Ib, rng)))
            res[c].append(r)
        print(k, {c: (round(res[c][-1]["rel2"][0], 3), res[c][-1]["cones"]) for c in E.CLASSES},
              f"{time.time()-t0:.0f}s", flush=True)
    mx = max(r["rel2"][0] if r["rel2"] else 0 for L in res.values() for r in L)
    thr2 = max(0.2, 1.25 * mx)
    out = dict(max_rel2=mx, THR2=thr2, per_class_max_rel2={c: max(r["rel2"][0] if r["rel2"] else 0 for r in L)
                                                            for c, L in res.items()},
               cones_seen={c: [r["cones"] for r in L] for c, L in res.items()}, nulls=res, time_s=time.time() - t0)
    print("THR2", thr2, out["per_class_max_rel2"], out["cones_seen"])
    json.dump(out, open(os.path.join(OUT, "eye_calib.json"), "w"), indent=1, default=float)
