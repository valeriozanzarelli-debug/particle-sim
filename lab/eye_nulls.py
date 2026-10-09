"""eye_nulls.py CLASS - L1 exam nulls: 25 fresh instances of one class, seed sequence (20261009+32, class index).
Writes lab/out/eye_nulls_CLASS.json."""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import eye_core as E

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
THR2 = json.load(open(os.path.join(OUT, "eye_calib.json")))["THR2"]
if __name__ == "__main__":
    c = sys.argv[1]
    t0 = time.time()
    Ib = np.load(os.path.join(OUT, "eye_nullbase.npy"))
    rng = np.random.default_rng([20261009 + 32, E.CLASSES.index(c)])
    L = []
    for k in range(25):
        r = E.strip(E.eye(E.make_null(c, Ib, rng), thr2=THR2))
        L.append(r)
        print(k, r["loops"], r["cavities"], round(r["rel2"][0], 3), r["cones"], r["cones_v2"], f"{time.time()-t0:.0f}s",
              flush=True)
    summ = dict(n=25, fp_cavities=float(np.mean([r["cavities"] >= 1 for r in L])),
                fp_cones_v1=float(np.mean([r["cones"] >= 1 for r in L])),
                fp_cones_v2=float(np.mean([r["cones_v2"] >= 1 for r in L])),
                fp_loops=float(np.mean([r["loops"] >= 1 for r in L])),
                max_rel2=float(max(r["rel2"][0] for r in L)), refused_v2=int(sum(r["refused_v2"] for r in L)))
    print(summ)
    json.dump(dict(cls=c, THR2=THR2, summary=summ, nulls=L, time_s=time.time() - t0),
              open(os.path.join(OUT, f"eye_nulls_{c}.json"), "w"), indent=1, default=float)
