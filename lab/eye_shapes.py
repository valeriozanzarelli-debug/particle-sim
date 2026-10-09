"""eye_shapes.py - L1 exams on drawn shapes: T2/Z3 (L=27), T4/Z2 (n=6), flat T2 (L=15), S2 (icosphere 3).
Grading of cone positions uses the drawn graph (hop distance on the orbit graph) ONLY after the eye has read I.
Writes lab/out/eye_shapes.json and eye_shapes.npz (for the figure)."""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.sparse.csgraph import shortest_path
import eye_core as E

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
THR2 = json.load(open(os.path.join(OUT, "eye_calib.json")))["THR2"]


def grade_cones(sites, rho, Hop, fixed, lo, hi):
    """each cone within 1 hop of a fixed orbit, distinct fixed orbits, rho in [lo, hi]"""
    if len(fixed) == 0:
        return dict(n=len(sites), ok=len(sites) == 0)
    dist = Hop[np.ix_(sites, fixed)] if len(sites) else np.zeros((0, len(fixed)))
    near = [int(np.argmin(d)) if d.min() <= 1 else -1 for d in dist]
    ok = (len(sites) == len(fixed) and all(k >= 0 for k in near) and len(set(near)) == len(fixed)
          and all(lo <= r <= hi for r in rho))
    return dict(n=len(sites), nearest_fixed=near, min_hops=[float(d.min()) for d in dist], rho=rho, ok=bool(ok))


if __name__ == "__main__":
    t0 = time.time()
    res, npz = {}, {}
    shapes = {}
    Hq, sz, fx, orb = E.t2z3(27); shapes["T2/Z3"] = (Hq, fx, (0.2, 0.5)); npz["t2z3_orbits_of_cover"] = None
    oid = np.empty(27 * 27, int)
    for k, o in enumerate(orb):
        oid[o] = k
    npz["t2z3_oid"] = oid; npz["t2z3_fixed"] = np.array(fx)
    Hq4, sz4, fx4, orb4 = E.t4z2(6); shapes["T4/Z2"] = (Hq4, fx4, (0.35, 0.7)); npz["t4z2_fixed"] = np.array(fx4)
    shapes["T2_flat"] = (E.hex_torus(15), [], None)
    A, X, F = E.icosphere(3); shapes["S2"] = (A, [], None); npz["ico_X"] = X
    npz["ico_deg5"] = np.where(A.sum(1) == 5)[0]
    for name, (T, fixed, tol) in shapes.items():
        t1 = time.time()
        I = E.mi_of_hopping(T)
        r = E.eye(I, thr2=THR2, keep=True)
        Hop = shortest_path((np.abs(T) > 1e-12).astype(float), unweighted=True, directed=False)
        out = E.strip(r)
        if fixed:
            out["grade_v1"] = grade_cones(np.array(r["cone_sites"], int), r["cone_rho"], Hop, list(fixed), *tol)
            out["grade_v2"] = grade_cones(np.array(r["cone_sites_v2"], int), r["cone_rho_v2"], Hop, list(fixed), *tol)
            out["rho_at_fixed_v1"] = [float(r["rho"][list(r["sites"]).index(f)]) for f in fixed]
            out["rho_at_fixed_v2"] = [float(r["rho_v2"][list(r["sites"]).index(f)]) for f in fixed]
        if name == "S2":
            d5 = npz["ico_deg5"]
            out["rho_at_deg5_v1"] = [float(r["rho"][x]) for x in d5]
            out["rho_at_deg5_v2"] = [float(r["rho_v2"][x]) for x in d5]
        key = name.replace("/", "").replace("2_", "2")
        npz[f"{key}_rho"] = r["rho"]; npz[f"{key}_rho_v2"] = r["rho_v2"]; npz[f"{key}_sites"] = r["sites"]
        out["time_s"] = time.time() - t1
        res[name] = out
        print(name, {k: v for k, v in out.items() if k not in ("cone_sites", "cone_sites_v2")}, flush=True)
    g = {}
    z3, z2, fl, s2 = res["T2/Z3"], res["T4/Z2"], res["T2_flat"], res["S2"]
    for v, ck in (("v1", "cones"), ("v2", "cones_v2")):
        g[v] = dict(T2Z3=bool(z3["cavities"] == 1 and z3["loops"] == 0 and z3[f"grade_{v}"]["ok"]),
                    T2Z3_parts=dict(cavities=z3["cavities"], loops=z3["loops"], cones=z3[ck]),
                    T4Z2_cones=bool(z2[f"grade_{v}"]["ok"] and z2["loops"] == 0),
                    T2_flat_cones0=bool(fl[ck] == 0), S2_cones0=bool(s2[ck] == 0))
    g["T2_flat_torus_recognised"] = bool(fl["loops"] == 2 and fl["cavities"] == 1)
    g["S2_secondary"] = bool(s2["cavities"] == 1 and s2["loops"] == 0)
    g["T4Z2_cavities_not_graded"] = z2["cavities"]
    res["grades"] = g
    res["THR2"] = THR2; res["time_s"] = time.time() - t0
    print(json.dumps(g, indent=1), res["time_s"])
    json.dump(res, open(os.path.join(OUT, "eye_shapes.json"), "w"), indent=1, default=float)
    npz = {k: v for k, v in npz.items() if v is not None}
    np.savez(os.path.join(OUT, "eye_shapes.npz"), **npz)
