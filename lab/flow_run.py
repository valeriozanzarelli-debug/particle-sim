"""flow_run.py CASE - L2: one normalised Ollivier-Ricci flow run, read before/after by the L1 eye and the Regge
spinor-holonomy check. CASE in a, b1, b2, c1, c2, z3 (see flow_prereg.json).
Writes lab/out/flow_CASE.json and flow_CASE.npz."""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.sparse.csgraph import shortest_path
import flow_core as F
import eye_core as E

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
THR2 = json.load(open(os.path.join(OUT, "eye_calib.json")))["THR2"]
S0 = 20261009
L = 15


def z3_lengths(m):
    i, j = np.divmod(np.arange(m.N), L)
    u, v = m.E[:, 0], m.E[:, 1]
    di = (i[v] - i[u] + 1) % L - 1; dj = (j[v] - j[u] + 1) % L - 1
    mi, mj = i[u] + di / 2, j[u] + dj / 2
    f = np.cos(2 * np.pi * mi / L) + np.cos(2 * np.pi * mj / L) + np.cos(2 * np.pi * (mi + mj) / L)
    return 1 + 0.2 * f


def read(m, w):
    r = E.eye(E.mi_of_hopping(F.hopping(m, w)), thr2=THR2, keep=True)
    hol, d = F.holonomy(m, w)
    return E.strip(r), hol, d, r["rho"], r["rho_v2"]


if __name__ == "__main__":
    case = sys.argv[1]
    t0 = time.time()
    if case == "a":
        m = F.hex_mesh(L); w0 = 1 + 0.4 * np.random.default_rng(S0 + 41).uniform(-1, 1, len(m.E))
    elif case in ("b1", "b2"):
        m = F.torus_delaunay(225, np.random.default_rng(S0 + 41 + int(case[1]))); w0 = np.ones(len(m.E))
    elif case in ("c1", "c2"):
        m = F.sphere_hull(225, np.random.default_rng(S0 + 43 + int(case[1]))); w0 = np.ones(len(m.E))
    elif case == "z3":
        m = F.hex_mesh(L); w0 = z3_lengths(m)
    print(case, m.kind, "V", m.N, "E", len(m.E), "T", len(m.T), "euler", m.euler, flush=True)
    eye0, hol0, d0, rho0, rho0v2 = read(m, w0)
    print("before eye", {k: eye0[k] for k in ("loops", "cavities", "rel1", "rel2", "cones", "cones_v2")}, "hol", hol0,
          flush=True)
    w, k0, k, hist, stop = F.flow(m, w0)
    eye1, hol1, d1, rho1, rho1v2 = read(m, w)
    print("after eye", {k_: eye1[k_] for k_ in ("loops", "cavities", "rel1", "rel2", "cones", "cones_v2")}, "hol", hol1,
          flush=True)
    res = dict(case=case, kind=m.kind, V=m.N, E=len(m.E), euler=m.euler, stop=stop, steps=len(hist) - 1,
               k_mean0=float(k0.mean()), k_std0=float(k0.std()), k_mean1=float(k.mean()), k_std1=float(k.std()),
               spread_ratio=float(k0.std() / max(k.std(), 1e-15)), w_cv0=float(w0.std() / w0.mean()),
               w_cv1=float(w.std() / w.mean()), eye0=eye0, eye1=eye1, hol0=hol0, hol1=hol1, hist=hist, THR2=THR2)
    vc0, vc1 = F.vertex_curv(m, k0), F.vertex_curv(m, k)
    if case == "z3":
        fixed = [0 * L + 0, 5 * L + 5, 10 * L + 10]
        A = np.zeros((m.N, m.N)); A[m.E[:, 0], m.E[:, 1]] = 1; A[m.E[:, 1], m.E[:, 0]] = 1
        Hop = shortest_path(A, unweighted=True, directed=False)
        near = np.min(Hop[:, fixed], axis=1) <= 1
        conc = lambda x: float(np.abs(x)[near].sum() / np.abs(x).sum())
        res["z3"] = dict(area_baseline=float(near.mean()), C_deficit0=conc(d0), C_deficit1=conc(d1),
                         C_kappa0=conc(vc0), C_kappa1=conc(vc1), maxdef0=float(np.nanmax(np.abs(d0))),
                         maxdef1=float(np.nanmax(np.abs(d1))))
        for v, key in (("v1", "cone_sites"), ("v2", "cone_sites_v2")):
            cs = eye1[key]
            res["z3"][f"eye_cones_at_fixed_{v}"] = bool(len(cs) == 3 and
                                                        len({int(np.argmin(Hop[c, fixed])) for c in cs}) == 3 and
                                                        all(Hop[c, fixed].min() <= 1 for c in cs))
        z = res["z3"]
        z["verdict"] = ("cone-like" if z["C_deficit1"] >= 0.5 and (z["eye_cones_at_fixed_v1"] or z["eye_cones_at_fixed_v2"])
                        else "washed out" if z["maxdef1"] <= z["maxdef0"] / 10 else "neither")
        print("z3", z, flush=True)
    res["time_s"] = time.time() - t0
    print({k_: v for k_, v in res.items() if k_ not in ("hist", "eye0", "eye1")}, flush=True)
    json.dump(res, open(os.path.join(OUT, f"flow_{case}.json"), "w"), indent=1, default=float)
    np.savez(os.path.join(OUT, f"flow_{case}.npz"), w0=w0, w=w, k0=k0, k=k, d0=d0, d1=d1, vc0=vc0, vc1=vc1,
             pos=m.pos, E=m.E, T=m.T, rho0=rho0, rho1=rho1, rho0v2=rho0v2, rho1v2=rho1v2)
