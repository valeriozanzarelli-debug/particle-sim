"""redo_r2.py N | diag - R2: the SAME Ollivier-Ricci flow of flow_core on the crumpled torus (case a), plus ONE term:
a pressure from a quantized flux n through the two torus cycles (see redo_prereg.json, R2).
  python redo_r2.py 0|1|2|4    -> lab/out/redo_r2_n{N}.json/.npz   (~2.5 min each)
  python redo_r2.py diag       -> lab/out/redo_r2_diag.json (step-size check of the instability, no flux)
FIX after the first diag run (declared): the eps = 0.1 run stopped after 1 step on the tolerance (eps*kappa ~ 2e-4 <
1e-3), so the diag now runs with the tolerance switched off. n = 4 diverged (see doc).
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import shortest_path
import flow_core as F
from flow_run import read, S0, L

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
G = 0.1


def wraps(m):
    """c[a][e] = 1 if edge e crosses the seam of cycle class a (a=0: index i, a=1: index j)"""
    i, j = np.divmod(np.arange(m.N), L)
    u, v = m.E[:, 0], m.E[:, 1]
    return [(np.abs(i[u] - i[v]) > 1).astype(int), (np.abs(j[u] - j[v]) > 1).astype(int)]


def loops(m, w, C):
    """for each class a and vertex v: length and edge list of the shortest closed loop through v with odd winding along a"""
    N = m.N
    out = []
    for a in range(2):
        c = C[a]
        u, v = m.E[:, 0], m.E[:, 1]
        rows = np.r_[u, u + N, v, v + N]
        cols = np.r_[v + N * c, v + N * (1 - c), u + N * c, u + N * (1 - c)]
        A = sp.csr_matrix((np.r_[w, w, w, w], (rows, cols)), shape=(2 * N, 2 * N))
        D, Pr = shortest_path(A, directed=False, indices=np.arange(N), return_predecessors=True)
        ell = D[np.arange(N), np.arange(N) + N]
        paths = []
        for s in range(N):
            t, ed = s + N, []
            while t != s:
                p = Pr[s, t]
                ed.append(m.eid[tuple(sorted((p % N, t % N)))])
                t = p
            paths.append(ed)
        out.append((ell, paths))
    return out


def pressure(m, w, C, Z0, gn2):
    P = np.zeros(len(m.E))
    if gn2 == 0:
        return P, None
    ls = loops(m, w, C)
    for ell, paths in ls:
        for s, ed in enumerate(paths):
            P[ed] += (L / ell[s]) ** 2
    return gn2 * P / Z0, [float(ls[0][0].min()), float(ls[1][0].min())]


def flux_flow(m, w0, gn2, eps=F.EPS, maxit=F.MAXIT, log=print, tol=F.TOL):
    C = wraps(m)
    w1 = np.ones(len(m.E))
    Z0 = sum(len(ed) for _, paths in loops(m, w1, C) for ed in paths) / len(m.E)
    t0 = time.time()
    w = np.array(w0, float); total = w.sum()
    k = F.ollivier(m, w); k0 = k.copy()
    hist, stop = [], "maxit"
    for it in range(maxit):
        P, sysl = pressure(m, w, C, Z0, gn2)
        hist.append(dict(it=it, k_std=float(k.std()), w_cv=float(w.std() / w.mean()), P_mean=float(P.mean()),
                         P_max=float(P.max()), systole=sysl, frac_long=float((w > 2).mean())))
        wn = w - eps * w * (k - P)
        wn = np.maximum(wn, 1e-6 * total / len(w))
        wn *= total / wn.sum()
        rel = float(np.max(np.abs(wn - w) / w))
        w = wn
        k = F.ollivier(m, w)
        log(f"  it {it:3d} k_std {hist[-1]['k_std']:.4f} w_cv {hist[-1]['w_cv']:.3f} P {hist[-1]['P_mean']:.3f}/"
            f"{hist[-1]['P_max']:.3f} sys {sysl} maxrel {rel:.2e} {time.time()-t0:.0f}s")
        if rel < tol:
            stop = "tol"; break
    hist.append(dict(it=len(hist), k_std=float(k.std()), w_cv=float(w.std() / w.mean()), frac_long=float((w > 2).mean())))
    return w, k0, k, hist, stop, Z0


if __name__ == "__main__":
    arg = sys.argv[1]
    t0 = time.time()
    m = F.hex_mesh(L)
    if arg == "diag":
        rng = np.random.default_rng(S0 + 99)
        w0 = 1 + 0.01 * rng.uniform(-1, 1, len(m.E))
        res = {}
        for eps, steps in ((0.5, 20), (0.1, 100)):
            w, k0, k, hist, stop, _ = flux_flow(m, w0, 0.0, eps=eps, maxit=steps, log=lambda s: None, tol=0.0)
            res[f"eps{eps}"] = dict(steps=len(hist) - 1, stop=stop, flow_time=eps * (len(hist) - 1),
                                    w_cv0=float(w0.std() / w0.mean()), w_cv1=float(w.std() / w.mean()),
                                    k_std0=float(k0.std()), k_std1=float(k.std()),
                                    w_cv_curve=[h["w_cv"] for h in hist])
            print(eps, {k_: v for k_, v in res[f"eps{eps}"].items() if k_ != "w_cv_curve"}, flush=True)
        res["instability_real"] = bool(all(r["w_cv1"] > r["w_cv0"] for r in res.values()))
        res["time_s"] = time.time() - t0
        json.dump(res, open(os.path.join(OUT, "redo_r2_diag.json"), "w"), indent=1)
        print(res["instability_real"], res["time_s"])
        sys.exit()
    n = int(arg)
    gn2 = G * n * n
    w0 = 1 + 0.4 * np.random.default_rng(S0 + 41).uniform(-1, 1, len(m.E))
    eye0, hol0, d0, _, _ = read(m, w0)
    w, k0, k, hist, stop, Z0 = flux_flow(m, w0, gn2)
    eye1, hol1, d1, _, _ = read(m, w)
    ex = dict(one_piece_loops=bool(eye1["pieces"] == 1 and eye1["loops"] == 2),
              spread=bool(k.std() <= k0.std() / 10), spinor=bool(hol1["N_inv"] == 2))
    res = dict(n=n, g=G, gn2=gn2, Z0=Z0, stop=stop, steps=len(hist) - 1, k_std0=float(k0.std()), k_std1=float(k.std()),
               spread_ratio=float(k0.std() / max(k.std(), 1e-15)), w_cv0=float(w0.std() / w0.mean()),
               w_cv1=float(w.std() / w.mean()), frac_long1=float((w > 2).mean()), w_median1=float(np.median(w)),
               eye0=eye0, eye1=eye1, hol0=hol0, hol1=hol1, exams=ex, all_pass=bool(all(ex.values())), hist=hist,
               time_s=time.time() - t0)
    print({k_: v for k_, v in res.items() if k_ not in ("hist", "eye0")}, flush=True)
    json.dump(res, open(os.path.join(OUT, f"redo_r2_n{n}.json"), "w"), indent=1, default=float)
    np.savez(os.path.join(OUT, f"redo_r2_n{n}.npz"), w0=w0, w=w, k0=k0, k=k, d1=d1, pos=m.pos, E=m.E)
