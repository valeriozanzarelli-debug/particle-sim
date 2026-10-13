"""agg2s_contact.py - production: exact Schroedinger kinetic + Galerkin field, cells a0/4, L = 48 a0 (N = 192), twists
000 and 111, states 1s, 2s, 3s (the 3s is the BLIND test of the Born tail / contact universality, agg2s_prereg.json).
Usage: timeout 590 python lab/agg2s_contact.py   (resumes; results lab/out/agg2s_contact.json, vectors not committed)"""
import os, sys, json, time
import numpy as np
import agg2s_core as c

JS = os.path.join(c.OUT, "agg2s_contact.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
DIV, N = 4, 192
h = c.A0 / DIV
t0 = time.time()
f = c.galerkin_field(N, h)
for tws in ("000", "111"):
    tw = tuple(int(x) for x in tws)
    at = c.TwistAtom(c.kin_schroedinger(N, h, tw), f, tw, sym="even")
    vecs = []
    for st in ("1s", "2s", "3s"):
        key = f"S{st}_{tws}"
        vf = os.path.join(c.OUT, f"agg2s_cvec_{key}.npy")
        if res.get(key, {}).get("converged"):
            vecs.append(np.load(vf))
            continue
        if time.time() - t0 > 400:
            sys.exit(0)
        if os.path.exists(vf):
            x0 = np.load(vf)
        elif st == "3s":
            X, Yg, Z = np.meshgrid(f.x, f.x, f.x, indexing="ij", sparse=True)
            r = np.sqrt(X ** 2 + Yg ** 2 + Z ** 2) / c.A0
            x0 = at.project(np.broadcast_to((1 - 2 * r / 3 + 2 * r ** 2 / 27) * np.exp(-r / 3), (N,) * 3).copy()).ravel()
        else:
            x0 = c.guess(at, st)
        Y = np.stack(vecs, 1) if vecs else None
        v, E, rr, it, ok = c.lobpcg1(at, x0, Y=Y, tol_res=1e-9, budget=max(60, 540 - (time.time() - t0)),
                                     log=lambda s: print(key, s, flush=True), every=10)
        np.save(vf, v)
        vv = v.reshape((N,) * 3)
        h_a0 = 1.0 / DIV
        nrm = float(np.sum(vv ** 2)) * h_a0 ** 3
        rho = float(vv[N // 2, N // 2, N // 2] ** 2 / nrm)
        F = np.fft.fftn(vv)
        ka = c.kgrid_axes(N, h_a0, tw)
        k2 = ka[0][:, None, None] ** 2 + ka[1][None, :, None] ** 2 + ka[2][None, None, :] ** 2
        w = np.abs(F) ** 2
        k4 = float(np.sum(w * k2 ** 2) / np.sum(w))
        n = int(st[0])
        res[key] = dict(E_eV=E * c.MC2, E_rel=E * c.MC2 / (-c.RY / n ** 2) - 1, residual=rr, converged=bool(ok), rho_c_a0=rho, k4_c_a0=k4,
                        h_a0=h_a0, L_a0=N * h_a0, iterations=res.get(key, {}).get("iterations", 0) + it)
        json.dump(res, open(JS, "w"), indent=1)
        print(key, res[key], flush=True)
        if not ok:
            sys.exit(0)
        vecs.append(v)
