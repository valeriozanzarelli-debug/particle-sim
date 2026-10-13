"""agg2s_posthoc_3s.py - POST-HOC diagnostic (declared, after the BORN_3s FAIL): is the 3s energy residual at L = 48 a0
(twist pair 000/111: +0.011 %) the box? Adds the twists 100 and 110 (1s, 2s, 3s; exact Schroedinger, cells a0/4) and
forms the full 8-twist average (E000 + 3 E100 + 3 E110 + E111)/8, which also cancels the EDGES of the image crystal.
Writes lab/out/agg2s_posthoc_3s.json"""
import os, sys, json, time
import numpy as np
import agg2s_core as c

JS = os.path.join(c.OUT, "agg2s_posthoc_3s.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
CT = json.load(open(os.path.join(c.OUT, "agg2s_contact.json")))
DIV, N = 4, 192
h = c.A0 / DIV
t0 = time.time()
f = c.galerkin_field(N, h)
for tws in ("100", "110"):
    tw = tuple(int(x) for x in tws)
    at = c.TwistAtom(c.kin_schroedinger(N, h, tw), f, tw, sym="even")
    vecs = []
    for st in ("1s", "2s", "3s"):
        key = f"S{st}_{tws}"
        vf = os.path.join(c.OUT, f"agg2s_cvec_{key}.npy")
        if res.get(key, {}).get("converged"):
            vecs.append(np.load(vf)); continue
        if time.time() - t0 > 380:
            sys.exit(0)
        base = np.load(os.path.join(c.OUT, f"agg2s_cvec_S{st}_000.npy")).reshape((N,) * 3)
        x0 = np.load(vf) if os.path.exists(vf) else at.project(base).ravel()
        v, E, rr, it, ok = c.lobpcg1(at, x0, Y=np.stack(vecs, 1) if vecs else None, tol_res=1e-9,
                                     budget=max(60, 540 - (time.time() - t0)), log=lambda s: None)
        np.save(vf, v)
        n = int(st[0])
        res[key] = dict(E_eV=E * c.MC2, E_rel=E * c.MC2 / (-c.RY / n ** 2) - 1, residual=rr, converged=bool(ok))
        json.dump(res, open(JS, "w"), indent=1)
        print(key, res[key], flush=True)
        if not ok:
            sys.exit(0)
        vecs.append(v)
if all(res.get(f"S{s}_{t}", {}).get("converged") for s in ("1s", "2s", "3s") for t in ("100", "110")):
    for s in ("1s", "2s", "3s"):
        n = int(s[0])
        E = [CT[f"S{s}_000"]["E_eV"], res[f"S{s}_100"]["E_eV"], res[f"S{s}_110"]["E_eV"], CT[f"S{s}_111"]["E_eV"]]
        full = (E[0] + 3 * E[1] + 3 * E[2] + E[3]) / 8
        rhoc = 0.5 * (CT[f"S{s}_000"]["rho_c_a0"] + CT[f"S{s}_111"]["rho_c_a0"])
        J = c.born_constants()["J_coef"]
        dEB = -2 * (4 * np.pi) ** 2 * rhoc * J * (1 / DIV) ** 3 * 2 * c.RY
        res[f"full_{s}"] = dict(E_full_rel=full / (-c.RY / n ** 2) - 1, E_pair_rel=0.5 * (E[0] + E[3]) / (-c.RY / n ** 2) - 1,
                                E_full_born_rel=(full + dEB) / (-c.RY / n ** 2) - 1, steps_rel=list(np.diff(E) / (c.RY / n ** 2)))
        print(s, res[f"full_{s}"], flush=True)
    json.dump(res, open(JS, "w"), indent=1)
