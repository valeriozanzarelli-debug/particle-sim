"""multi_e_g1v2.py - gate G1' (pre-registrazione v2, lab/out/multi_e_prereg_v2.json).
  Hooke omega = 1/2, seme 12: crescita SVM a K = 100, 200, 300 con 3 giri a ogni tappa, ricostruzione: |E(300) - 2| < 1e-8
  MC (seme 5, 3 casi come v1): quadratura indipendente = formula a 1e-12 relativo  E  |z| < 3 con 1e8 campioni.
  H, Li_nonint, Be_nonint: PASS della v1 (stesso codice degli elementi).
Run: python lab/multi_e_g1v2.py"""
import time
import numpy as np
from scipy.integrate import quad
from multi_e_ecg import System, SVM, random_A, JS
from multi_e_core import load_json, save_json

res = load_json(JS)
out = res.get("G1_v2", {})
t0 = time.time()
if "Hooke" not in out:
    s = System(2, 1, 0.0, 1.0, omega=0.5, name="Hooke")
    sv = SVM(s, np.random.default_rng(12), 1.0, 20.0, 30)
    rows = []
    for Kt in (100, 200, 300):
        while len(sv.A) < Kt:
            sv.add()
        for sw in range(3):
            sv.refine_sweep()
        cs = float(np.linalg.cond(sv.Smat))
        Er = sv.rebuild()
        rows.append(dict(K=Kt, E_minus_2=Er - 2.0, cond_S=cs, sec=time.time() - t0))
        print(rows[-1], flush=True)
    err = rows[-1]["E_minus_2"]
    out["Hooke"] = dict(rows=rows, err=err, PASS=bool(abs(err) < 1e-8))
    res = load_json(JS); res["G1_v2"] = out; save_json(JS, res)
if "MC" not in out:
    rng = np.random.default_rng(5)
    N = 3
    recs = []
    for t in range(3):
        A = random_A(rng, N, 2.0); B = random_A(rng, N, 2.0)
        Ci = np.linalg.inv(A + B)
        Lc = np.linalg.cholesky(Ci / 2)
        s1 = np.zeros(2); s2 = np.zeros(2); n = 0
        for ch in range(40):
            z = rng.normal(size=(N, 3, 10 ** 7 // 4))
            r = np.einsum("ij,jkn->ikn", Lc, z)
            x = np.stack([1 / np.linalg.norm(r[0] - r[1], axis=0), 1 / np.linalg.norm(r[2], axis=0)])
            s1 += x.sum(1); s2 += (x ** 2).sum(1); n += x.shape[1]
        mean = s1 / n; se = np.sqrt((s2 / n - mean ** 2) / n)
        g = np.array([Ci[0, 0] + Ci[1, 1] - 2 * Ci[0, 1], Ci[2, 2]])
        an = (2 / np.sqrt(np.pi)) * g ** -0.5
        qd = np.array([(2 / np.sqrt(np.pi)) * quad(lambda tt: (1 + 2 * tt * tt * s) ** -1.5, 0, np.inf, epsabs=0, epsrel=1e-13,
                                                   limit=200)[0] for s in g / 2])
        recs.append(dict(analytic=an.tolist(), quadrature=qd.tolist(), quad_rel=(qd / an - 1).tolist(), mc=mean.tolist(),
                         z=((mean - an) / se).tolist(), n=n))
        print(recs[-1], flush=True)
    mq = max(abs(x) for r_ in recs for x in r_["quad_rel"]); mz = max(abs(x) for r_ in recs for x in r_["z"])
    out["MC"] = dict(recs=recs, max_quad_rel=mq, max_abs_z=mz, PASS=bool(mq < 1e-12 and mz < 3))
    res = load_json(JS); res["G1_v2"] = out; save_json(JS, res)
v1 = res["G1"]
out["PASS"] = bool(out["Hooke"]["PASS"] and out["MC"]["PASS"] and all(v1[k]["PASS"] for k in ("H", "Li_nonint", "Be_nonint")))
res = load_json(JS); res["G1_v2"] = out; save_json(JS, res)
print("G1':", "PASS" if out["PASS"] else "FAIL", "Hooke", out["Hooke"]["err"], "MC quad", out["MC"]["max_quad_rel"], "z", out["MC"]["max_abs_z"])
