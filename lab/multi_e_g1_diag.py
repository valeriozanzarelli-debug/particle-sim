"""multi_e_g1_diag.py - DIAGNOSTICA POST HOC (dichiarata, non rigiudica il gate) dei due FAIL di G1:
  (1) Hooke omega = 1/2 a K = 200: errore 2.7e-8 contro 1e-8. E' convergenza della base o un errore della macchina?
      -> E(K) per K fino a 400 (stesso seme, stessa ricetta), piu' il pavimento numerico dovuto al condizionamento di S.
  (2) Monte Carlo <1/r12>: un elemento a -1.34e-4 contro 1e-4. E' rumore statistico o la formula?
      -> errore standard del MC per elemento + controllo deterministico indipendente:
         1/r = (2/sqrt pi) int_0^inf exp(-t^2 r^2) dt  => <1/|u|> = (2/sqrt pi) int_0^inf (1 + 2 t^2 sigma2)^{-3/2} dt,
         sigma2 = varianza per componente di u = w^T r (= w^T C^-1 w / 2), integrata con quadratura adattiva.
Run: python lab/multi_e_g1_diag.py hooke | mc"""
import sys, time
from multi_e_ecg import *                    # noqa
from multi_e_ecg import System, SVM, random_A, JS
from multi_e_core import load_json, save_json
import numpy as np
from scipy.integrate import quad

res = load_json(JS)
D = res.get("G1_diag_posthoc", {})
what = sys.argv[1] if len(sys.argv) > 1 else "mc"

if what == "mc":
    rng = np.random.default_rng(5)                       # stesso seme e stessa sequenza di G1(d)
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
        an = np.array([(2 / np.sqrt(np.pi)) * (Ci[0, 0] + Ci[1, 1] - 2 * Ci[0, 1]) ** -0.5, (2 / np.sqrt(np.pi)) * Ci[2, 2] ** -0.5])
        sig2 = np.array([(Ci[0, 0] + Ci[1, 1] - 2 * Ci[0, 1]) / 2, Ci[2, 2] / 2])
        quadv = np.array([(2 / np.sqrt(np.pi)) * quad(lambda tt: (1 + 2 * tt * tt * s) ** -1.5, 0, np.inf, epsabs=0, epsrel=1e-13)[0]
                          for s in sig2])
        recs.append(dict(analytic=an.tolist(), quadrature=quadv.tolist(), quad_rel=(quadv / an - 1).tolist(),
                         mc=mean.tolist(), mc_rel=(mean / an - 1).tolist(), mc_rel_se=(se / an).tolist(),
                         z_scores=((mean - an) / se).tolist()))
        print(recs[-1], flush=True)
    D["MC"] = dict(recs=recs, max_quad_rel=max(abs(x) for r in recs for x in r["quad_rel"]),
                   max_abs_z=max(abs(x) for r in recs for x in r["z_scores"]))
    res["G1_diag_posthoc"] = D; save_json(JS, res)

if what == "hooke":
    s = System(2, 1, 0.0, 1.0, omega=0.5, name="Hooke")
    rng = np.random.default_rng(12)
    sv = SVM(s, rng, 1.0, 20.0, 30)
    rows = []
    t0 = time.time()
    for Kt in (100, 200, 300, 400):
        while len(sv.A) < Kt:
            sv.add()
        for sw in range(3):
            sv.refine_sweep()
        cs = float(np.linalg.cond(sv.Smat))
        Er = sv.rebuild()                                  # ricostruita da zero: elimina errori accumulati
        rows.append(dict(K=Kt, E_minus_2=Er - 2.0, cond_S=cs, sec=time.time() - t0))
        print(rows[-1], flush=True)
        D["Hooke"] = dict(rows=rows, note="stesso seme di G1; ogni K: crescita + 3 giri; post hoc")
        res["G1_diag_posthoc"] = D; save_json(JS, res)
