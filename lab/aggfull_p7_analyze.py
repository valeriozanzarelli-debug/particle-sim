"""aggfull_p7_analyze.py - verdicts of PHASE 7 (lab/out/aggfull_p7_prereg.json) -> lab/out/aggfull_p7_verdict.json"""
import os, sys, json, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_core import dirac_binding, ALPHA       # read only
import aggfull_l1 as l1
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "aggfull_p7_prereg.json")))
E = P["exams"]
G = json.load(open(os.path.join(OUT, "aggfull_p7_glue_err.json")))
C = json.load(open(os.path.join(OUT, "aggfull_p7_chain.json")))
TH = (1 / ALPHA) / 4096
V, rep = [], {}


def check(n, ok, msg):
    V.append((n, bool(ok), msg)); print(f"[{'PASS' if ok else 'FAIL'}] {n}: {msg}", flush=True)


def dpair_err(B1, B2, n):
    return ((B2 / 4 - B1) / 3) / ((dirac_binding(n, -1, 2 * ALPHA) / 4 - dirac_binding(n, -1, ALPHA)) / 3) - 1


D4 = {1: ALPHA ** 4 / 8, 2: 5 * ALPHA ** 4 / 128}       # Z^4 coefficient of the Dirac binding (mc^2)
fits = {}
for v in ("c1", "c2"):
    tf = [x for x in P["config"]["chain"]["variants"] if x[0] == v][0][1]
    m = np.tan(TH * tf) / (TH * tf)
    for n in (1, 2):
        Z = np.array(P["config"]["chain"]["Z"], float)
        y = np.array([C[f"{v}_Z{int(z)}_n{n}"]["B"] / m - dirac_binding(n, -1, z * ALPHA) for z in Z])
        A = np.stack([Z ** 3, Z ** 4], 1); e3, e4 = np.linalg.lstsq(A, y, rcond=None)[0]
        A2 = np.stack([Z ** 2, Z ** 3, Z ** 4, Z ** 6], 1); c4 = np.linalg.lstsq(A2, y, rcond=None)[0]
        two = dpair_err(C[f"{v}_Z1_n{n}"]["B"], C[f"{v}_Z2_n{n}"]["B"], n)
        fits[f"{v}_{n}s"] = dict(e3=e3, e4=e4, z3_part=(e3 / 3) / D4[n], z4_part=e4 / D4[n],
                                 two_Z=two, full_fit=c4.tolist(), z4_part_full=c4[2] / D4[n])
rep["chain_fits"] = fits
two_chain = fits["c1_1s"]["two_Z"]
check("P1 chain alone = glue plateau (1s, theta, J0 = 64)", abs(two_chain - G["th1_J64_1s"]) <= E["P1"]["tol"],
      f"chain {100 * two_chain:+.2f} %, glue {100 * G['th1_J64_1s']:+.2f} %")
check("P2 chain cell theta/4: 1s, 2s within +-1 % at J0 = 24-64", all(abs(G[f"th4_J{j}_{n}s"]) <= E["P2"]["tol"] for j in (24, 32, 48, 64) for n in (1, 2)),
      "; ".join(f"J{j}: " + " / ".join(f"{100 * G[f'{t}_J{j}_1s']:+.2f}" for t in ("th1", "th2", "th4")) + " (1s, theta/1,2,4), "
                + " / ".join(f"{100 * G[f'{t}_J{j}_2s']:+.2f}" for t in ("th1", "th2", "th4")) + " (2s)" for j in (24, 32, 48, 64)))
f = fits["c1_1s"]
# (fixed after the first run: e4 is already an ERROR (y = B/m - B_Dirac), its part is e4/D4, not e4/D4 - 1; declared in s.16)
recon = f["z3_part"] + f["z4_part"]
check("P3 term by term (theta, 1s): Z^3 part + Z^4 part = two-Z plateau; Z^4 residual <= 2 %",
      abs(recon - f["two_Z"]) <= E["P3"]["tol_point"] and abs(f["z4_part"]) <= E["P3"]["tol_res"],
      f"two-Z {100 * f['two_Z']:+.2f} % = Z^3 term {100 * f['z3_part']:+.2f} % + Z^4 residual {100 * f['z4_part']:+.2f} % "
      f"(sum {100 * recon:+.2f}); with Z^2, Z^6 free: Z^4 residual {100 * f['z4_part_full']:+.2f} %; 2s: Z^3 {100 * fits['c1_2s']['z3_part']:+.2f} %, "
      f"Z^4 {100 * fits['c1_2s']['z4_part']:+.2f} % (two-Z {100 * fits['c1_2s']['two_Z']:+.2f} %)")
r = fits["c1_2s"]["e3"] / fits["c1_1s"]["e3"]
check("P4 e3(2s)/e3(1s) = 1/8 (contact-like defect)", abs(r * 8 - 1) <= E["P4"]["tol"], f"{r:.4f} (1/8 = 0.125)")
r5 = fits["c2_1s"]["e3"] / fits["c1_1s"]["e3"]
check("P5 |e3(theta/2)| <= |e3(theta)|/4", abs(r5) <= 0.25, f"ratio {r5:.4f}; 2s {fits['c2_2s']['e3'] / fits['c1_2s']['e3']:.4f}")
k1 = {}
for rec in pickle.load(open(os.path.join(OUT, "aggfull_relearnT.pkl"), "rb")):
    B = rec["B"]; h = B * TH; offs = np.asarray(rec["offs"], float); K = rec["K"]
    g = np.real(0.5 * (K[0, 0] + K[1, 1]))
    D = l1.dirs13()
    c4 = np.array([np.sum(g * (offs @ d) ** 4) / 24 for d in D])
    dirac = -h * (1 / 8) / h ** 4                         # G = tau T(k/h), tau = h, T = -k^4/8 (m = 1)
    idx = {tuple(o): i for i, o in enumerate(np.asarray(rec["offs"]).astype(int))}
    odd = max(float(np.max(np.abs(K[:, :, i] - K[:, :, idx[tuple(-np.asarray(o))]]))) for o, i in idx.items() if tuple(-np.asarray(o)) in idx)
    k1[B] = dict(c4_over_dirac_min=float((c4 / dirac).min()), c4_over_dirac_max=float((c4 / dirac).max()), parity_odd_max=odd)
    for rk in pickle.load(open(os.path.join(OUT, "aggfull_relearn.pkl"), "rb")):
        if rk["B"] == B:
            gk = np.real(0.5 * (rk["K"][0, 0] + rk["K"][1, 1])); ok_ = np.asarray(rk["offs"], float)
            ck = np.array([np.sum(gk * (ok_ @ d) ** 4) / 24 for d in D]) / dirac
            k1[B].update(K_c4_min=float(ck.min()), K_c4_max=float(ck.max()), K_c4_mean=float(ck.mean()))
rep["K1"] = k1
check("K1 level T = K with c4 pinned: c4 within 0.5 % of Dirac in 13 directions, parity-odd <= 1e-10",
      all(abs(x["c4_over_dirac_min"] - 1) <= 0.005 and abs(x["c4_over_dirac_max"] - 1) <= 0.005 and x["parity_odd_max"] <= 1e-10 for x in k1.values()),
      "; ".join(f"B={B}: T c4/Dirac {x['c4_over_dirac_min']:.4f}-{x['c4_over_dirac_max']:.4f}, odd {x['parity_odd_max']:.1e}; "
                f"K {x['K_c4_min']:.3f}-{x['K_c4_max']:.3f} (mean {x['K_c4_mean']:.3f})" for B, x in k1.items()))
rep["glue"] = G
rep["verdicts"] = [dict(exam=n, passed=o, detail=m) for n, o, m in V]
json.dump(rep, open(os.path.join(OUT, "aggfull_p7_verdict.json"), "w"), indent=1, default=float)
print(f"{sum(o for _, o, _ in V)} PASS / {sum(not o for _, o, _ in V)} FAIL")
