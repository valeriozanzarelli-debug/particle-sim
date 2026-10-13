"""aggfull_p5_analyze.py - verdicts of PHASE 5 (lab/out/aggfull_p5_prereg.json) -> lab/out/aggfull_p5_verdict.json"""
import os, json
import numpy as np
import aggfull_core as ac
from zz_core import dirac_formula
OUT = ac.OUT
P = json.load(open(os.path.join(OUT, "aggfull_p5.json")))
A = json.load(open(os.path.join(OUT, "aggfull_atom.json")))
AK = json.load(open(os.path.join(OUT, "aggfull_atomK.json")))
AT = json.load(open(os.path.join(OUT, "aggfull_atomT.json")))
V, rep = [], {}
def check(n, ok, msg):
    V.append((n, bool(ok), msg)); print(f"[{'PASS' if ok else 'FAIL'}] {n}: {msg}", flush=True)
s = [P[f"s{i}"] for i in range(1, 7)]
check("M2c Z = 10 level (auto_level rule) passes 3e-5 on an independent check, seeds 1-6",
      all(x["w"] is not None and x["exam_ok"] for x in s),
      "w = " + str([x["w"] for x in s]) + "; independent energy " + ", ".join(f"{x['exam_check']['energy_rel']:.2e}" for x in s))
fl = lambda f, ax: np.roll(np.flip(f, axis=ax), 1, axis=ax)
def proj(v):
    v = 0.5 * (v + fl(v, 0)); v = 0.5 * (v + fl(v, 1)); return 0.5 * (v - fl(v, 2))
N = 256
k1 = 2 * np.pi * np.fft.fftfreq(N)
KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij", sparse=True)
hi = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2) > 1.6
c1 = {}
for job, R in (("e8_2p_c", A), ("eK8_2p_c", AK), ("eT8_2p_c", AT)):
    v = np.load(os.path.join(OUT, f"aggfull_vec_{job}.npy")).reshape((N,) * 3); v /= np.linalg.norm(v)
    j = float(np.sum((v - proj(v)) ** 2))
    w = np.abs(ac.fftn(v)) ** 2; w /= w.sum()
    c1[job] = dict(junk=j, tail=float(w[hi].sum()), residual=R[job]["residual"], E=R[job]["E_eV"])
    del v, w
rep["C1"] = c1
check("C1 clean solver: no wrong-parity junk, no high-k tail, residual < 2e-9",
      all(x["junk"] <= 1e-20 and x["tail"] <= 1e-8 and x["residual"] < 2e-9 for x in c1.values()),
      "; ".join(f"{k}: junk {x['junk']:.1e}, tail {x['tail']:.1e}, residual {x['residual']:.1e}" for k, x in c1.items()))
d12 = -(1 - dirac_formula(2, 1)) * ac.MC2 + ac.RY / 4
d32 = -(1 - dirac_formula(2, -2)) * ac.MC2 + ac.RY / 4
f3 = {}
for nm, q in (("orig", A["q8_2p_c"]), ("K", AK["qK8_2p_c"]), ("T", AT["qT8_2p_c"])):
    ev = np.array(q["dV_eigen_eV"]); R = q["kinetic_remainder_eV"]
    f3[nm] = dict(R=R, e12=(R + ev[:2].mean()) / d12 - 1, e32=(R + ev[2:].mean()) / d32 - 1)
rep["C2"] = f3
check("C2 F3 of the 2p on the whole cube with the T kinetic within 5 %", abs(f3["T"]["e12"]) <= 0.05 and abs(f3["T"]["e32"]) <= 0.05,
      "; ".join(f"{k}: 2p1/2 {100 * x['e12']:+.2f} %, 2p3/2 {100 * x['e32']:+.2f} % (R {x['R']:.4e} eV)" for k, x in f3.items())
      + f"; exact mass-velocity 2p -2.64e-5 eV")
pairs = (("e8_2p", A), ("eK8_2p", AK))
rel = {j: c1[j + "_c"]["E"] / R[j]["E_eV"] - 1 for j, R in pairs}
check("C3 clean 2p energies equal the old within 1e-6", all(abs(x) <= 1e-6 for x in rel.values()),
      "; ".join(f"{j}: {x:+.2e}" for j, x in rel.items()) + " (the junk norm^2 3.8e-6 dilutes the old Rayleigh quotient by the same amount)")
rep["verdict"] = V
json.dump(rep, open(os.path.join(OUT, "aggfull_p5_verdict.json"), "w"), indent=1, default=float)
print(f"=== {sum(v[1] for v in V)} PASS / {sum(not v[1] for v in V)} FAIL")
