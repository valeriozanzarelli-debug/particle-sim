"""aggfull_p4_f3.py - F3P of PHASE 4 as amended (lab/out/aggfull_p4_amend.json): alpha^4 shifts of 2p1/2 and 2p3/2 at
a0/8 with three kinetics (original chain A, K = amendments 1-2, T = phase 4), kinetic remainder on |kappa| < 0.8 and on
the whole cube. Writes lab/out/aggfull_p4_f3.json"""
import os, json, pickle
import numpy as np
import aggfull_core as ac
import aggfull_atom as aa
from zz_core import dirac_formula

OUT = ac.OUT
N, B = 256, 512
h = B * ac.TH_A
A = json.load(open(os.path.join(OUT, "aggfull_atom.json")))
AK = json.load(open(os.path.join(OUT, "aggfull_atomK.json")))
recs = {"orig": [l for l in ac.load_chain_A() if l["B"] == B][0],
        "K": [l for l in pickle.load(open(os.path.join(OUT, "aggfull_relearn.pkl"), "rb")) if l["B"] == B][0],
        "T": [l for l in pickle.load(open(os.path.join(OUT, "aggfull_relearnT.pkl"), "rb")) if l["B"] == B][0]}
src = {"orig": ("e8_2p", A["q8_2p"]), "K": ("eK8_2p", AK["qK8_2p"]), "T": ("eK8_2p", AK["qK8_2p"])}
k1 = 2 * np.pi * np.fft.fftfreq(N)
KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij", sparse=True)
inside = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2) < 0.8
d12 = -(1 - dirac_formula(2, 1)) * ac.MC2 + ac.RY / 4
d32 = -(1 - dirac_formula(2, -2)) * ac.MC2 + ac.RY / 4
res = dict(dirac_shift_2p12=d12, dirac_shift_2p32=d32)
for nm, rec in recs.items():
    job, q = src[nm]
    v = np.load(os.path.join(OUT, f"aggfull_vec_{job}.npy")).reshape((N,) * 3)
    w = np.abs(ac.fftn(v)) ** 2; del v
    w /= w.sum()
    T, _, _ = aa.kinetic(rec, N, h)
    mk, _ = aa.kinetic_mass(rec, h)
    Rk = (T - aa.kgrid2(N, h) / (2 * mk)) * ac.MC2
    del T
    R_all = float(np.sum(w * Rk))
    R_in = float(np.sum(w[inside] * Rk[inside]) / w[inside].sum())
    del Rk
    ev = np.array(q["dV_eigen_eV"])
    lo, up = ev[:2].mean(), ev[2:].mean()
    out = dict(vector=job, m_kin=mk, R_all_eV=R_all, R_in_eV=R_in, weight_outside=float(1 - w[inside].sum()))
    for tag, R in (("in", R_in), ("all", R_all)):
        s12, s32 = R + lo, R + up
        out[f"shift12_{tag}"], out[f"shift32_{tag}"] = s12, s32
        out[f"err12_{tag}"], out[f"err32_{tag}"] = s12 / d12 - 1, s32 / d32 - 1
    res[nm] = out
    print(nm, {k: (f"{x:.4e}" if isinstance(x, float) else x) for k, x in out.items()}, flush=True)
ok = {nm: abs(res[nm]["err12_in"]) <= 0.05 and abs(res[nm]["err32_in"]) <= 0.05 for nm in recs}
res["F3P_pass_by_kinetic"] = ok
res["F3P_PASS"] = any(ok.values())
json.dump(res, open(os.path.join(OUT, "aggfull_p4_f3.json"), "w"), indent=1, default=float)
print("F3P", ok, "->", "PASS" if res["F3P_PASS"] else "FAIL")
