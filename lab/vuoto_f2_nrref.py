"""vuoto_f2_nrref.py - EXAM REFERENCE ONLY (not the rule): the Bethe sum of non-relativistic hydrogen (Schroedinger,
finite differences, atomic units) with the SAME cutoff definition used on the rule:
    ln k0(K) = ln(K/Ry) - int_0^K S(k) dk / (K S(K)),   S(k) = sum_m |<m|v|n>|^2 dE/(dE + k)
and the full-D version ln(K/Ry) - int S / D (D = 2/n^3 sum rule) whose K -> inf limit is the textbook ln k0
(2.984128556 for 1S, 2.811769893 for 2S). Two grids for the discretization error. Writes lab/out/vuoto_f2_nrref.json.
Run: timeout 590 python lab/vuoto_f2_nrref.py"""
import os, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
KRY = [30, 100, 300, 1000, 3000]


def nr_curve(n, h, R):
    N = int(R / h)
    r = h * np.arange(1, N + 1)

    def H(l):
        main = 1 / h**2 + l * (l + 1) / (2 * r**2) - 1 / r
        off = -0.5 / h**2 * np.ones(N - 1)
        return sp.diags([off, main, off], [-1, 0, 1], format="csc")
    E, V = spl.eigsh(H(0), k=n, sigma=-0.6, which="LM")
    o = np.argsort(E); E = E[o]; V = V[:, o]
    En = E[n - 1]
    u = r * V[:, n - 1]                                   # sum u^2 = <r^2> (V normalized as a vector)
    A = (H(1) - En * sp.identity(N, format="csc")).tocsc()
    uv = A @ u
    D = float(uv @ (A @ uv))
    Ry = 0.5
    Ks = Ry * np.array(KRY, float)
    kg = np.unique(np.concatenate([np.geomspace(1e-7, Ks[-1], 300), Ks]))
    S = np.array([uv @ (A @ spl.spsolve((A + k * sp.identity(N, format="csc")).tocsc(), uv)) for k in kg])
    lk = np.log(kg)
    cum = kg[0] * (uv @ uv) + np.concatenate([[0], np.cumsum(0.5 * (kg[1:] * S[1:] + kg[:-1] * S[:-1]) * np.diff(lk))])
    out = {}
    for K, kr in zip(Ks, KRY):
        i = int(np.argmin(np.abs(kg - K)))
        out[str(kr)] = {"lnk0_KS": float(np.log(K / Ry) - cum[i] / (K * S[i])), "D_K_over_sumrule": float(K * S[i] / (2.0 / n**3)),
                        "lnk0_fullD": float(np.log(K / Ry) - cum[i] / (2.0 / n**3))}
    return {"E_n": float(En), "D_over_sumrule": D / (2.0 / n**3), "K": out}


R = {}
for n, Rbox in ((1, 60.0), (2, 120.0)):
    for h in (0.02, 0.01):
        R["%d|%g" % (n, h)] = nr_curve(n, h, Rbox)
        print(n, h, json.dumps(R["%d|%g" % (n, h)]["K"]), flush=True)
json.dump(R, open(os.path.join(OUT, "vuoto_f2_nrref.json"), "w"), indent=1)
