"""locale_f4_pair.py - phase 4 (3), exam X1 of lab/out/locale_f4_prereg.json: is the pair 0.8811 / 0.8822 (twisted
Wilson h = 0.1, nucleus on the corner, reduced model of phase 3) the 1s spin doublet?  Same model as
lab/locale_f3_valleys.py (48^3, theta = 0.5, Z = 60, Ritz), initial spinors spin-down and spin-x; spin expectation
<S> = <psi| (1/2) I_chi (x) sigma |psi> / <psi|psi> of each level. Writes lab/out/locale_f4_pair.json (exit 3 = again)."""
import os, sys, json, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locale_core import RealWalk
from locale_cands import fermion
from zz_core import ALPHA

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "locale_f4_pair.json")
T0 = time.time()
R = json.load(open(OUT)) if os.path.exists(OUT) else {}
NB, TH, Z = 48, 0.5, 60
ZA = Z * ALPHA
AB = 1 / (ZA * TH)
NPER, OM = 700, np.arange(0.80, 0.9501, 0.005)
SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0]).astype(complex)]
SPIN = [0.5 * np.kron(np.eye(2), s) for s in SIG]
SPINORS = {"down": np.array([0, 1, 0, 1]) / np.sqrt(2), "x": np.array([1, 1, 1, 1]) / 2.0}


def run(spname):
    w = fermion("F2", 3, TH, h=0.1)
    x = np.arange(NB) - NB // 2
    X, Y, Zc = np.meshgrid(x - 0.5, x - 0.5, x - 0.5, indexing="ij")
    r = np.sqrt(X ** 2 + Y ** 2 + Zc ** 2)
    rw = RealWalk(w, NB, 3, V=-ZA / r * w.T)
    Xi, Yi, Zi = np.meshgrid(x, x, x, indexing="ij")
    cls = (Xi % 2 == 0) & (Yi % 2 == 0) & (Zi % 2 == 0)
    psi0 = SPINORS[spname].astype(complex)[:, None, None, None] * (np.exp(-r / AB) * cls)[None]
    psi0 /= np.linalg.norm(psi0)
    acc = np.zeros((len(OM),) + psi0.shape, complex)
    hann = 0.5 - 0.5 * np.cos(2 * np.pi * (np.arange(NPER) + 0.5) / NPER)
    psi = psi0.copy()
    for t in range(NPER):
        acc += (hann[t] * np.exp(1j * OM * t))[:, None, None, None, None] * psi[None]
        psi = rw.period(psi)
    F = acc.reshape(len(OM), -1); del acc
    G = F.conj() @ F.T
    ev, U = np.linalg.eigh(G)
    keep = ev > 1e-11 * ev.max()
    Q = (U[:, keep] / np.sqrt(ev[keep])).T @ F
    UQ = np.array([rw.period(q.reshape(psi0.shape)).ravel() for q in Q])
    lam, Vv = np.linalg.eig(Q.conj() @ UQ.T)
    out = []
    for i in range(len(lam)):
        v = Vv[:, i] @ Q; nv = np.linalg.norm(v); v /= nv
        res = float(np.linalg.norm(Vv[:, i] @ UQ / nv - lam[i] * v))
        E = float(-np.angle(lam[i]) / (TH * w.T))
        wgt = float(abs(np.vdot(v, psi0.ravel())) ** 2)
        if 0.85 <= E <= 0.95 and wgt >= 0.005:
            vv = v.reshape(4, -1)
            S = [float(np.real(np.einsum("an,ab,bn->", vv.conj(), Sm, vv))) for Sm in SPIN]
            out.append(dict(E=E, weight=wgt, resid=res, S=S))
    return sorted(out, key=lambda d: d["E"])


for sp in ("down", "x"):
    if sp in R:
        continue
    if time.time() - T0 > 30:
        print("rerun"); sys.exit(3)
    R[sp] = run(sp)
    print(sp, [(round(l["E"], 5), round(l["weight"], 3), [round(s, 3) for s in l["S"]], f"{l['resid']:.0e}") for l in R[sp]], flush=True)
    with open(OUT, "w") as f:
        json.dump(R, f, indent=1)
print("ALL DONE")
