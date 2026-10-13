"""locale_f5_split.py - phase 5, T1 of lab/out/locale_f5_prereg.json: scaling with the lattice step of the 1s spin-doublet
split (broken time reversal) in the reduced model (Z = 60, nucleus on the corner, twisted Wilson h = 0.1, spin-up start on
one copy, Ritz on the 2-tick operator). theta in {0.7, 0.5, 0.35, 0.25}; box 48^3 (64^3 for 0.25); periods 350/theta.
Checkpoints of the long filter runs go to the scratchpad (never to the repo). Writes lab/out/locale_f5_split.json."""
import os, sys, json, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locale_core import RealWalk
from locale_cands import fermion
from zz_core import ALPHA

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "locale_f5_split.json")
CK = os.environ.get("LOCALE_CK", "/tmp/claude-0/-home-claude-particle-sim/b264be98-f15d-5b13-bc01-3d629f2d0f9f/scratchpad")
T0 = time.time()
R = json.load(open(OUT)) if os.path.exists(OUT) else {}
Z = 60
ZA = Z * ALPHA
SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0]).astype(complex)]
SPIN = [0.5 * np.kron(np.eye(2), s) for s in SIG]
CONF = {0.7: (48, 0.005), 0.5: (48, 0.005), 0.35: (48, 0.005), 0.25: (64, 0.01)}


def setup(th):
    NB, step = CONF[th]
    w = fermion("F2", 3, th, h=0.1)
    x = np.arange(NB) - NB // 2
    X, Y, Zc = np.meshgrid(x - 0.5, x - 0.5, x - 0.5, indexing="ij")
    r = np.sqrt(X ** 2 + Y ** 2 + Zc ** 2)
    rw = RealWalk(w, NB, 3, V=-ZA / r * w.T)
    Xi, Yi, Zi = np.meshgrid(x, x, x, indexing="ij")
    cls = (Xi % 2 == 0) & (Yi % 2 == 0) & (Zi % 2 == 0)
    AB = 1 / (ZA * th)
    psi0 = (np.array([1, 0, 1, 0]) / np.sqrt(2)).astype(complex)[:, None, None, None] * (np.exp(-r / AB) * cls)[None]
    psi0 /= np.linalg.norm(psi0)
    Egrid = np.arange(0.80, 0.9701, step)
    OM = Egrid * th * w.T                      # per-period phase
    return w, rw, psi0, OM, int(round(350 / th))


def run(th):
    w, rw, psi0, OM, NPER = setup(th)
    ck = os.path.join(CK, f"locale_f5_split_{th}.npz")
    if os.path.exists(ck):
        d = np.load(ck); psi, acc, t0 = d["psi"], d["acc"], int(d["t"])
    else:
        psi, acc, t0 = psi0.copy(), np.zeros((len(OM),) + psi0.shape, complex), 0
    hann = 0.5 - 0.5 * np.cos(2 * np.pi * (np.arange(NPER) + 0.5) / NPER)
    for t in range(t0, NPER):
        if time.time() - T0 > 500:
            np.savez(ck, psi=psi, acc=acc, t=t); print(f"theta {th}: checkpoint at period {t}/{NPER}", flush=True)
            sys.exit(3)
        acc += (hann[t] * np.exp(1j * OM * t))[:, None, None, None, None] * psi[None]
        psi = rw.period(psi)
    if time.time() - T0 > 400:
        np.savez(ck, psi=psi, acc=acc, t=NPER); print("checkpoint before Ritz", flush=True); sys.exit(3)
    F = acc.reshape(len(OM), -1); del acc
    G = F.conj() @ F.T
    ev, U = np.linalg.eigh(G)
    keep = ev > 1e-11 * ev.max()
    Q = (U[:, keep] / np.sqrt(ev[keep])).T @ F
    UQ = np.array([rw.period(q.reshape(psi0.shape)).ravel() for q in Q])
    lam, Vv = np.linalg.eig(Q.conj() @ UQ.T)
    lev = []
    for i in range(len(lam)):
        v = Vv[:, i] @ Q; nv = np.linalg.norm(v); v /= nv
        res = float(np.linalg.norm(Vv[:, i] @ UQ / nv - lam[i] * v))
        E = float(-np.angle(lam[i]) / (th * w.T))
        wgt = float(abs(np.vdot(v, psi0.ravel())) ** 2)
        if 0.80 <= E <= 0.97 and wgt >= 0.005:
            vv = v.reshape(4, -1)
            S = [float(np.real(np.einsum("an,ab,bn->", vv.conj(), Sm, vv))) for Sm in SPIN]
            lev.append(dict(E=E, weight=wgt, resid=res, S=S))
    lev.sort(key=lambda d: -d["weight"])
    pair = None
    for i in range(len(lev)):
        for j in range(i + 1, len(lev)):
            a, b = np.array(lev[i]["S"]), np.array(lev[j]["S"])
            if np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-30) <= -0.9:
                pair = (i, j); break
        if pair:
            break
    out = dict(levels=lev, periods=NPER, box=CONF[th][0])
    if pair:
        out["doublet"] = [lev[pair[0]], lev[pair[1]]]
        out["split"] = abs(lev[pair[0]]["E"] - lev[pair[1]]["E"])
    if os.path.exists(ck):
        os.remove(ck)
    return out


for th in (0.7, 0.5, 0.35, 0.25):
    key = str(th)
    if key in R:
        continue
    R[key] = run(th)
    print(th, "split", R[key].get("split"), [(round(l["E"], 5), round(l["weight"], 3)) for l in R[key]["levels"][:4]], flush=True)
    json.dump(R, open(OUT, "w"), indent=1)
    if time.time() - T0 > 300:
        print("rerun"); sys.exit(3)
ths = [float(k) for k in R if "split" in R[k]]
if len(ths) >= 3:
    sp = [R[str(t)]["split"] for t in ths]
    R["fit"] = dict(exponent=float(np.polyfit(np.log(ths), np.log(sp), 1)[0]), thetas=ths, splits=sp)
    json.dump(R, open(OUT, "w"), indent=1)
    print("fit", R["fit"])
print("ALL DONE")
