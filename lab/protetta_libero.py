"""protetta_libero.py - free level of the composite photon (exams L1-L4 of lab/out/protetta_prereg.json).
The rule's light walker vs a pair of R-Weyl walkers of the rule at k/2 +- q.  Output lab/out/protetta_libero.json.
Run: python3 lab/protetta_libero.py   (seconds)"""
import json, os
import numpy as np
from scipy.linalg import null_space
from zz_core import PAULI, SPIN1, T_axis, quasi_energies, SEED

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
ORDERS = {"R": [(0, 1, 2), (2, 1, 0)],
          "C": [(0, 1, 2), (2, 1, 0), (1, 2, 0), (0, 2, 1), (2, 0, 1), (1, 0, 2)]}


def W(k, rule, spin):
    """period operator of the massless walker (R-Weyl spin 1/2: exp(-i k_j sigma_j); light spin 1: exp(-i k_j S_j))"""
    M = np.eye(2 if spin == 0.5 else 3, dtype=complex)
    for order in ORDERS[rule]:
        for j in order:
            M = T_axis(k[j], j, +1, spin) @ M
    return M


# triplet = symmetric subspace of C2 x C2, generators (sigma x 1 + 1 x sigma)/2 ; intertwiner V: V S_k = T_k V
I2 = np.eye(2)
P = null_space(np.array([[1, 0, 0, 0], [0, 1, -1, 0], [0, -1, 1, 0], [0, 0, 0, 1]]) - np.eye(4))   # sym. projector range
SW = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
Ps = null_space(SW - np.eye(4))                       # 4x3 isometry onto the symmetric subspace
Tg = [Ps.conj().T @ (np.kron(s, I2) + np.kron(I2, s)) / 2 @ Ps for s in PAULI]
# solve V S_k - T_k V = 0 (linear in V)
A = np.vstack([np.kron(np.eye(3), Sk.T) - np.kron(Tk, np.eye(3)) for Sk, Tk in zip(SPIN1, Tg)])
v = null_space(A)
V = v[:, 0].reshape(3, 3)
V = V / np.sqrt(np.trace(V.conj().T @ V).real / 3)
res = {"intertwiner_unitarity": float(np.max(np.abs(V.conj().T @ V - np.eye(3)))), "null_dim": int(v.shape[1])}


def phases(M):
    return np.sort(np.angle(np.linalg.eigvals(M)))


rng = np.random.default_rng(SEED)
L1 = {}
for rule in ("R", "C"):
    es, eo = 0.0, 0.0
    for _ in range(2000):
        k = rng.uniform(-np.pi, np.pi, 3)
        Wl = W(k, rule, 1)
        Wh = W(k / 2, rule, 0.5)
        Wc = Ps.conj().T @ np.kron(Wh, Wh) @ Ps          # triplet block (the singlet is invariant: exp(-i k S) on it = 1)
        d = np.abs(np.exp(1j * phases(Wl)) - np.exp(1j * phases(Wc)))
        es = max(es, float(np.max(d)))
        eo = max(eo, float(np.max(np.abs(V @ Wl @ V.conj().T - Wc))))
    L1[rule] = {"max_spectrum_diff": es, "max_operator_diff": eo, "PASS": bool(es < 1e-12 and eo < 1e-12)}
res["L1"] = L1


def Om_plus(k, rule, spin):
    """largest positive period phase Omega(k) (e^{-i Omega})"""
    return np.max(-np.angle(np.linalg.eigvals(W(k, rule, spin))))


def vgroup(f, k, n, h=1e-6):
    return (f(k + h * n) - f(k - h * n)) / (2 * h)


L2 = {}
dirs = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1), "xy": (1, 1, 0), "xyz": (1, 1, 1)}
for rule in ("R", "C"):
    nt = len(ORDERS[rule])
    worst = 0.0
    rows = {}
    for name, d in dirs.items():
        n = np.array(d, float) / np.linalg.norm(d)
        for ka in (1e-3, 0.05, 0.2, 0.5):
            k = ka * n
            vl = vgroup(lambda kk: Om_plus(kk, rule, 1), k, n) / nt
            vw = vgroup(lambda kk: Om_plus(kk, rule, 0.5), k / 2, n) / nt
            worst = max(worst, abs(vl - vw))
            rows["%s|%g" % (name, ka)] = {"v_light": vl, "v_weyl_half": vw}
    L2[rule] = {"max_diff": worst, "PASS": bool(worst < 1e-9), "rows": rows}
res["L2"] = L2

# L3 non-collinear pair: k along x, q along z (and along y)
L3 = {}
ka = 0.05
for qdir in (2, 1):
    for r in (0.01, 0.02, 0.05):
        q = np.zeros(3); q[qdir] = r * ka
        n = np.array([1.0, 0, 0])
        f = lambda kk: (Om_plus(kk / 2 + q, "R", 0.5) + Om_plus(kk / 2 - q, "R", 0.5)) / 2
        v = vgroup(f, ka * n, n, h=1e-7)
        L3["q%s|%g" % ("zyx"[2 - qdir], r)] = {"v": v, "ratio": (1 - v) / (2 * r**2)}
res["L3"] = {"rows": L3, "PASS": bool(all(abs(x["ratio"] - 1) < 0.05 for x in L3.values()))}

# L4 flipping (massive) constituents: rest quasi-energy of the q=0 pair
L4 = {}
for th in (0.05, 0.0125):
    om, _ = quasi_energies(np.zeros(3), th)
    E = float(np.max(om))                      # positive band at rest, per tick
    L4["%g" % th] = {"pair_rest_over_2theta": 2 * E / (2 * th)}
res["L4"] = {"rows": L4, "PASS": bool(all(abs(x["pair_rest_over_2theta"] - 1) < 1e-3 for x in L4.values()))}
json.dump(res, open(os.path.join(OUT, "protetta_libero.json"), "w"), indent=1)
print(json.dumps({k: (v if k != "L2" else {r: {"max_diff": v[r]["max_diff"], "PASS": v[r]["PASS"]} for r in v})
                  for k, v in res.items()}, indent=1))
