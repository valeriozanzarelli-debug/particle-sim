"""senzadoppi_r1.py - exams E1 (theorem W3 = 0, any number of components, any lattice), E2 (census of nodes and
valleys of F0, B1, B2, PA), E3 (search over all 2-tick axis orders for L and R).  Output lab/out/senzadoppi_r1.json.
Run: OMP_NUM_THREADS=2 timeout 590 python lab/senzadoppi_r1.py   (checkpointed; rerun until exit 0)"""
import os
import sys
import json
import time
import itertools
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import senzadoppi_core as s

OUT = os.path.join(os.path.dirname(__file__), "out", "senzadoppi_r1.json")
CK = os.path.join(os.environ.get("SD_SCRATCH", "/tmp/claude-0/-home-claude-particle-sim/9de435bc-26c2-5d88-b70d-07efb8033e5b/scratchpad"),
                  "senzadoppi_r1_ck.json")
T0 = time.time()
BUDGET = 540


def load():
    if os.path.exists(CK):
        return json.load(open(CK))
    return {}


def save(d):
    json.dump(d, open(CK, "w"), indent=1)


# ------------------------------------------------------------------ E1: W3 of products
def w3(Ufun, n=36):
    """(1/24 pi^2) int eps tr(U^-1 dU)^3, central differences on an n^3 grid (periodic)"""
    g = (np.arange(n) + 0.5) / n * 2 * np.pi - np.pi
    K = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T
    h = 2 * np.pi / n
    U = Ufun(K)
    Ui = np.conj(np.transpose(U, (0, 2, 1)))
    A = []
    e = 1e-5
    for i in range(3):
        dK = np.zeros(3)
        dK[i] = e
        D = (Ufun(K + dK) - Ufun(K - dK)) / (2 * e)
        A.append(Ui @ D)
    tot = 0
    for (i, j, k), sg in (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1),
                          ((1, 0, 2), -1)):
        tot = tot + sg * np.trace(A[i] @ A[j] @ A[k], axis1=1, axis2=2)
    return float(np.real(np.sum(tot)) * h**3 / (24 * np.pi**2)), float(np.imag(np.sum(tot)) * h**3 / (24 * np.pi**2))


def rand_unitary(rng, d):
    Z = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def random_rule(rng, d, nfac):
    """list of factor functions K -> (N,d,d), each depending on ONE combination k.v"""
    facs = []
    dirs = [np.array(v, float) for v in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1], [1, -1, 1], [1, 1, -1],
                                         [-1, 1, 1], [1, 1, 0], [0, 1, -1], [2, 0, 0])]
    for _ in range(nfac):
        kind = rng.integers(3)
        v = dirs[rng.integers(len(dirs))]
        if kind == 0:     # shift: eigen +1 part of a random Hermitian with spectrum +-1, 0 moves by lam*v
            V = rand_unitary(rng, d)
            lam = rng.choice([-1, 0, 1], size=d)
            facs.append(("shift", v, V, lam))
        elif kind == 1:   # coin
            facs.append(("coin", None, rand_unitary(rng, d), None))
        else:             # link: S_v^dag C S_v with random projector
            V = rand_unitary(rng, d)
            q = rng.integers(0, 2, size=d)
            facs.append(("link", v, V, q, rand_unitary(rng, d)))

    def U(K):
        N = K.shape[0]
        M = np.broadcast_to(np.eye(d, dtype=complex), (N, d, d)).copy()
        for f in facs:
            if f[0] == "shift":
                ph = np.exp(-1j * np.outer(K @ f[1], f[3]))
                F = np.einsum("ab,nb,cb->nac", f[2], ph, f[2].conj())
            elif f[0] == "coin":
                F = np.broadcast_to(f[2], (N, d, d))
            else:
                ph = np.exp(-1j * np.outer(K @ f[1], f[3]))
                S = np.einsum("ab,nb,cb->nac", f[2], ph, f[2].conj())
                F = np.conj(np.transpose(S, (0, 2, 1))) @ f[4][None] @ S
            M = F @ M
        return M
    return U


def control(M):
    def U(K):
        a = np.stack([M - np.sum(np.cos(K), 1), np.sin(K[:, 0]), np.sin(K[:, 1]), np.sin(K[:, 2])], 1)
        a /= np.linalg.norm(a, axis=1)[:, None]
        return (a[:, 0, None, None] * s.I2[None] + 1j * np.einsum("ni,ijk->njk", a[:, 1:], s.PAULI))
    return U


# ------------------------------------------------------------------ rules
F0 = s.Rule([((0, 1, 2), (0, 1, 2)), ((2, 1, 0), (2, 1, 0))], name="F0")
B1 = s.Rule([((0, 1, 2), (0, 1, 2))], name="B1")
B2 = s.Rule([((0, 1, 2), (2, 1, 0))], name="B2")


class OneCoin(s.Rule):
    """translations of a multi-tick period, ONE coin of angle T*theta per period (+ optional Wilson around it)"""

    def __init__(self, ticks, h=0.0, vecs=None, name=""):
        super().__init__(ticks, name=name)
        self.hh = h
        self.vecs = vecs if vecs is not None else [2 * e for e in np.eye(3)] + [-2 * e for e in np.eye(3)]

    def bloch(self, K, theta):
        K = np.atleast_2d(np.asarray(K, float))
        N = K.shape[0]
        W = np.broadcast_to(np.eye(4, dtype=complex), (N, 4, 4)).copy()
        for n in range(self.T):
            W = self.tick_trans(K, n) @ W
        C = s.coin(self.T * theta)[None]
        if self.hh:
            W = s.wilson_links(K, self.hh / 2, vecs=self.vecs[::-1]) @ C @ s.wilson_links(K, self.hh / 2, vecs=self.vecs) @ W
        else:
            W = C @ W
        return W


class WilsonEach(s.Rule):
    """rule with a coin each tick and Wilson split around each coin on given link vectors (mirrored on odd ticks)"""

    def __init__(self, ticks, h, vecs, name=""):
        super().__init__(ticks, name=name)
        self.hh, self.vecs = h, vecs

    def bloch(self, K, theta):
        K = np.atleast_2d(np.asarray(K, float))
        N = K.shape[0]
        W = np.broadcast_to(np.eye(4, dtype=complex), (N, 4, 4)).copy()
        C = s.coin(theta)[None]
        for n in range(self.T):
            W = self.tick_trans(K, n) @ W
            v = self.vecs if n % 2 == 0 else self.vecs[::-1]
            W = s.wilson_links(K, self.hh / 2, vecs=v[::-1]) @ C @ s.wilson_links(K, self.hh / 2, vecs=v) @ W
        return W


PA = OneCoin([((0, 1, 2), (0, 1, 2)), ((2, 1, 0), (2, 1, 0))], name="PA")
VEC_A = [2 * e for e in np.eye(3)] + [-2 * e for e in np.eye(3)]
VEC_B = [np.array(d, float) for d in itertools.product((1, -1), repeat=3)]
PAW = OneCoin([((0, 1, 2), (0, 1, 2)), ((2, 1, 0), (2, 1, 0))], h=0.3, vecs=VEC_A, name="PAW")
B2W = WilsonEach([((0, 1, 2), (2, 1, 0))], 0.3, VEC_B, name="B2W")
RULES = {"F0": F0, "B1": B1, "B2": B2, "PA": PA, "PAW": PAW, "B2W": B2W}


def canon(pr):
    """canonical representative of a 2-tick order pattern (R0, L0, R1, L1) under global axis relabeling and
    cyclic exchange of the two ticks (both give the same spectrum)"""
    best = None
    for sg in s.PERMS:
        rel = [tuple(sg[j] for j in p) for p in pr]
        for cand in (tuple(rel), (rel[2], rel[3], rel[0], rel[1])):
            if best is None or cand < best:
                best = cand
    return best


def main():
    d = load()
    rng = np.random.default_rng(20261013)
    # ---------------- E1
    if "E1" not in d:
        res = []
        for i in range(24):
            dim = 2 if i < 12 else 4
            U = random_rule(rng, dim, int(rng.integers(3, 8)))
            r, im = w3(U)
            res.append({"dim": dim, "W3": r, "imag": im})
        ctrl = {"M2": w3(control(2.0))[0], "M4": w3(control(4.0))[0]}
        # the frozen rule and the candidates (R part 2x2 and full 4x4 at theta = 0,3)
        rules_w3 = {}
        for nm, R in RULES.items():
            rules_w3[nm + "_R"] = w3(lambda K, R=R: R.bloch(K, 0.0)[:, 2:, 2:])[0] if getattr(R, "hh", 0) == 0 else \
                w3(lambda K, R=R: R.bloch(K, 0.0))[0]
            rules_w3[nm + "_4x4"] = w3(lambda K, R=R: R.bloch(K, 0.3))[0]
        ok = all(abs(x["W3"]) < 0.01 for x in res) and abs(abs(ctrl["M2"]) - 1) < 0.01 and abs(ctrl["M4"]) < 0.01
        d["E1"] = {"random": res, "control": ctrl, "rules": rules_w3, "max_abs_random": max(abs(x["W3"]) for x in res),
                   "PASS": bool(ok)}
        save(d)
        print("E1", d["E1"]["max_abs_random"], ctrl, rules_w3, ok, flush=True)
    # ---------------- E2
    if "E2" not in d:
        cen = {}
        for nm in ("F0", "B1", "B2", "PA"):
            R = RULES[nm]
            nodes = s.weyl_nodes(R if nm != "PA" else F0, n_starts=250)
            nc = len(s.copy_shifts(R))
            n0 = sum(1 for x in nodes if x[1] == 1)
            npi = len(nodes) - n0
            chir0 = sum(x[2] for x in nodes if x[1] == 1)
            chirpi = sum(x[2] for x in nodes if x[1] == -1)
            vs = s.valley_scan(R, 0.1, n=32)
            cen[nm] = {"copies": nc, "nodes0_per_copy": n0 / nc, "nodespi_per_copy": npi / nc,
                       "chirality_sum_0": chir0, "chirality_sum_pi": chirpi, "valleys_theta0.1": vs,
                       "nodes": [[list(np.round(x[0] / np.pi, 4)), x[1], x[2]] for x in nodes]}
            print("E2", nm, nc, n0 / nc, npi / nc, chir0, chirpi, vs["E0"], vs["Epi"], flush=True)
        pred = {"F0": (4, 4, 0.1, 0.0), "B1": (2, 2, 0.0, None), "B2": (2, 2, 0.1, 0.1), "PA": (4, 4, 0.1, 0.1)}
        passes = {}
        for nm, (a, b, e0, epi) in pred.items():
            c = cen[nm]
            ok = c["nodes0_per_copy"] == a and c["nodespi_per_copy"] == b and abs(c["valleys_theta0.1"]["E0"] - e0) < 1e-6
            if epi is not None:
                ok = ok and abs(c["valleys_theta0.1"]["Epi"] - epi) < 1e-6
            passes[nm] = bool(ok)
        cen["PASS"] = passes
        d["E2"] = cen
        save(d)
    # ---------------- E3 search
    if "E3_done" not in d:
        pats = sorted({canon(p) for p in itertools.product(s.PERMS, repeat=4)})
        d.setdefault("E3", {})
        d["E3_n_classes"] = len(pats)
        for p in pats:
            key = "|".join(s.PNAME[q] for q in p)
            if key in d["E3"]:
                continue
            if time.time() - T0 > BUDGET:
                save(d)
                print("checkpoint", len(d["E3"]), "/", len(pats), flush=True)
                sys.exit(3)
            R = s.Rule([(p[0], p[1]), (p[2], p[3])], name=key)
            vs = s.valley_scan(R, 0.1, n=20)
            iso = s.isotropy(R, ndir=120)
            d["E3"][key] = {"E0": vs["E0"], "Epi": vs["Epi"], "copies": vs["n_copies"], "iso24": iso[0],
                            "vmax24": iso[3]}
        d["E3_done"] = True
        save(d)
    # summary of E3
    E3 = d["E3"]
    best0 = max(E3.values(), key=lambda v: v["E0"])["E0"]
    good = {k: v for k, v in E3.items() if v["Epi"] >= 0.1 - 1e-6 and v["iso24"] <= 0.02}
    d["E3_summary"] = {"n_classes": len(E3), "max_E0": best0,
                       "n_E0_gt_theta": sum(1 for v in E3.values() if v["E0"] > 0.1 + 1e-6),
                       "n_Epi_ge_theta": sum(1 for v in E3.values() if v["Epi"] >= 0.1 - 1e-6),
                       "n_Epi_ge_theta_and_iso_le_2pc": len(good), "those": good,
                       "n_E0_zero": sum(1 for v in E3.values() if v["E0"] < 1e-6),
                       "best_iso": min(E3.items(), key=lambda kv: kv[1]["iso24"]),
                       "PASS": bool(best0 <= 0.1 + 1e-6 and len(good) == 0)}
    print("E3", json.dumps({k: v for k, v in d["E3_summary"].items() if k != "those"}), flush=True)
    print("E3 good", good)
    json.dump(d, open(OUT, "w"), indent=1)
    print("time", time.time() - T0)


if __name__ == "__main__":
    main()
