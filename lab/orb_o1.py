"""orb_o1.py - O1: holonomy from walking balls on orbifolds. Writes lab/out/orb_o1.json.
(a) 4-D zig-zag split-step walk on n^4 torus (spin matrices alpha_1..4, flip exp(-i theta beta)).
    Orbifold operators: K3 = T4/Z2 (x -> -x), control Z2 on one complex coordinate (x1,x2 -> -), control reflection.
(b) hexagonal Wilson-Dirac walker exp(-i H dt): T2/Z3, T6/Z3 = (T2)^3 with simultaneous Z3, control Z3 on one T2.
For each: S solved from the walk (linear algebra), compared with the ball's own rotation lifts, holonomy around every
fixed point, invariant spinors."""
import json, time
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from orb_core import *

t0 = time.time()
rng = np.random.default_rng(20261009)
R = {}

# ------------------------------------------------------------------ (a) 4-D zig-zag walk
AL = [kron(SX, SX), kron(SX, SY), kron(SX, SZ), kron(SY, I2)]
BETA = kron(SZ, I2)
EIG = []
for a in AL:
    w, v = np.linalg.eigh(a); EIG.append((np.round(w).astype(int), v))


def walk4(n, theta):
    def W(psi):
        p = psi.reshape(4, n, n, n, n)
        for j in range(4):
            lam, V = EIG[j]
            c = np.tensordot(V.conj().T, p, axes=(1, 0))
            for q in range(4):
                c[q] = np.roll(c[q], lam[q], axis=j)
            p = np.tensordot(V, c, axes=(1, 0))
        U = sla.expm(-1j * theta * BETA)
        p = np.tensordot(U, p, axes=(1, 0))
        return p.reshape(4, -1)
    return W


def perm4(n, signs, shift=(0, 0, 0, 0)):
    """(P psi)(x) = psi(g^-1 x), g x = signs*x + shift"""
    def P(psi):
        p = psi.reshape(4, n, n, n, n)
        idx = np.indices((n,) * 4)
        src = [((signs[d] * (idx[d] - shift[d])) % n) for d in range(4)]   # g^-1 x
        return p[:, src[0], src[1], src[2], src[3]].reshape(4, -1)
    return P


def chir_inv(S):
    w, v = np.linalg.eig(S)
    g5 = AL[0] @ AL[1] @ AL[2] @ AL[3]
    return [round(float(np.real(v[:, q].conj() @ g5 @ v[:, q] / (v[:, q].conj() @ v[:, q]))), 4)
            for q in range(4) if abs(w[q] - 1) < 1e-8]


def rot_lift(planes):
    """ball's own rotation exp(phi/2 alpha_a alpha_b) for each (a, b, phi)"""
    S = np.eye(4, dtype=complex)
    for a, b, phi in planes:
        S = S @ (np.cos(phi / 2) * np.eye(4) + np.sin(phi / 2) * AL[a] @ AL[b])
    return S


cases4 = {
    "K3: T4/Z2 x->-x (angles pi,pi)": ((-1, -1, -1, -1), [(0, 1, np.pi), (2, 3, np.pi)]),
    "control: Z2 on one complex coord (pi,0)": ((-1, -1, 1, 1), [(0, 1, np.pi)]),
    "control: reflection x1->-x1": ((-1, 1, 1, 1), None),
}
for n in (4, 6):
    W = walk4(n, theta=0.3)
    for name, (signs, planes) in cases4.items():
        basis = solve_spinor_operator(W, perm4(n, signs), 4, n ** 4, rng)
        d = {"n": n, "dim_solution_space": len(basis)}
        if planes is not None:
            lifts = {"+": rot_lift(planes), "-": -rot_lift(planes)}
            d["lifts"] = {}
            for k, S in lifts.items():
                ok = in_span(S, basis)
                # holonomy around every fixed point f: O_f psi(x) = S psi(2f - x)
                fixed = [f for f in np.ndindex(*(2,) * 4)]
                fps = []
                for f in fixed:
                    sh = tuple((n // 2) * f[q] * 2 if signs[q] < 0 else 0 for q in range(4))
                    # g x = signs x + 2 f (n/2) ; acting coordinates only
                    sh = tuple(((n // 2) * f[q] * (1 - signs[q])) % n for q in range(4))
                    fps.append(sh)
                fps = sorted(set(fps))
                res = []
                psi = rng.normal(size=(4, n ** 4)) + 1j * rng.normal(size=(4, n ** 4))
                for sh in fps:
                    P = perm4(n, signs, sh)
                    O = lambda x: S @ P(x)
                    res.append(float(np.linalg.norm(W(O(psi)) - O(W(psi))) / np.linalg.norm(psi)))
                ev = np.linalg.eigvals(S)
                d["lifts"][k] = dict(S_commutes_with_walk=bool(ok), n_fixed_point_loops=len(fps),
                                     max_residual_all_fixed_points=max(res), S_order2_is_identity=bool(
                                         np.allclose(S @ S, np.eye(4))),
                                     eigenvalues=[[round(float(z.real), 6), round(float(z.imag), 6)] for z in ev],
                                     invariant_spinors=n_invariant(S),
                                     chirality_of_invariant=chir_inv(S))
            # loophole: extra U(1) phase allowed
            ev = np.linalg.eigvals(lifts["+"])
            vals, cnt = np.unique(np.round(np.angle(ev), 6), return_counts=True)
            d["with_extra_U1_phase_max_invariant"] = int(cnt.max())
        if signs == (-1, 1, 1, 1):
            Wm = walk4(n, theta=0.0)
            d["dim_solution_space_massless"] = len(solve_spinor_operator(Wm, perm4(n, signs), 4, n ** 4, rng))
        R.setdefault(name, []).append(d)
        print(n, name, d.get("dim_solution_space"), {k: v["invariant_spinors"] for k, v in d.get("lifts", {}).items()},
              f"{time.time()-t0:.1f}s", flush=True)

# orbifold dynamics check (K3, n = 6): packet in the invariant sector stays invariant
n = 6; W = walk4(n, 0.3); P = perm4(n, (-1, -1, -1, -1)); S = rot_lift([(0, 1, np.pi), (2, 3, np.pi)])
O = lambda x: S @ P(x)
idx = np.indices((n,) * 4); r2 = sum((((idx[d] - 1 + n // 2) % n) - n // 2) ** 2 for d in range(4))
phi = np.exp(-r2 / 2.0).reshape(-1)[None, :] * np.array([1, 0, 0, 0])[:, None] * np.exp(1j * 0.8 * idx[0].reshape(-1))[None, :]
psi = (phi + O(phi)) / 2; psi /= np.linalg.norm(psi)
leak = []
for t in range(60):
    psi = W(psi); leak.append(float(np.linalg.norm(psi - O(psi))))
R["K3_packet_60_ticks"] = dict(max_leak_out_of_orbifold_sector=max(leak), norm=float(np.linalg.norm(psi)))

# ------------------------------------------------------------------ (b) hexagonal walkers
dt = 0.5


def hex_case(L, ntor, act, angles_lift):
    """ntor = 1 or 3 tori; act = tuple of bools: which tori are rotated."""
    N = L * L
    H1 = hex_dirac(L)
    # gammas for ntor tori (Jordan-Wigner), 2^ntor components
    if ntor == 1:
        H = H1
        D = 2
    else:
        D = 8
        Gx = [kron(SX, I2, I2), kron(SZ, SX, I2), kron(SZ, SZ, SX)]
        Gy = [kron(SY, I2, I2), kron(SZ, SY, I2), kron(SZ, SZ, SY)]
        G7 = kron(SZ, SZ, SZ)
        H = sp.csr_matrix((D * N ** 3, D * N ** 3), dtype=complex)
        Is = sp.eye(N, format="csr")
        for t in range(3):
            for e, nv in zip(HEX_E, HEX_N):
                T = hex_shift(L, e)
                ops = [Is, Is, Is]; ops[t] = (T - T.T) / (2j) * (2 / 3)
                H = H + sp.kron(nv[0] * Gx[t] + nv[1] * Gy[t], sp.kron(sp.kron(ops[0], ops[1]), ops[2]))
                ops = [Is, Is, Is]; ops[t] = sp.eye(N) - (T + T.T) / 2
                H = H + sp.kron(G7, sp.kron(sp.kron(ops[0], ops[1]), ops[2]))
        H = H.tocsr()
    nv_ = H.shape[0] // D

    def W(psi):  # exp(-i H dt) to 12th order Taylor (unitary to 1e-12 for ||H dt|| < 1.5)
        v = psi.reshape(-1); out = v.copy(); term = v.copy()
        for k in range(1, 25):
            term = (-1j * dt / k) * (H @ term); out = out + term
        return out.reshape(D, nv_)

    def Pmaker(centers):
        Ps = [hex_rot_perm(L, c) if a else sp.eye(N, format="csr") for a, c in zip(act, centers)]
        Pf = Ps[0]
        for q in Ps[1:]:
            Pf = sp.kron(Pf, q, format="csr")
        return lambda psi: (Pf @ psi.reshape(D, nv_).T).T
    fixed1 = hex_fixed_points(L)
    basis = solve_spinor_operator(W, Pmaker([(0, 0)] * ntor), D, nv_, rng, nrand=1)
    out = {"L": L, "tori": ntor, "rotated": act, "fixed_points_T2": fixed1, "dim_solution_space": len(basis), "lifts": {}}
    for k, sgn in (("+", 1), ("-", -1)):
        S = np.eye(1, dtype=complex)
        for t in range(ntor):
            phi = angles_lift[t]
            S = np.kron(S, sla.expm(-1j * phi / 2 * SZ))
        S = sgn * S
        ok = in_span(S, basis)
        # holonomy around each fixed point (all combinations for the rotated tori)
        import itertools
        combos = list(itertools.product(*[fixed1 if a else [(0, 0)] for a in act]))
        if len(combos) > 9:
            combos = combos[:: max(1, len(combos) // 9)]
        psi = rng.normal(size=(D, nv_)) + 1j * rng.normal(size=(D, nv_))
        res = []
        for cs in combos:
            P = Pmaker(list(cs)); O = lambda x: S @ P(x)
            res.append(float(np.linalg.norm(W(O(psi)) - O(W(psi))) / np.linalg.norm(psi)))
        S3 = np.linalg.matrix_power(S, 3)
        out["lifts"][k] = dict(S_commutes_with_walk=bool(ok), S_cubed_is_identity=bool(np.allclose(S3, np.eye(D))),
                               fixed_point_loops_tested=len(combos), max_residual=max(res),
                               invariant_spinors=n_invariant(S))
    ev = np.linalg.eigvals(S)
    _, cnt = np.unique(np.round(np.angle(ev), 6), return_counts=True)
    out["with_extra_U1_phase_max_invariant"] = int(cnt.max())
    return out


w = 2 * np.pi / 3
hcases = {
    "T2/Z3": (6, 1, (True,), (w,)),
    "T6/Z3 (w,w,w)": (3, 3, (True, True, True), (w, w, -2 * w)),   # spin lift with angle sum 0 (S^3 = 1)
    "control: T6, Z3 on one T2 (w,0,0)": (3, 3, (True, False, False), (w, 0, 0)),
}
for name, (L, nt, act, ang) in hcases.items():
    d = hex_case(L, nt, act, ang)
    R[name] = d
    print(name, d["dim_solution_space"], {k: (v["S_commutes_with_walk"], v["invariant_spinors"], v["S_cubed_is_identity"])
                                          for k, v in d["lifts"].items()}, f"{time.time()-t0:.1f}s", flush=True)
# T2/Z3 with the other angle choice (-4pi/3 = same rotation) for completeness
d = hex_case(6, 1, (True,), (-2 * w,)); R["T2/Z3 (lift via -4pi/3)"] = d
print("T2/Z3 alt", {k: (v["invariant_spinors"], v["S_cubed_is_identity"]) for k, v in d["lifts"].items()})

# ------------------------------------------------------------------ grading
def inv(name, n=None):
    v = R[name] if n is None else R[name][-1]
    return max(l["invariant_spinors"] for l in v["lifts"].values() if l["S_commutes_with_walk"])
G = {
    "K3 invariant spinors (walk commutes, all 16 loops)": inv("K3: T4/Z2 x->-x (angles pi,pi)", 1),
    "K3 residual max": max(l["max_residual_all_fixed_points"] for d in R["K3: T4/Z2 x->-x (angles pi,pi)"] for l in d["lifts"].values()),
    "T6/Z3 invariant": inv("T6/Z3 (w,w,w)"),
    "T2/Z3 invariant": max(inv("T2/Z3"), inv("T2/Z3 (lift via -4pi/3)")),
    "control one-coord Z2 invariant": inv("control: Z2 on one complex coord (pi,0)", 1),
    "control one-coord Z3 invariant": inv("control: T6, Z3 on one T2 (w,0,0)"),
    "control reflection solution dim (massive)": [d["dim_solution_space"] for d in R["control: reflection x1->-x1"]],
}
G["PASS"] = bool(G["K3 invariant spinors (walk commutes, all 16 loops)"] >= 1 and G["K3 residual max"] < 1e-10 and
                 G["T6/Z3 invariant"] >= 1 and G["T2/Z3 invariant"] == 0 and G["control one-coord Z2 invariant"] == 0
                 and G["control one-coord Z3 invariant"] == 0 and max(G["control reflection solution dim (massive)"]) == 0)
R["grading"] = G
R["seconds"] = time.time() - t0
print(json.dumps(G, indent=1))
json.dump(R, open(f"{OUT}/orb_o1.json", "w"), indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
