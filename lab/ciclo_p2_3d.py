"""ciclo_p2_3d.py - phase 2 (3-D, k-space): the loop on the free 3-D zig-zag from the existing chain-A level B = 16,
with blocking factors b = 2, 3, phi (generator transported as b G(kappa/b)). Also exam A2 (automatic descriptor).
Run: python lab/ciclo_p2_3d.py <2|3|phi|A2>     -> lab/out/ciclo_p2_3d_<scheme>.json"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, pickle, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zz_core as zz
import agg_zz_core as ac
ac.WK = 1

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A0 = 1 / zz.ALPHA
TH = A0 / 16 / 256
PHI = (1 + 5 ** 0.5) / 2
KS = zz.ALPHA                      # 1s momentum scale (units of m)
E1S = zz.ALPHA ** 2 / 2
NK = 6000
SEED = 20261013


class Stencil:
    def __init__(self, g, offs):
        self.g, self.offs = g, offs.astype(float)

    def __call__(self, kv):
        return np.real(np.exp(1j * kv @ self.offs.T) @ self.g)


def fine_G(kv, B):
    """fine walk's exact positive-band kinetic phase over B ticks, at coarse wave vectors kv (cells of size B)"""
    k = kv / B
    z = np.zeros((1, 3))
    Wm = zz.bloch_batch(np.concatenate([k, z]), TH)
    lam = np.linalg.eigvals(Wm)
    om = -np.angle(lam) / 2
    op = np.sort(om, axis=1)[:, -1]
    return (op[:-1] - op[-1]) * B


def sample(rng, n, B):
    """75 % like the 1s (cells of B fine cells: k_s = alpha m h), 25 % uniform in the cube"""
    h = B * TH
    kg = ac.sample_state_k(rng, int(0.75 * n), 3, KS * h)
    kg = kg[np.all(np.abs(kg) < np.pi, axis=1)]
    ku = np.pi * (2 * rng.random((n - int(0.75 * n), 3)) - 1)
    return kg, ku


def c2_axis(G, offs):
    return np.array([-0.5 * np.sum(G * offs[:, a] ** 2) for a in range(3)])


def learn(prev, b, B, w, rng):
    """least squares for a real even stencil with G(0) = 0 and the axis c2 transported (b c2_prev / b^2 ... exactly as
    b G_prev(kappa/b) gives at small kappa)"""
    offs = ac.offsets(w, 3)
    # parity: g_j = g_{-j}: parametrize by representatives
    key = {}
    rep = []
    for i, o in enumerate(offs):
        t = tuple(o); nt = tuple(-o)
        if nt in key:
            key[t] = key[nt]
        else:
            key[t] = len(rep); rep.append(i)
    P = len(rep)
    M = np.zeros((len(offs), P))
    for i, o in enumerate(offs):
        M[i, key[tuple(o)]] = 1.0
    kg, ku = sample(rng, NK, B)
    kap = np.concatenate([kg, ku])
    tgt = b * prev(kap / b)
    sw = np.ones(len(kap)); sw[len(kg):] = 0.03
    sw = sw / (np.abs(tgt) + 1e-3 * np.median(np.abs(tgt[:len(kg)])))
    A = np.cos(kap @ offs.T.astype(float)) @ M
    # constraints: G(0) = 0; c2 per axis = transported value
    c2t = c2_axis(prev.g, prev.offs) / b          # b * c2_prev * (1/b)^2
    Ceq = [np.ones(len(offs)) @ M]
    for a in range(3):
        Ceq.append((-0.5 * offs[:, a] ** 2) @ M)
    Ceq = np.array(Ceq); deq = np.concatenate([[0.0], c2t])
    Aw = A * sw[:, None]
    K = np.block([[Aw.T @ Aw, Ceq.T], [Ceq, np.zeros((4, 4))]])
    rhs = np.concatenate([Aw.T @ (tgt * sw), deq])
    x = np.linalg.lstsq(K, rhs, rcond=1e-14)[0][:P]
    return Stencil(M @ x, offs)


def check(prev, cand, b, B, rng):
    kg, ku = sample(rng, NK, B)
    tp = b * prev(kg / b); tc = cand(kg)
    rel = float(np.mean(np.abs(tc - tp)) / np.mean(np.abs(tp)))
    up = b * prev(ku / b); uc = cand(ku)
    bz = float(np.percentile(np.abs(uc - up) / np.maximum(np.abs(up), 1e-12), 90))
    neg = float(max(0.0, -uc.min()) / np.mean(np.abs(tp)))
    return dict(energy_rel=rel, bz_rel90=bz, neg=neg, ok=rel <= 3e-4 and bz <= 0.2 and neg <= 0.01)


def accumulated(st, B, rng):
    kg, _ = sample(rng, 40000, B)
    Gf = fine_G(kg, B)
    Gc = st(kg)
    tau = B * TH
    dE = float(np.mean(Gc - Gf) / tau)          # first-order 1s shift (energy units m c^2)
    return dict(acc_rel=float(np.mean(np.abs(Gc - Gf)) / np.mean(np.abs(Gf))), dE1s_rel=dE / E1S)


def run(scheme):
    b = {"2": 2.0, "3": 3.0, "phi": PHI}[scheme]
    st0 = pickle.load(open(os.path.join(OUT, "agg_zz_3d_A.pkl"), "rb"))["levels"][0]
    g0 = np.real(0.5 * (st0["K"][0, 0] + st0["K"][1, 1]))
    prev = Stencil(g0, st0["offs"])
    B = float(st0["B"])
    rng = np.random.default_rng(SEED)
    levels = [dict(B=B, w=2, taps=125, **accumulated(prev, B, rng))]
    t0 = time.time()
    while B < 600:
        Bn = B * b
        got = None
        tried = []
        for w in range(1, 6):
            cand = learn(prev, b, Bn, w, rng)
            ck = check(prev, cand, b, Bn, rng)
            tried.append(dict(w=w, **ck))
            print(f"  [{scheme}] B={Bn:8.2f} w={w}: energy {ck['energy_rel']:.2e} bz {ck['bz_rel90']:.2e} neg {ck['neg']:.1e} {'OK' if ck['ok'] else ''}", flush=True)
            if ck["ok"]:
                got = (w, cand)
                break
        if got is None:
            levels.append(dict(B=Bn, flag="no admissible w", tried=tried))
            break
        w, cand = got
        acc = accumulated(cand, Bn, rng)
        taps = int(np.count_nonzero(np.abs(cand.g) > 0))
        levels.append(dict(B=Bn, w=w, taps=taps, tried=tried, **acc))
        print(f"  [{scheme}] level B={Bn:.2f} (cell a0/{A0 / (Bn * TH):.2f}) w={w}: acc {acc['acc_rel']:.2e} dE1s {acc['dE1s_rel']:+.2e}", flush=True)
        prev, B = cand, Bn
    res = dict(scheme=scheme, b=b, levels=levels, seconds=time.time() - t0)
    json.dump(res, open(os.path.join(OUT, f"ciclo_p2_3d_{scheme}.json"), "w"), indent=1, default=float)


def exam_A2():
    """automatic descriptor: average the positive-band projector of the fine walk over the 1s momentum density"""
    rng = np.random.default_rng(SEED)
    kv = ac.sample_state_k(rng, 20000, 3, KS * TH)       # fine cells
    Wm = zz.bloch_batch(kv, TH)
    lam, V = np.linalg.eig(Wm)
    om = -np.angle(lam) / 2
    o = np.argsort(-om, axis=1)[:, :2]
    U = np.take_along_axis(V, o[:, None, :], axis=2)                 # (N, 4, 2)
    Uq = np.linalg.qr(U)[0]
    R = np.einsum("nab,ncb->ac", Uq, Uq.conj()) / len(kv)
    e, E = np.linalg.eigh(R)
    top = E[:, np.argsort(e)[::-1][:2]]
    Lg = np.zeros((4, 2), complex)                         # large component (L + R)/sqrt 2 (x) spin
    Lg[0, 0] = Lg[2, 0] = 1 / np.sqrt(2); Lg[1, 1] = Lg[3, 1] = 1 / np.sqrt(2)
    ov = float(np.real(np.trace(top.conj().T @ Lg @ Lg.conj().T @ top)) / 2)
    res = dict(eigs=np.sort(e)[::-1].tolist(), overlap=ov)
    json.dump(res, open(os.path.join(OUT, "ciclo_p2_3d_A2.json"), "w"), indent=1)
    print(res)


if __name__ == "__main__":
    if sys.argv[1] == "A2":
        exam_A2()
    else:
        run(sys.argv[1])
