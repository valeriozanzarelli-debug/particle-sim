"""ciclo_p3_dev.py - phase 3, DEVELOPMENT measurements (declared, not exams): can the cheap level be glued to the base?
(1) block variables of phase 2 glued by Galerkin to a refined region -> spurious states?
(2) smoothed-aggregation prolongator (block variables smoothed s times by the level's own rule, Loewdin): its Galerkin
rule vs the LEARNED phase-2 rule and the bound states vs the base.   -> lab/out/ciclo_p3_dev.json"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ciclo_p3_core as p

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def lowdin(X):
    S = X.conj().T @ X
    e, U = np.linalg.eigh(S)
    return X @ (U / np.sqrt(e)) @ U.conj().T


def SA(W, H, s):
    D = np.real(np.diag(H))
    rho = np.linalg.eigvalsh(H / np.sqrt(D)[:, None] / np.sqrt(D)[None, :]).max()
    X = W.conj().T.copy()
    for _ in range(s):
        X = X - (4 / 3 / rho) * (H @ X) / D[:, None]
    return lowdin(X)


def disp(row, i, k):
    return np.array([np.real(sum(row[i + j] * np.exp(1j * kk * j) for j in range(-8, 9))) for kk in k])


def main():
    hy = p.hierarchy(2)
    N = p.N0
    Hb = hy["Hb"].toarray(); V0 = hy["V0"]
    eb = np.linalg.eigvalsh(Hb + np.diag(V0))[:3]
    res = dict(base_bound=eb.tolist(), block_glue={}, sa={})
    for lev in (0, 1):
        A = p.Adaptive(hy, lev)
        ref = np.zeros(A.nc, bool); ref[A.nc // 2 - 5:A.nc // 2 + 5] = True
        A2 = p.Adaptive(hy, lev, ref)
        L = hy["levels"][lev]; H = L["H"].toarray(); i = H.shape[0] // 2
        Gal = (L["W"] @ hy["Hb"] @ L["W"].conj().T).toarray()
        res["block_glue"][str(int(L["B"]))] = dict(
            coarse_only=np.linalg.eigvalsh(A.H0.toarray())[:3].tolist(),
            ten_cells_refined=np.linalg.eigvalsh(A2.H0.toarray())[:3].tolist(),
            learned_row=H[i, i - 3:i + 4].real.tolist(), galerkin_blocks_row=Gal[i, i - 3:i + 4].real.tolist())
    W1 = hy["levels"][0]["W"].toarray(); W2 = hy["levels"][1]["W"].toarray(); W2r = W2 @ W1.conj().T
    k = np.linspace(0, np.pi / 6, 5)[1:]
    for s in (1, 2, 3):
        P1 = SA(W1, Hb, s); G1 = P1.conj().T @ Hb @ P1
        P = P1 @ SA(W2r, G1, s); G2 = P.conj().T @ Hb @ P
        for lev, (PP, GG) in enumerate(((P1, G1), (P, G2))):
            H = hy["levels"][lev]["H"].toarray(); n = H.shape[0]; i = n // 2; B = N // n
            e = np.linalg.eigvalsh(GG + PP.conj().T @ np.diag(V0) @ PP)[:3]
            res["sa"][f"s{s}_B{B}"] = dict(disp_ratio_vs_learned=(disp(GG[i], i, k) / disp(H[i], i, k)).tolist(),
                                          bound_rel=(np.abs(e - eb) / np.abs(eb)).tolist(),
                                          tails=[float(np.abs(PP[(i * B + B // 2 + m * B) % N, i])) for m in (1, 2, 3, 4)])
    json.dump(res, open(os.path.join(OUT, "ciclo_p3_dev.json"), "w"), indent=1)
    print(json.dumps(res, indent=1)[:3000])


if __name__ == "__main__":
    main()
