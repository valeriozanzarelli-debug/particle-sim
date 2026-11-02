"""sc_esporta_stati.py - SCENE_CALCOLATE: export of the helium and H2 states for the program (src/atomo/calcolate.rs).

    python lab/sc_esporta_stati.py          # after lab/sc_elio_base.py and lab/h2v_wf.py (bases in lab/out/, not in git)

Writes assets/scene_calcolate/{elio,h2}/:
    manifest.json       R grid, sizes, the reference numbers of the exams, the vibrational packet
    stati.f64           per R: n_termini x 7 float64 LE (the 6 columns of h2v_core.WF._G, then the weight): the normalised
                        Psi(r1, r2) = sum_k w_k exp(-(G0 r1.r1 + G1 r1.r2 + G2 r2.r2 + G3 z1 + G4 z2 + G5))
    rho1.f32            one-electron density rho1(rho_cyl, z) at the R of parts A and B (integrates to 1), nr x nz
    fotografie.f32      per R: n_foto x 2 x 3 float32 LE, independent samples of |Psi_R|^2 (Metropolis, 200 moves)
    prova_psi.json      200 fixed points and Psi there, from Python (exam S0)
Nothing is computed here that the studies did not compute: same recipe as h2v_core.py.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("OMP_NUM_THREADS", "1")
import h2v_core as C  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DST = os.path.join(REPO, "assets", "scene_calcolate")
N_FOTO = 4000
RMAX, NR, NZ = 6.0, 121, 241


def foto(psi, n, R, rng):
    X = rng.normal(0, 0.8, (n, 2, 3))
    X[:, :, 2] += rng.choice([-R / 2, R / 2], (n, 2))
    p = psi(X) ** 2
    for _ in range(200):
        Y = X + rng.normal(0, 0.35, X.shape)
        q = psi(Y) ** 2
        ok = rng.random(n) * p < q
        X[ok] = Y[ok]
        p[ok] = q[ok]
    return X


def stato(wf):
    """the normalised state: columns of G and weights"""
    p = wf.pairs()
    norm = float(np.sum(p["ww"] * p["S"]))
    wf.psi(np.zeros((1, 2, 3)))
    return wf._G.T, wf.wa / np.sqrt(norm), norm


def tabella(wf):
    W, A, Z = wf._gauss1()
    r = np.linspace(0, RMAX, NR)
    z = np.linspace(-RMAX, RMAX, NZ)
    T = np.zeros((NR, NZ))
    for i in range(0, len(W), 20000):
        w, a, c = W[i:i + 20000], A[i:i + 20000], Z[i:i + 20000]
        T += (w[:, None] * np.exp(-a[:, None] * r[None] ** 2)).T @ np.exp(-a[:, None] * (z[None] - c[:, None]) ** 2)
    return T / 2


def scrivi(nome, wfs, Rs, i_rho, extra):
    d = os.path.join(DST, nome)
    os.makedirs(d, exist_ok=True)
    rng = np.random.default_rng(2026)
    with open(os.path.join(d, "stati.f64"), "wb") as fs, open(os.path.join(d, "fotografie.f32"), "wb") as ff:
        for wf, R in zip(wfs, Rs):
            G, w, norm = stato(wf)
            fs.write(np.column_stack([G, w]).astype("<f8").tobytes())
            ff.write(foto(wf.psi, N_FOTO, R, rng).astype("<f4").tobytes())
            print(nome, "R", R, "termini", len(w), "norma", norm, flush=True)
    T = tabella(wfs[i_rho])
    T.astype("<f4").tofile(os.path.join(d, "rho1.f32"))
    # S0: 200 fixed points
    wf = wfs[i_rho]
    X = np.random.default_rng(7).normal(0, 1.0, (200, 2, 3))
    _, _, norm = stato(wf)
    json.dump(dict(R=Rs[i_rho], X=X.reshape(-1).tolist(), psi=(wf.psi(X) / np.sqrt(norm)).tolist()), open(os.path.join(d, "prova_psi.json"), "w"))
    man = dict(formato="particle-sim/scena-stato-v1", nome=nome, R=list(Rs), i_rho=i_rho, n_termini=[len(w.wa) for w in wfs], n_foto=N_FOTO,
               rho=dict(rmax=RMAX, nr=NR, nz=NZ, massimo=float(T.max())), **extra)
    json.dump(man, open(os.path.join(d, "manifest.json"), "w"), indent=1)


def main():
    if "h2" not in sys.argv:
        he = C.WF(0.0, os.path.join(C.OUT, "sc_he_basis.npz"))
        scrivi("elio", [he], [0.0], 0, dict(Z=[2], rif=dict(r_medio=0.929472, r12_medio=1.42207), energia_Eh=float(he.E_el)))
    if "elio" in sys.argv:
        return
    wfs = [C.WF(R) for R in C.RGRID]
    V = C.vib()
    c0, c1 = V["chi"][:, 0], V["chi"][:, 1]
    x, dx = V["x"], V["dx"]
    mids = np.array([(a + b) / 2 for a, b in zip(C.RGRID[:-1], C.RGRID[1:])])
    idx = np.searchsorted(mids, x)
    # |chi(R, t)|^2 = dens0 + dens1 cos(w t): the weights of the R grid are linear in cos(w t) (tails lumped at the ends, declared)
    d0, d1 = 0.5 * (c0 ** 2 + c1 ** 2), c0 * c1
    w0 = np.bincount(idx, d0 * dx, len(C.RGRID))
    w1 = np.bincount(idx, d1 * dx, len(C.RGRID))
    period_au = 2 * np.pi / (V["E"][1] - V["E"][0])
    wf14 = wfs[C.RGRID.index(1.4)]
    pa = wf14.prob_z(C.images(-np.inf, 0, 0, np.inf)) * 2
    num = wf14.prob_z(C.images(-np.inf, -0.35, 0, np.inf))
    den = wf14.prob_z(C.images(-np.inf, -0.35, -np.inf, np.inf))
    scrivi("h2", wfs, list(C.RGRID), C.RGRID.index(1.4), dict(
        Z=[1, 1], w0=w0.tolist(), w1=w1.tolist(), R0=float(np.sum(x * d0) * dx), R1=float(np.sum(x * d1) * dx),
        periodo_fs=float(period_au * C.AU_FS), abs_z=[w.abs_z() for w in wfs],
        rif=dict(P_opposti=pa, P_altro=num / den), energia_Eh=float(wf14.E_el + 1 / 1.4)))
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
