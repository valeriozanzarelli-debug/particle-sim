"""pre_p4.py - P4: the "two sheets" (space doubling) of the zig-zag rule. Is it the lattice fermion doubling?
 a) Bloch symmetries of the one-tick U1(k) (fixed order) and of the two-tick alternating walk.
 b) Dirac points and the winding (2-D chirality) of each cone, R and L blocks, theta = 0.
 c) Wilson-like fix: coin angle theta(k) = theta + r (2 - cos kx - cos ky) (k-space; quasi-local). Vacuum MI on
    a 16x16 torus -> geo_core blind reconstruction: how many pieces?
 d) Sheet labels: vacuum MI of sheet B shifted by one site vs sheet A; Chern number of the filled bands.
Pre-registered in lab/out/pre_prereg.json. Writes lab/out/pre_p4.json (~1 min).
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import geo_core as G
from zz_core import T_axis, bloch_batch, PAULI

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
PI = np.pi
X4 = np.block([[np.zeros((2, 2)), np.eye(2)], [np.eye(2), np.zeros((2, 2))]]).astype(complex)


def U1(kx, ky, theta=0.5, r=0.0):
    TL = T_axis(ky, 1, -1, 0.5) @ T_axis(kx, 0, -1, 0.5)
    TR = T_axis(ky, 1, +1, 0.5) @ T_axis(kx, 0, +1, 0.5)
    T = np.zeros((4, 4), complex); T[:2, :2] = TL; T[2:, 2:] = TR
    th = theta + r * (2 - np.cos(kx) - np.cos(ky))
    Cm = np.cos(th) * np.eye(4) - 1j * np.sin(th) * X4
    return Cm @ T


def logm_u(Ub):
    w, V = np.linalg.eig(Ub)
    return V @ np.diag(1j * np.log(w)) @ np.linalg.inv(V)      # H with U = exp(-i H)


def winding(block, k0, q=1e-4):
    """2-D chirality of the cone at k0 for a 2x2 block U_b(k) = phase * exp(-i d(q).sigma)"""
    # bug fixed after the first run: sqrt(det) has a branch ambiguity at omega = pi; use the fixed phase of U(k0)
    ph = np.trace(block(*k0)) / 2
    ph = ph / abs(ph)

    def dvec(qx, qy):
        Ub = block(k0[0] + qx, k0[1] + qy)
        Hm = logm_u(Ub / ph)
        return np.array([np.real(np.trace(Hm @ s)) / 2 for s in PAULI])
    dx = (dvec(q, 0) - dvec(-q, 0)) / (2 * q)
    dy = (dvec(0, q) - dvec(0, -q)) / (2 * q)
    # orientation of the map q -> d (projected on the plane spanned by dx, dy)
    nrm = np.cross(dx, dy)
    return dx.tolist(), dy.tolist(), float(np.sign(nrm[np.argmax(np.abs(nrm))])) if np.linalg.norm(nrm) > 1e-6 else 0.0, nrm.tolist()


def vacuum_bloch(n, theta, r):
    """real-space vacuum correlation matrix on an n x n torus from the Bloch projectors"""
    ks = 2 * PI * np.arange(n) / n
    P = np.zeros((n, n, 4, 4), complex)
    for i, kx in enumerate(ks):
        for j, ky in enumerate(ks):
            U = U1(kx, ky, theta, r)
            Hs = (U - U.conj().T) / 2j
            w, V = np.linalg.eigh(Hs)
            Vo = V[:, w > 0]
            P[i, j] = Vo @ Vo.conj().T
    # C[a,(x,y); b,(x',y')] = 1/N sum_k e^{i k.(r - r')} P_ab(k)
    c = np.fft.ifft2(P, axes=(0, 1))            # c[dx,dy] = 1/N sum_k e^{+i k d} P(k)
    Ns = n * n
    xs = np.arange(n)
    DX = (xs[:, None, None, None] - xs[None, None, :, None]) % n
    DY = (xs[None, :, None, None] - xs[None, None, None, :]) % n
    blk = c[DX, DY]                              # (x, y, x', y', 4, 4)
    C = blk.reshape(Ns, Ns, 4, 4).transpose(2, 0, 3, 1).reshape(4 * Ns, 4 * Ns)
    return C


def chern(theta, r, nk=48):
    ks = 2 * PI * np.arange(nk) / nk
    V = np.empty((nk, nk, 4, 2), complex)
    for i, kx in enumerate(ks):
        for j, ky in enumerate(ks):
            U = U1(kx, ky, theta, r)
            w, vv = np.linalg.eigh((U - U.conj().T) / 2j)
            V[i, j] = vv[:, w > 0]
    tot = 0.0
    for i in range(nk):
        for j in range(nk):
            a, b, c, d = V[i, j], V[(i + 1) % nk, j], V[(i + 1) % nk, (j + 1) % nk], V[i, (j + 1) % nk]
            l = (np.linalg.det(a.conj().T @ b) * np.linalg.det(b.conj().T @ c) *
                 np.linalg.det(c.conj().T @ d) * np.linalg.det(d.conj().T @ a))
            tot += np.angle(l)
    return tot / (2 * PI)


if __name__ == "__main__":
    t0 = time.time()
    rng = np.random.default_rng(1)
    res = {}
    K = rng.uniform(-PI, PI, (200, 2))
    d_pp = max(np.abs(U1(kx + PI, ky + PI) - U1(kx, ky)).max() for kx, ky in K)
    d_p0 = max(np.abs(U1(kx + PI, ky) + U1(kx, ky)).max() for kx, ky in K)
    K3 = np.c_[K, np.zeros(len(K))]
    W = bloch_batch(K3, 0.5)
    d_alt_pp = float(np.abs(bloch_batch(K3 + [PI, PI, 0], 0.5) - W).max())
    d_alt_p0 = float(np.abs(bloch_batch(K3 + [PI, 0, 0], 0.5) - W).max())
    res["P4a"] = dict(U1_pipi_minus_U1=float(d_pp), U1_pi0_plus_U1=float(d_p0))
    res["P4b"] = dict(alt2_pipi_minus_W=d_alt_pp, alt2_pi0_minus_W=d_alt_p0)
    print(res["P4a"], res["P4b"], flush=True)
    # cones, theta = 0
    cones = {}
    for nm, k0 in (("(0,0)", (0, 0)), ("(pi,pi)", (PI, PI)), ("(pi,0)", (PI, 0)), ("(0,pi)", (0, PI))):
        for blkname, sl in (("R", slice(2, 4)), ("L", slice(0, 2))):
            f = lambda kx, ky, sl=sl: U1(kx, ky, 0.0)[sl, sl]
            ev = np.linalg.eigvals(f(*k0))
            dx, dy, wsign, nrm = winding(f, k0)
            cones[f"{blkname}{nm}"] = dict(omega=[float(-np.angle(e)) for e in ev], winding=wsign, normal=nrm)
    res["P4d_cones"] = cones
    print(json.dumps(cones), flush=True)
    # gaps at the special points with Wilson
    gaps = {}
    for r in (0.0, 0.5):
        for nm, k0 in (("(0,0)", (0, 0)), ("(pi,pi)", (PI, PI)), ("(pi,0)", (PI, 0))):
            om = np.sort(-np.angle(np.linalg.eigvals(U1(*k0, 0.5, r))))
            gaps[f"r{r}{nm}"] = om.tolist()
    res["gaps"] = gaps
    # vacuum MI: check Bloch construction against geo_core (r = 0), then Wilson
    n = 16
    shp = G.torus(n, n)
    Cg, gap, nocc = G.vacuum(G.rule_unitary(shp, 0.5))
    Cb = vacuum_bloch(n, 0.5, 0.0)
    res["bloch_vs_geo_core_C"] = float(np.abs(Cb - Cg).max())
    print("Bloch vs geo_core vacuum:", res["bloch_vs_geo_core_C"], flush=True)
    for r in (0.0, 0.5):
        C = vacuum_bloch(n, 0.5, r)
        I = G.mutual_information(C, n * n)
        nc, rs = G.reconstruct_from_I(I, maxdim=2)
        res[f"P4c_wilson_r{r}"] = dict(pieces=int(nc), betti=[[int(b) for b in q["betti"]] for q in rs],
                                       D=[float(q["dim"]) for q in rs],
                                       I_between_sheets_max=float(np.max(I[np.ix_(
                                           np.where((shp.xy.sum(1) % 2) == 0)[0], np.where((shp.xy.sum(1) % 2) == 1)[0])])))
        print(f"Wilson r={r}", res[f"P4c_wilson_r{r}"], f"{time.time()-t0:.0f}s", flush=True)
        if r == 0.0:
            # sheet B shifted by (1,0) vs sheet A
            idx = shp.index
            A = np.where((shp.xy.sum(1) % 2) == 0)[0]
            Bm = np.array([idx[(shp.xy[i, 0] + 1) % n, shp.xy[i, 1]] for i in A])
            res["P4d_sheet_shift_maxdiff"] = float(np.abs(I[np.ix_(A, A)] - I[np.ix_(Bm, Bm)]).max())
            print("sheet shift diff", res["P4d_sheet_shift_maxdiff"], flush=True)
    # POST-HOC: Wilson strength scan (the pre-registered r = 0.5 left two geometric clusters)
    for r in (1.0, 2.0):
        C = vacuum_bloch(n, 0.5, r)
        I = G.mutual_information(C, n * n)
        nc, rs = G.reconstruct_from_I(I, maxdim=2)
        res[f"posthoc_wilson_r{r}"] = dict(pieces=int(nc), betti=[[int(b) for b in q["betti"]] for q in rs],
                                           D=[float(q["dim"]) for q in rs])
        print(f"post-hoc Wilson r={r}", res[f"posthoc_wilson_r{r}"], flush=True)
    res["chern"] = {f"r{r}": float(chern(0.5, r)) for r in (0.0, 0.5)}
    print("chern", res["chern"])
    wR = [cones[f"R{p}"]["winding"] for p in ("(0,0)", "(pi,pi)", "(pi,0)", "(0,pi)")]
    wL = [cones[f"L{p}"]["winding"] for p in ("(0,0)", "(pi,pi)", "(pi,0)", "(0,pi)")]
    w0 = res["P4c_wilson_r0.5"]
    res["verdict"] = {
        "P4a_sheet_copy_and_pi_doubler": bool(d_pp < 1e-12 and d_p0 < 1e-12),
        "P4b_alternation_keeps_sheets": bool(d_alt_pp < 1e-12),
        "P4c_wilson_one_piece_b01=1,2": bool(w0["pieces"] == 1 and w0["betti"][0][:2] == [1, 2]),
        "P4d_same_winding_R_all_L_opposite": bool(len(set(wR)) == 1 and len(set(wL)) == 1 and wR[0] == -wL[0]),
        "P4d_sheets_identical_by_translation": bool(res["P4d_sheet_shift_maxdiff"] < 1e-10),
    }
    res["time_s"] = time.time() - t0
    print(json.dumps(res["verdict"], indent=1), res["time_s"])
    with open(os.path.join(OUT, "pre_p4.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
