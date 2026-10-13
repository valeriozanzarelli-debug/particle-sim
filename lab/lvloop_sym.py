"""lvloop_sym.py - which spatial symmetries does the INTERACTING zig-zag rule really have?

Free dispersion of the rule is invariant under all 48 cubic operations (LIMITI_LORENTZ E6). But a loop sees
the full dynamics (eigenvectors, vertices, lattice-scale momenta), not only the free spectrum. Test: put a
generic, asymmetric set of point potentials (the rule's own charge phase e^{-iV dt}, allowed by zz_core) on a
periodic box and compare the spectrum of the period operator for the configuration r_i and for g(r_i), for each
of the 48 signed permutations g. Equal spectra (to 1e-10) <=> g is a symmetry of the dynamics (translations are
exact, so only relative positions matter).

Rules compared (the first is the frozen one; the others are VARIANTS for comparison, not the rule):
  'rule'   : order x,y,z on even ticks, z,y,x on odd ticks (zz_core), period 2 ticks
  'fixed'  : x,y,z every tick, period 1
  'six'    : the 6 permutations in turn xyz,zyx,yzx,xzy,zxy,yxz (the cure "ordine a giro"), period 6
Usage: python lvloop_sym.py  -> lab/out/lvloop_sym.json   (about 1 min, 1 core)
"""
import os, sys, json, time, itertools
import numpy as np
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zz_core as zz

ORDERS = {"rule": [(0, 1, 2), (2, 1, 0)], "fixed": [(0, 1, 2)],
          "six": [(0, 1, 2), (2, 1, 0), (1, 2, 0), (0, 2, 1), (2, 0, 1), (1, 0, 2)]}


def roll_matrix(n, j, s):
    """sparse matrix R with (R f)[i] = f[i - s e_j]  (= np.roll(f, s, axis=j)) on an n^3 periodic box"""
    N = n**3
    idx = np.arange(N).reshape(n, n, n)
    src = np.roll(idx, s, axis=j).ravel()          # new[i] = old[src[i]]
    return sp.csr_matrix((np.ones(N), (np.arange(N), src)), shape=(N, N))


def shift_op(n, j, h):
    lam, P = zz.BASIS_HALF[j]
    N = n**3
    M = sp.csr_matrix((2 * N, 2 * N), dtype=complex)
    for q in range(2):
        proj = np.outer(P[:, q], P[:, q].conj())
        M = M + sp.kron(sp.csr_matrix(proj), roll_matrix(n, j, h * lam[q]), format="csr")
    return M


def tick_op(n, order, theta, V):
    N = n**3
    TL, TR = sp.identity(2 * N, complex, format="csr"), sp.identity(2 * N, complex, format="csr")
    for j in order:
        TL = shift_op(n, j, -1) @ TL
        TR = shift_op(n, j, +1) @ TR
    T = sp.block_diag([TL, TR], format="csr")
    ct, st = np.cos(theta), np.sin(theta)
    I = sp.identity(2 * N, complex, format="csr")
    Cm = sp.bmat([[ct * I, -1j * st * I], [-1j * st * I, ct * I]], format="csr")
    ph = sp.diags(np.tile(np.exp(-1j * V.ravel()), 4))      # dt = a = 1 here: phase e^{-iV}
    return ph @ Cm @ T


def period_op(n, rule, theta, V):
    W = sp.identity(4 * n**3, complex, format="csr")
    for order in ORDERS[rule]:
        W = tick_op(n, order, theta, V) @ W
    return W.toarray()


def check_against_zz(n=6, theta=0.37, seed=1):
    rng = np.random.default_rng(seed)
    V = rng.random((n, n, n))
    wk = zz.Walk("dirac", n, 1.0, theta=theta, V=V / 1.0)    # zz phase = exp(-i V a), a = 1
    psi = rng.normal(size=(4, n, n, n)) + 1j * rng.normal(size=(4, n, n, n))
    out = wk.tick(wk.tick(psi.copy()))
    W = period_op(n, "rule", theta, V)
    return float(np.max(np.abs(W @ psi.ravel() - out.ravel())))


def spectrum(W):
    ph = np.sort(np.angle(np.linalg.eigvals(W)))
    return ph


def spec_dist(a, b):
    # both sorted angles in (-pi, pi]; compare as multisets on the circle via sorted lists (generic, no wrap issue
    # when no eigenvalue sits at +-pi within 1e-8; checked)
    return float(np.max(np.abs(a - b)))


def signed_perms():
    out = []
    for p in itertools.permutations(range(3)):
        for s in itertools.product((1, -1), repeat=3):
            G = np.zeros((3, 3), int)
            for i in range(3):
                G[i, p[i]] = s[i]
            out.append(G)
    return out


def config_potential(n, pts, vals):
    V = np.zeros((n, n, n))
    for r, v in zip(pts, vals):
        V[tuple(np.mod(r, n))] += v
    return V


def main():
    t0 = time.time()
    res = {"check_vs_zz_core": check_against_zz()}
    theta = 0.37
    pts = np.array([[0, 0, 0], [1, 2, 0], [0, 1, 3], [2, 0, 1]])
    vals = [0.9, 0.55, 0.31, 0.73]
    G48 = signed_perms()
    for n in (6, 8):
        for rule in ("rule", "fixed", "six"):
            if n == 8 and rule == "six":
                pass
            W0 = period_op(n, rule, theta, config_potential(n, pts, vals))
            s0 = spectrum(W0)
            rows = []
            for G in G48:
                Vg = config_potential(n, pts @ G.T, vals)
                d = spec_dist(s0, spectrum(period_op(n, rule, theta, Vg)))
                rows.append({"G": G.tolist(), "dist": d})
            sym = [r for r in rows if r["dist"] < 1e-9]
            # describe the symmetric subgroup: which axis permutations appear
            perms = sorted({tuple(int(np.argmax(np.abs(np.array(r["G"])[i]))) for i in range(3)) for r in sym})
            res[f"{rule}_n{n}"] = {"n_sym": len(sym), "axis_perms_in_group": perms,
                                   "min_dist_nonsym": float(min([r["dist"] for r in rows if r["dist"] >= 1e-9], default=0)),
                                   "max_dist": float(max(r["dist"] for r in rows))}
            print(rule, n, len(sym), perms, "%.1fs" % (time.time() - t0), flush=True)
            if n == 8 and rule == "rule":
                break
        if n == 8:
            break
    # free check: with V = 0 everywhere, the spectrum is trivially invariant; with a single impurity too (translations)
    # pair-separation test, rule, n=8: separations along the three axes
    n = 8
    sep = {}
    for lab, d in (("x2", (2, 0, 0)), ("y2", (0, 2, 0)), ("z2", (0, 0, 2)),
                   ("xy", (1, 1, 0)), ("yz", (0, 1, 1)), ("xz", (1, 0, 1)),
                   ("x1", (1, 0, 0)), ("y1", (0, 1, 0)), ("z1", (0, 0, 1))):
        V = config_potential(n, [np.zeros(3, int), np.array(d)], [0.8, 0.8])
        sep[lab] = spectrum(period_op(n, "rule", theta, V))
    res["pair_rule_n8"] = {f"{a}-{b}": spec_dist(sep[a], sep[b]) for a, b in
                           (("x2", "z2"), ("x2", "y2"), ("z2", "y2"), ("xy", "yz"), ("xy", "xz"),
                            ("x1", "z1"), ("x1", "y1"))}
    print(res["pair_rule_n8"])
    res["time_s"] = time.time() - t0
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "lvloop_sym.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
