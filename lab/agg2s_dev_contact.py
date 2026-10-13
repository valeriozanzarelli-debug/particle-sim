"""agg2s_dev_contact.py - DEVELOPMENT exploration (declared): how much of the s-state contact (|psi(0)|^2, <p^4>) does a
coarse cell lose, and is the loss UNIVERSAL (the same per unit coarse contact density for 1s, 2s, 3s)? Exact Schroedinger
kinetic k^2/2 + Galerkin field, cells a0/2, a0/4, a0/8. Also the "Born tail" of the nucleus cell: beyond the cube the
coarse rule's own equation gives psi(k) -> 8 pi m Z alpha psi(0)/k^4 (Coulomb input + learned mass only).
Writes lab/out/agg2s_dev_contact.json.  Usage: timeout 590 python lab/agg2s_dev_contact.py DIV N L_a0_tag"""
import os, sys, json, time
import numpy as np
import agg2s_core as c

JS = os.path.join(c.OUT, "agg2s_dev_contact.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}


def cube_I_coef(nsamp=2_000_000, seed=1):
    """I(h) = int_{outside cube [-pi/h, pi/h]^3} d^3k / ((2 pi)^3 k^4) = coef * h ; coef = <max|n_i|>_sphere * 4 pi / (pi (2 pi)^3)"""
    rng = np.random.default_rng(seed)
    n = rng.normal(size=(nsamp, 3))
    n /= np.linalg.norm(n, axis=1)[:, None]
    m = np.abs(n).max(1).mean()
    return float(m * 4 * np.pi / (np.pi * (2 * np.pi) ** 3)), float(m)


def contact(v, N, h_a0):
    """rho_c = |psi_c(0)|^2 (a0^-3), <k^4>_c (a0^-4) of a grid vector (nucleus at index N//2), h in a0"""
    v = v.reshape((N,) * 3)
    nrm = np.sum(v ** 2) * h_a0 ** 3
    rho = v[N // 2, N // 2, N // 2] ** 2 / nrm
    F = np.fft.fftn(v)
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=h_a0)
    k2 = k1[:, None, None] ** 2 + k1[None, :, None] ** 2 + k1[None, None, :] ** 2
    w = np.abs(F) ** 2
    k4 = float(np.sum(w * k2 ** 2) / np.sum(w))
    return float(rho), k4


if __name__ == "__main__":
    div, N = int(sys.argv[1]), int(sys.argv[2])
    states = sys.argv[3].split(",")
    h = c.A0 / div
    coef, mmax = cube_I_coef()
    res["I_coef"] = coef
    f = c.galerkin_field(N, h, cache=False)
    T = c.kin_schroedinger(N, h, (0, 0, 0))
    at = c.TwistAtom(T, f, (0, 0, 0), sym="even")
    Y = None
    vecs = []
    for st in states:
        key = f"d{div}_N{N}_{st}"
        if st == "3s":
            x = f.x
            X, Yg, Z = np.meshgrid(x, x, x, indexing="ij", sparse=True)
            r = np.sqrt(X ** 2 + Yg ** 2 + Z ** 2) / c.A0
            g = np.broadcast_to((1 - 2 * r / 3 + 2 * r ** 2 / 27) * np.exp(-r / 3), (N,) * 3).copy()
            x0 = (g / np.linalg.norm(g)).ravel()
        else:
            x0 = c.guess(at, st)
        v, E, rr, it, ok = c.solve_state(at, x0, Y=None if not vecs else np.stack(vecs, 1), tol_res=1e-9, budget=250)
        vecs.append(v)
        n = int(st[0])
        rho, k4 = contact(v, N, 1.0 / div)
        I = coef / div                                   # a0 units (Z = 1, m = 1)
        rho_born = rho / (1 - 8 * np.pi * I) ** 2
        k4_born = k4 + (8 * np.pi) ** 2 * rho_born * I
        ex_rho = 1 / (np.pi * n ** 3)
        ex_k4 = -3 / n ** 4 + 8 / n ** 3
        out = dict(E_rel=E * c.MC2 / (-c.RY / n ** 2) - 1, residual=rr, rho_c=rho, rho_exact=ex_rho, rho_ratio=rho / ex_rho,
                   rho_born_ratio=rho_born / ex_rho, k4_c=k4, k4_exact=ex_k4, k4_ratio=k4 / ex_k4, k4_born_ratio=k4_born / ex_k4,
                   k4_missing_per_rho=(ex_k4 - k4) / rho, rho_missing_per_rho=(ex_rho - rho) / rho)
        res[key] = out
        json.dump(res, open(JS, "w"), indent=1)
        print(key, " ".join(f"{k}={v:.5g}" for k, v in out.items()), flush=True)
