"""surf_f2_kramers.py - POST-HOC (declared B6 in surf_f2_amend.json), questions of the coordinator after phase 2:
(1) full spin polarisation of the exact levels (total and per valley), and how <Sigma_z> = 0.93 lives with v_phi = 0
    (one-tick flux, its Hodge split, the Gordon spin term);
(2) time reversal of the rule: Theta = (Ph C) (i sigma_y on L and R) K, with Theta W2 Theta^-1 = W2^-1 and Theta^2 = -1
    (checked on a small random-potential lattice and on the 48^3 levels);
(3) Ritz from the spin-DOWN seed.
usage: python surf_f2_kramers.py small | levels | ritzdown      -> lab/out/surf_f2_kramers.json (merged)"""
import json
import os
import pickle
import sys
import time
import numpy as np
from surf_core import Walk3D, net_flux, SEED, OUT
from surf_f2_valli import walk, seed, valley_split, sigma_z, CK, N3, A3, X, Y, Z, R, ZA, SIG

JF = os.path.join(OUT, "surf_f2_kramers.json")
ISY = np.array([[0, 1], [-1, 0]], complex)          # i sigma_y


def T(psi):
    """i sigma_y K on the L pair and on the R pair (no L<->R swap)"""
    c = psi.conj()
    return np.concatenate([np.einsum("ab,b...->a...", ISY, c[0:2]), np.einsum("ab,b...->a...", ISY, c[2:4])])


def PhC(psi, theta, V):
    L, Rr = psi[0:2], psi[2:4]
    out = np.concatenate([np.cos(theta) * L - 1j * np.sin(theta) * Rr, np.cos(theta) * Rr - 1j * np.sin(theta) * L])
    return out if V is None else out * np.exp(-1j * V)


def Theta(psi, theta, V):
    return PhC(T(psi), theta, V)


def W2(W, psi):
    W.tick_count = 0
    F1, F2 = [0, 0, 0], [0, 0, 0]
    psi, s = W.tick(psi)
    for (j, A, B) in s:
        F1[j] = F1[j] + net_flux(j, A, B)
    psi1 = psi
    psi, s = W.tick(psi)
    for (j, A, B) in s:
        F2[j] = F2[j] + net_flux(j, A, B)
    W.tick_count = 0
    return psi, psi1, F1, F2


def load():
    return json.load(open(JF)) if os.path.exists(JF) else {}


def save(d):
    json.dump(d, open(JF, "w"), indent=1, default=float)


def small():
    n, th = 4, 0.37
    rng = np.random.default_rng(SEED + 7)
    V = rng.normal(size=(n, n, n))
    W = Walk3D(n, th, V=V)
    err = 0.0
    for _ in range(10):
        x = rng.normal(size=(4, n, n, n)) + 1j * rng.normal(size=(4, n, n, n))
        a, *_ = W2(W, Theta(x, th, V))
        b, *_ = W2(W, a)          # W2 Theta W2 x
        y, *_ = W2(W, x)
        err = max(err, float(np.abs(W2(W, Theta(y, th, V))[0] - Theta(x, th, V)).max()))
    x = rng.normal(size=(4, n, n, n)) + 1j * rng.normal(size=(4, n, n, n))
    sq = Theta(Theta(x, th, V), th, V) / x
    D = 4 * n**3
    M = np.zeros((D, D), complex)
    for k in range(D):
        e = np.zeros(D, complex); e[k] = 1
        M[:, k] = W2(W, e.reshape(4, n, n, n))[0].ravel()
    lam, U = np.linalg.eig(M)
    ph = np.sort(np.angle(lam))
    pair_err = float(np.max(np.abs(ph[0::2] - ph[1::2])))
    gaps = float(np.min(ph[2::2] - ph[1:-1:2]))
    fl = 0.0
    for k in range(D):
        v = U[:, k].reshape(4, n, n, n); v = v / np.linalg.norm(v)
        _, _, F1, F2 = W2(W, v)
        fl = max(fl, max(float(np.abs(F1[j] + F2[j]).max()) for j in range(3)))
    # plain T (no PhC) fails
    y, *_ = W2(W, x)
    plain = float(np.abs(W2(W, T(y))[0] - T(x)).max())
    d = load()
    d["small"] = dict(n=n, theta=theta_s(th), W2_Theta_W2_minus_Theta=err, Theta_sq_real=float(np.mean(sq).real), Theta_sq_imag_max=float(np.abs(sq.imag).max()),
                      plain_T_error=plain, kramers_pair_max_dphase=pair_err, min_gap_between_pairs=gaps,
                      max_net_two_tick_flux_all_eigvecs=fl)
    save(d); print(d["small"])


def theta_s(t):
    return t


def spin_vec(p):
    nrm = float(np.sum(np.abs(p)**2))
    return [sum(float(np.real(np.sum(np.einsum("a...,ab,b...->...", p[q:q + 2].conj(), SIG[k], p[q:q + 2]))))
                for q in (0, 2)) / nrm for k in range(3)]


def vphi_from_flux(F, rho, per_tick=1.0):
    Fc = [0.5 * (F[j] + np.roll(F[j], 1, axis=j)) / per_tick for j in range(3)]
    rc = np.sqrt(X**2 + Y**2)
    m = (R > 1) & (R < 6) & (np.abs(Z) < 0.5 * R)
    J = (-Y * Fc[0] + X * Fc[1]) / rc
    return float(np.sum((J * R / rc)[m]) / np.sum(rho[m]) / ZA)


def hodge(F):
    """F = grad(phi) + S on the periodic link lattice; div on cells: sum_j F_j(x) - F_j(x - e_j)"""
    div = sum(F[j] - np.roll(F[j], 1, axis=j) for j in range(3))
    k = 2 * np.pi * np.fft.fftfreq(N3)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    lap = sum(2 * np.cos(Kc) - 2 for Kc in (KX, KY, KZ)); lap[0, 0, 0] = 1
    phi = np.real(np.fft.ifftn(np.fft.fftn(div) / lap))
    G = [np.roll(phi, -1, axis=j) - phi for j in range(3)]
    return G, [F[j] - G[j] for j in range(3)]


def gordon_spin_vphi(p):
    """continuum Gordon spin term j_s = curl(psibar Sigma psi) / (2m) (c = hbar = m = 1), central differences,
    with beta = L<->R swap in the chiral basis; v_phi from it, same window and normalisation."""
    L, Rr = p[0:2], p[2:4]
    Mv = [np.real(np.einsum("a...,ab,b...->...", L.conj(), SIG[k], Rr) + np.einsum("a...,ab,b...->...", Rr.conj(), SIG[k], L))
          for k in range(3)]
    d = lambda f, ax: (np.roll(f, -1, axis=ax) - np.roll(f, 1, axis=ax)) / (2 * A3)
    jx = d(Mv[2], 1) - d(Mv[1], 2); jy = d(Mv[0], 2) - d(Mv[2], 0)
    rho = np.sum(np.abs(p)**2, 0) / A3**3
    rc = np.sqrt(X**2 + Y**2)
    m = (R > 1) & (R < 6) & (np.abs(Z) < 0.5 * R)
    J = 0.5 * (-Y * jx + X * jy) / rc / A3**3
    return float(np.sum((J * R / rc)[m]) / np.sum(rho[m]) / ZA)


def levels():
    W = walk()
    V = -ZA / R * A3
    d = load(); out = {}
    for nm in ("0.89930", "0.89960"):
        p = np.load(os.path.join(CK, f"f2_level_{nm}.npy")); p = p / np.linalg.norm(p)
        q = Theta(p, A3, V); q = q / np.linalg.norm(q)
        p2, p1, F1, F2 = W2(W, p)
        lam = np.vdot(p, p2)
        q2, *_ = W2(W, q)
        r = dict(omega=float(-np.angle(lam) / (2 * A3)), spin=spin_vec(p),
                 valleys={k: dict(weight=v["weight"], spin=spin_vec(v["comp"])) for k, v in valley_split(p).items()},
                 Theta_psi_resid=float(np.linalg.norm(q2.ravel() - lam * q.ravel())),
                 Theta_psi_overlap_with_psi=abs(complex(np.vdot(p, q))), Theta_psi_spin=spin_vec(q))
        other = "0.89960" if nm == "0.89930" else "0.89930"
        po = np.load(os.path.join(CK, f"f2_level_{other}.npy")); po = po / np.linalg.norm(po)
        r["Theta_psi_overlap_with_other_level"] = abs(complex(np.vdot(po, q)))
        # doublet {p, q}: 2x2 spin matrices, largest polarisation reachable inside the doublet
        B = [p, q]
        Sm = np.array([[[sum(complex(np.sum(np.einsum("a...,ab,b...->...", B[a][s:s + 2].conj(), SIG[k], B[b][s:s + 2])))
                             for s in (0, 2)) for b in range(2)] for a in range(2)] for k in range(3)])
        avec = np.array([np.real(np.trace(Sm[k] @ np.array([[0, 1], [1, 0]]))) / 2 for k in range(3)]), \
            np.array([np.real(np.trace(Sm[k] @ np.array([[0, -1j], [1j, 0]]))) / 2 for k in range(3)]), \
            np.array([np.real(np.trace(Sm[k] @ np.array([[1, 0], [0, -1]]))) / 2 for k in range(3)])
        Amat = np.array(avec)              # 3 (pauli in doublet) x 3 (Sigma_k)
        sv = np.linalg.svd(Amat, compute_uv=False)
        r["doublet_trace_spin"] = [float(np.real(np.trace(Sm[k]))) / 2 for k in range(3)]
        r["doublet_max_polarisation"] = float(sv[0])
        # currents of combinations inside the doublet
        combos = {"psi": p, "Theta_psi": q, "(psi+Theta_psi)/sqrt2": (p + q) / np.sqrt(2),
                  "(psi+i Theta_psi)/sqrt2": (p + 1j * q) / np.sqrt(2)}
        cur = {}
        for cn, c in combos.items():
            c2, c1, G1, G2 = W2(W, c)
            rho = np.sum(np.abs(c)**2, 0)
            cur[cn] = dict(max_net_two_tick=max(float(np.abs(G1[j] + G2[j]).max()) for j in range(3)),
                           max_one_tick=max(float(np.abs(G1[j]).max()) for j in range(3)),
                           vphi_tick1=vphi_from_flux(G1, rho), vphi_tick2=vphi_from_flux(G2, rho),
                           vphi_two_tick=vphi_from_flux([G1[j] + G2[j] for j in range(3)], rho, 2.0))
        r["currents"] = cur
        rho = np.sum(np.abs(p)**2, 0)
        Gg, Ss = hodge(F1)
        r["hodge_tick1"] = dict(vphi_gradient_part=vphi_from_flux(Gg, rho), vphi_solenoidal_part=vphi_from_flux(Ss, rho),
                                norm_gradient=float(sum(np.sum(g**2) for g in Gg)), norm_solenoidal=float(sum(np.sum(s**2) for s in Ss)))
        Gg2, Ss2 = hodge(F2)
        r["hodge_tick2"] = dict(vphi_gradient_part=vphi_from_flux(Gg2, rho), vphi_solenoidal_part=vphi_from_flux(Ss2, rho))
        r["gordon_spin_vphi_full"] = gordon_spin_vphi(p)
        p0 = valley_split(p)["0"]["comp"]; p0 = p0 / np.linalg.norm(p0)
        r["gordon_spin_vphi_valley0"] = gordon_spin_vphi(p0)
        out[nm] = r
        print(nm, json.dumps({k: v for k, v in r.items() if k not in ("valleys",)}, default=float)[:1800], flush=True)
        print("   valleys", {k: (round(v["weight"], 3), [round(s, 3) for s in v["spin"]]) for k, v in r["valleys"].items()}, flush=True)
    sd = seed()
    out["radial_seed_gordon_spin_vphi"] = gordon_spin_vphi(sd)
    d["levels"] = out; save(d)


def ritzdown(budget=540.0):
    """same as surf_f2_valli.ritz, seed = T(spin-up radial 1s) = spin-down 1s; own checkpoint and output names"""
    t_start = time.time()
    f = os.path.join(CK, "f2_ritzdown.pkl")
    Tt = 1200
    cand = np.arange(0.82, 0.94 + 1e-9, 0.005)
    W = walk()
    s0 = T(seed()); s0 = s0 / np.linalg.norm(s0)
    if os.path.exists(f):
        S = pickle.load(open(f, "rb")); W.tick_count = S["t"]
    else:
        S = dict(t=0, psi=s0.copy(), acc=np.zeros((len(cand),) + s0.shape, complex))
    psi, acc = S["psi"], S["acc"]
    while S["t"] < Tt and time.time() - t_start < budget:
        t = S["t"]
        if t % 2 == 0:
            w = np.sin(np.pi * (t + 1.0) / Tt)**2
            ph = w * np.exp(1j * cand * t * A3)
            for k in range(len(cand)):
                acc[k] += ph[k] * psi
        psi, _ = W.tick(psi)
        S["t"] += 1
    S["psi"] = psi
    if S["t"] < Tt:
        pickle.dump(S, open(f, "wb")); print("ritzdown checkpoint", S["t"], flush=True); sys.exit(3)
    B = acc.reshape(len(cand), -1).T
    Q, _ = np.linalg.qr(B)
    UQ = np.zeros_like(Q)
    for k in range(Q.shape[1]):
        UQ[:, k] = W2(walk(), Q[:, k].reshape(4, N3, N3, N3))[0].ravel()
    lam, C = np.linalg.eig(Q.conj().T @ UQ)
    om = -np.angle(lam) / (2 * A3)
    lev = []
    up = {nm: np.load(os.path.join(CK, f"f2_level_{nm}.npy")) for nm in ("0.89930", "0.89960")}
    for i in np.argsort(om):
        if not 0.82 < om[i] < 0.94:
            continue
        v = Q @ C[:, i]; v /= np.linalg.norm(v)
        u = W2(walk(), v.reshape(4, N3, N3, N3))[0].ravel()
        l2 = np.vdot(v, u)
        e = dict(omega=float(om[i]), resid=float(np.linalg.norm(u - l2 * v)), w_seed=float(abs(np.vdot(v, s0.ravel()))**2),
                 spin=spin_vec(v.reshape(4, N3, N3, N3)))
        e["overlap_with_spin_up_levels"] = {nm: abs(complex(np.vdot(x.ravel() / np.linalg.norm(x), v)))**2 for nm, x in up.items()}
        lev.append(e); print(e, flush=True)
    if os.path.exists(f):
        os.remove(f)
    d = load(); d["ritz_spin_down"] = lev; save(d)


if __name__ == "__main__":
    {"small": small, "levels": levels, "ritzdown": ritzdown}[sys.argv[1]]()
