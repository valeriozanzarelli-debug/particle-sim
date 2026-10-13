"""aggfull_l1.py - PHASE 3 of AGGLOMERATORE COMPLETO: the three problems of the level-1 rule found by the helium study
(docs/MOLTI_ELETTRONI.md on master): (1) systematic mass error, (2) wrong-sign taps (kernel not a probability below
beta_min), (3) the band carries only ~40 % of the k^4 (the '-0.28 mEh' piece).
Learners compared on the SAME samples (agg_zz_core.learn_h data: fine zig-zag walk, b ticks, pauli2 descriptor):
  U  : agg_zz_core.learn_h as it is (relative least squares, uniform share weight 0.03)
  T  : Taylor-pinned: the spin-blind band G(kappa) is fitted with EXACT linear constraints taken from the fine walk itself:
       G(0) = 0, curvature (-2 c2) and quartic (24 c4) along the 3 axes, 6 face and 4 body diagonals
       (c2, c4 from the fine Bloch band at tiny kappa: the teacher, not Dirac); the traceless (spin) part as in learn_h
  M  : T + Markov sign constraint: every off-site tap g_j <= 0 (the hop amplitudes of -G are probabilities: e^{-beta G}
       is then positive at every beta, a 'Metzler' generator)
Usage: python aggfull_l1.py dev  (development numbers)  |  python aggfull_l1.py run (production, after the prereg)"""
import os, sys, json, time, pickle
import numpy as np
import scipy.fft as sfft
from scipy.optimize import lsq_linear
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggfull_core as ac
import agg_zz_core as azc
import zz_core as zz

A0, ALPHA = ac.A0, zz.ALPHA
THETA = {"A": A0 / 16 / 256, "B": A0 / 16 / 32, "A2": A0 / 16 / 512}
OUT = ac.OUT


def dirs13():
    D = [np.eye(3)[i] for i in range(3)]
    D += [np.array(v, float) / np.sqrt(2) for v in ([1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1])]
    D += [np.array(v, float) / np.sqrt(3) for v in ([1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1])]
    return np.array(D)


def fine_band_phase(prev, kap, b):
    """spin-blind positive-band motion phase per coarse step of the level below, at kappa (1/coarse cell), relative to k=0"""
    z = np.zeros((1, 3))
    T = azc.conv_k(prev.mult_k(np.concatenate([kap, z]) / b, b), prev.desc, "pauli2")
    ph = azc.band_phase(T, "pauli2")
    g = -np.angle(np.exp(1j * (ph[:-1] - ph[-1])))
    return np.where(g < -0.5, g + 2 * np.pi, g)


def taylor(prev, b, D):
    """c2(d), c4(d): G(kappa d) = c2 kappa^2 + c4 kappa^4 + c6 kappa^6 + ..., from the level below at tiny kappa"""
    kk = np.array([0.004, 0.006, 0.008, 0.01, 0.013, 0.016, 0.02, 0.025])
    out = []
    for d in D:
        g = fine_band_phase(prev, kk[:, None] * d[None], b)
        V = np.stack([kk ** 2, kk ** 4, kk ** 6, kk ** 8], 1)
        c = np.linalg.lstsq(V, g, rcond=None)[0]
        out.append(c[:2])
    return np.array(out)


def taylor6(prev, b, D, kmax=0.12, npts=24, nfit=7):
    """c2, c4, c6 (phase 4): G(kappa d) = sum_n c_2n kappa^2n fitted with nfit even powers on kappa in [kmax/30, kmax]"""
    kk = np.linspace(kmax / 30, kmax, npts)
    out = []
    for d in D:
        g = fine_band_phase(prev, kk[:, None] * d[None], b)
        V = np.stack([kk ** (2 * n) for n in range(1, nfit + 1)], 1)
        sc = np.abs(V).max(0)
        c = np.linalg.lstsq(V / sc, g, rcond=None)[0] / sc
        out.append(c[:3])
    return np.array(out)


def samples(prev, b, ks, rng, nk, frac_uniform, w_uniform):
    """the learn_h data: state-like + uniform + k = 0; returns kap, Gs (nk, 2, 2), weights"""
    ng = int(nk * (1 - frac_uniform))
    kg = azc.sample_state_k(rng, ng, 3, -ks)
    kg = kg[np.all(np.abs(kg) < np.pi, axis=1)]
    ng = len(kg)
    ku = np.pi * (2 * rng.random((int(nk * frac_uniform), 3)) - 1)
    kap = np.concatenate([kg, ku, np.zeros((1, 3))])
    T = azc.conv_k(prev.mult_k(kap / b, b), prev.desc, "pauli2")
    lam0 = np.linalg.eigvals(T[-1])
    phi0 = float(np.angle(np.mean(lam0)))
    lam, V = np.linalg.eig(T * np.exp(-1j * phi0))
    lg = np.real(1j * np.log(lam / np.abs(lam)))
    lg = np.where(lg < -0.5, lg + 2 * np.pi, lg)
    Gs = V @ (lg[..., :, None] * np.linalg.inv(V))
    sw = np.ones(len(kap)); sw[ng:-1] = w_uniform
    g = np.linalg.norm(Gs, axis=(-2, -1))
    sw = sw / (g + 1e-3 * np.median(g[:ng]))
    return kap, Gs, sw, phi0


def fit(kap, Gs, sw, phi0, w, mode, cons=None, D=None):
    """mode U (= learn_h), T (Taylor-pinned spin-blind part), M (T + off-site taps <= 0)"""
    offs = azc.offsets(w, 3)
    nt = len(offs)
    if mode == "U":
        Phi = np.exp(1j * kap @ offs.T) * sw[:, None]
        rhs = (Gs * sw[:, None, None]).reshape(len(kap), 4)
        sol = np.linalg.lstsq(Phi, rhs, rcond=1e-13)[0]
        G = sol.T.reshape(2, 2, nt)
        jidx = {tuple(o): i for i, o in enumerate(offs)}
        neg = np.array([jidx[tuple(-o)] for o in offs])
        G = 0.5 * (G + np.conj(np.swapaxes(G, 0, 1))[:, :, neg])
        return azc.GenRule(G, offs, 3, 2, phi0)
    # spin-blind part: even real stencil g_j = g_{-j}: unknowns on half the offsets (o > -o lexicographically, plus 0)
    half = [i for i, o in enumerate(offs) if tuple(o) >= tuple(-o)]
    ho = offs[half]
    mult = np.array([1.0 if not np.any(o) else 2.0 for o in ho])
    A = np.cos(kap @ ho.T) * mult[None, :]                          # G_s(kap) = sum_half mult g cos(kap.j)
    y = np.real(0.5 * (Gs[:, 0, 0] + Gs[:, 1, 1]))
    Aw, yw = A * sw[:, None], y * sw
    # exact constraints: G(0) = 0; sum g (d.j)^2 = -2 c2(d); sum g (d.j)^4 = 24 c4(d)
    Crow, crhs = [mult.copy()], [0.0]
    for d, cc in zip(D, cons):
        dj = ho @ d
        Crow.append(mult * dj ** 2); crhs.append(-2 * cc[0])
        Crow.append(mult * dj ** 4); crhs.append(24 * cc[1])
        if len(cc) > 2:                                   # phase 4: k^6 pinned too
            Crow.append(mult * dj ** 6); crhs.append(-720 * cc[2])
    Cm, cv = np.array(Crow), np.array(crhs)
    if mode == "T":
        # null-space solution of min |Aw x - yw| s.t. C x = c
        xp = np.linalg.lstsq(Cm, cv, rcond=None)[0]
        _, s, Vt = np.linalg.svd(Cm)
        r = int(np.sum(s > 1e-10 * s[0]))
        N = Vt[r:].T
        z = np.linalg.lstsq(Aw @ N, yw - Aw @ xp, rcond=1e-13)[0]
        x = xp + N @ z
    else:
        # constraints as heavy rows (relative to the data scale) + bounds g_j <= 0 off site
        scale = np.linalg.norm(Aw) / np.linalg.norm(Cm) * 1e4
        Ab = np.concatenate([Aw, Cm * scale]); yb = np.concatenate([yw, cv * scale])
        lo = np.full(len(ho), -np.inf); hi = np.zeros(len(ho))
        i0 = [i for i, o in enumerate(ho) if not np.any(o)][0]
        hi[i0] = np.inf
        x = lsq_linear(Ab, yb, bounds=(lo, hi), method="bvls", tol=1e-15, max_iter=20000).x
    gs = np.zeros(nt)
    jidx = {tuple(o): i for i, o in enumerate(offs)}
    for i, o in zip(half, ho):
        gs[i] = x[list(half).index(i)]
        gs[jidx[tuple(-o)]] = x[list(half).index(i)]
    # traceless (spin) part: as learn_h on Gs - G_s 1
    Gt = Gs - y[:, None, None] * np.eye(2)[None]
    Phi = np.exp(1j * kap @ offs.T) * sw[:, None]
    sol = np.linalg.lstsq(Phi, (Gt * sw[:, None, None]).reshape(len(kap), 4), rcond=1e-13)[0]
    G = sol.T.reshape(2, 2, nt)
    neg = np.array([jidx[tuple(-o)] for o in offs])
    G = 0.5 * (G + np.conj(np.swapaxes(G, 0, 1))[:, :, neg])
    G = 0.5 * (G + G[:, :, neg])          # parity: G(-k) = G(k) (spin is axial: unchanged), with Hermiticity -> G_j Hermitian
    G = G - 0.5 * (G[0, 0] + G[1, 1])[None, None] * np.eye(2)[:, :, None]
    G = G + gs[None, None] * np.eye(2)[:, :, None]
    return azc.GenRule(G, offs, 3, 2, phi0)


# ------------------------------------------------------------------------------------------- measurements
def band_of(rule):
    g = np.real(0.5 * (rule.K[0, 0] + rule.K[1, 1]))
    offs = np.asarray(rule.offs, float)
    jidx = {tuple(o): i for i, o in enumerate(rule.offs)}
    neg = np.array([jidx[tuple(-o)] for o in rule.offs])
    g = 0.5 * (g + g[neg])                                         # parity-even (symmetrised, as DEBT-8)
    return g, offs


def gband(g, offs, kap):
    return np.cos(kap @ offs.T) @ g - g.sum()


def mass_k4(g, offs, prev, b, h_lc, cfine):
    """m_kin (mean over 13 directions, from c2) relative to the fine band's own, and k^4 fraction along axis/face/body
    (relative to the FINE band's c4; and the fine c4 relative to Dirac -k^4/8)"""
    D = dirs13()
    c2 = np.array([-0.5 * np.sum(g * (offs @ d) ** 2) for d in D])
    c4 = np.array([np.sum(g * (offs @ d) ** 4) / 24 for d in D])
    m_rel = np.mean(cfine[:, 0] / c2)                   # m_kin / m_fine
    m_abs = np.mean(1 / (2 * c2 * h_lc))               # G = tau k^2/2m, tau = h_lc, kap = k h_lc  ->  c2 = 1/(2 m h_lc)
    fr = {nm: float(np.mean(c4[s] / cfine[s, 1])) for nm, s in (("axis", slice(0, 3)), ("face", slice(3, 9)), ("body", slice(9, 13)))}
    fr["c2_aniso"] = float((c2.max() - c2.min()) / c2.mean())          # cubic check of the mass (13 directions)
    return float(m_rel), float(m_abs), fr


def kernel_neg(g, offs, n, betas, tau_au):
    """negative weight of K_beta = IFFT e^{-beta (H(k) - H(0))}, H = G/tau (a.u.), on an n^3 grid"""
    A = np.zeros((n, n, n))
    idx = tuple(np.asarray(offs)[:, i].astype(int) % n for i in range(3))
    np.add.at(A, idx, g)
    Hk = np.real(sfft.fftn(A).conj()) / tau_au
    H0 = g.sum() / tau_au
    out = {}
    for be in betas:
        Kx = np.real(sfft.ifftn(np.exp(-be * (Hk - H0))))
        s = Kx.sum()
        out[f"{be:g}"] = float(-Kx[Kx < 0].sum() / s)
    return out


def measure(rule, prev, b, B, th, cfine, rng, kfac_list=((1 / A0, "Z1"), (10 / A0, "Z10"))):
    h_lc = B * th
    tau_au = h_lc * ALPHA ** 2
    g, offs = band_of(rule)
    m_rel, m_abs, fr = mass_k4(g, offs, prev, b, h_lc, cfine)
    off = np.any(offs != 0, axis=1)
    wrong = int(np.sum(g[off] > 0))
    # parity of the RAW learned rule (before any symmetrisation): ||G(k) - G(-k)|| / ||G(k) - G(0)|| at kh = 0.01, 0.1, 1
    odd = {}
    z = rule.gen_at(np.zeros((1, 3)))[0]
    for kh in (0.01, 0.1, 1.0):
        kv = dirs13() * kh
        gp, gm = rule.gen_at(kv), rule.gen_at(-kv)
        odd[f"{kh:g}"] = float(np.max(np.linalg.norm(gp - gm, axis=(1, 2)) / np.linalg.norm(0.5 * (gp + gm) - z, axis=(1, 2))))
    out = dict(parity_odd_rel=odd, m_over_mfine_minus_1=m_rel - 1, m_kin=m_abs, k4frac=fr, wrong_sign_taps=wrong,
               largest_wrong_over_diag=float(g[off][g[off] > 0].max() / abs(g[~off][0])) if wrong else 0.0)
    out["kernel_neg"] = kernel_neg(g, offs, 64, [1e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 1e-3, 1e-2], tau_au)
    cand = azc.Level("learned", 3, "pauli2", 2, B, th, rule=rule, b=b, w=None)
    for kf, nm in kfac_list:
        ck = azc.kcheck(prev, cand, b, rng, 1000, kf * th, sig=max(1.0, 2.0 / (kf * th * B)), k0_cells=0.5 * kf * th * B,
                        energy_mode="state")
        ok3 = ck["packet"] <= 1e-3 and ck["energy_rel"] <= 3e-5 and ck["rest_rel"] <= 3e-5 and ck["bz_rel90"] <= 0.2 and ck["neg_gen"] <= 0.01
        out[f"check_{nm}"] = {k: v for k, v in ck.items()}
        out[f"check_{nm}"]["ok_budget_3e-5"] = bool(ok3)
    return out


def band_k4_on_state(g, offs, h_lc, Zeff, cfine_mean_c4, rng, n=200000):
    """<G/tau - k^2/2m_kin> on a 1s(Zeff) momentum density (a.u., Eh), and the same for Dirac (-alpha^2 k^4/8)"""
    k = azc.sample_state_k(rng, n, 3, Zeff)            # a.u.
    kap = k * h_lc / A0
    G = gband(g, offs, kap)
    D = dirs13()
    c2 = np.mean([-0.5 * np.sum(g * (offs @ d) ** 2) for d in D])
    e = G / h_lc * A0 ** 2                              # Eh (x 1/alpha^2 handled: mc^2 -> Eh is / alpha^2 = A0^2)
    k2 = np.sum(k ** 2, 1)
    ekin = c2 * (k2 * (h_lc / A0) ** 2) / h_lc * A0 ** 2
    R = float(np.mean(e - ekin))
    Rd = float(np.mean(-ALPHA ** 2 * k2 ** 2 / 8))
    return R, Rd


def prod_one(chain, obs, w, mode, seed, res):
    """one production fit: chain theta, observable Z (1s of size a0/Z), stencil half-width w, learner mode, sampling seed"""
    th = THETA[chain]
    prev = azc.Level("fine3", 3, "dirac4", 4, 1, th, period=2)
    b = B = 16
    h_lc = B * th
    D = dirs13()
    key_f = f"{chain}_fine"
    if key_f not in res:
        cf = taylor(prev, b, D)
        m_fine = float(np.mean(1 / (2 * cf[:, 0] * h_lc)))
        res[key_f] = dict(c2c4=cf.tolist(), m_fine=m_fine, tan_law=float(np.tan(th) / th),
                          c4_over_dirac=float(np.mean(cf[:, 1] / (-1 / (8 * h_lc ** 3 * m_fine ** 3)))))
    cf = np.array(res[key_f]["c2c4"])
    rng = np.random.default_rng(seed)
    kap, Gs, sw, phi0 = samples(prev, b, obs / A0 * th * B, rng, 8000, 0.25, 1e-3)
    r = fit(kap, Gs, sw, phi0, w, mode, cons=cf, D=D)
    out = measure(r, prev, b, B, th, cf, np.random.default_rng(seed + 100),
                  kfac_list=((obs / A0, f"Z{obs:g}"),) + (((10 / A0, "Z10"),) if obs != 10 else ()))
    g, offs = band_of(r)
    off = np.any(offs != 0, axis=1)
    out["wrong_sign_taps_tol"] = int(np.sum(g[off] > 1e-12 * abs(g[~off][0])))
    R, Rd = band_k4_on_state(g, offs, h_lc, 1.6875, None, np.random.default_rng(9))
    out["R_band_He_Zeff1.6875_mEh"], out["R_dirac_He_mEh"] = 1e3 * R, 1e3 * Rd
    # constraint residuals (feasibility of M): mass and quartic of the fitted band vs the fine walk's
    c2 = np.array([-0.5 * np.sum(g * (offs @ d) ** 2) for d in D]); c4 = np.array([np.sum(g * (offs @ d) ** 4) / 24 for d in D])
    out["max_c2_rel_dev"] = float(np.max(np.abs(c2 / cf[:, 0] - 1))); out["max_c4_rel_dev"] = float(np.max(np.abs(c4 / cf[:, 1] - 1)))
    return out, r


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    JS = os.path.join(OUT, "aggfull_l1.json")
    res = json.load(open(JS)) if os.path.exists(JS) else {}
    rules = {}
    # production list (lab/out/aggfull_l1_prereg.json)
    JOBS = [("A", 10, 4, m, s) for s in (1, 2, 3) for m in ("U", "T")]                     # M1, M2, K1 (Z = 10 observable)
    JOBS += [(c, 1, 2, "T", 1) for c in ("A2", "B")] + [("A", 1, 2, "T", 1), ("A", 1, 2, "U", 1)]   # M1 (theta law), level-1 shape
    JOBS += [(c, 1, 2, "U", 1) for c in ("A2", "B")]
    JOBS += [("A", 2, w, "M", s) for w in (2, 3, 4) for s in (1, 2, 3)]                    # P1, P2 (helium scale)
    JOBS += [("A", 2, 4, "T", s) for s in (1, 2, 3)]
    t0 = time.time()
    for chain, obs, w, m, seed in JOBS:
        key = f"{chain}_Z{obs}_w{w}_{m}_s{seed}"
        if key in res:
            continue
        t1 = time.time()
        out, r = prod_one(chain, obs, w, m, seed, res)
        out["sec"] = time.time() - t1
        res[key] = out
        if m == "M" and w == 4 and seed == 1:
            pickle.dump(dict(K=r.K, offs=r.offs, phi0=r.phi0, B=16, b=16, w=4, C=2, desc="pauli2", chain="A", learner="M"),
                        open(os.path.join(OUT, "aggfull_l1_M.pkl"), "wb"))
        if m == "T" and obs == 10 and seed == 1:
            pickle.dump(dict(K=r.K, offs=r.offs, phi0=r.phi0, B=16, b=16, w=4, C=2, desc="pauli2", chain="A", learner="T"),
                        open(os.path.join(OUT, "aggfull_l1_T.pkl"), "wb"))
        json.dump(res, open(JS, "w"), indent=1, default=float)
        ck = out.get("check_Z10", {})
        print(f"{key}: odd(kh=.01) {out['parity_odd_rel']['0.01']:.1e} m/mf-1 {out['m_over_mfine_minus_1']:+.2e} k4 {out['k4frac']['axis']:.4f}/{out['k4frac']['face']:.4f}/{out['k4frac']['body']:.4f} "
              f"wrong {out['wrong_sign_taps_tol']} negK(1e-6) {out['kernel_neg']['1e-06']:.1e} Z10 energy {ck.get('energy_rel', float('nan')):.2e} "
              f"ok3e-5 {ck.get('ok_budget_3e-5')} R_He {out['R_band_He_Zeff1.6875_mEh']:.4f}/{out['R_dirac_He_mEh']:.4f} "
              f"dev c2 {out['max_c2_rel_dev']:.1e} c4 {out['max_c4_rel_dev']:.1e} ({out['sec']:.0f}s)", flush=True)
        if time.time() - t0 > 480:
            print("checkpoint: run again", flush=True)
            sys.exit(0)
    print("L1 ALL DONE", flush=True)
