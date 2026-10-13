"""aggfull_p6.py - AGGLOMERATO_COMPLETO phase 6: the coarse level learned WITH THE FIELD ON, in the language of the
LARGE COMPONENT (the same language as the fine rule and as the zoom chain of nucleo1s), docs/AGGLOMERATO_COMPLETO.md s.14.

Why: NUCLEO_1S phase 8 glued its zoom chain (large component g = (u_L + u_R)/sqrt2 of the zig-zag walk) to the
Taylor-pinned level T (aggfull_relearnT.pkl), which is a free band + the bare potential: no field-motion coupling. The
coupling of phase 2 (aggfull_learn.pkl) was learned in the Foldy-Wouthuysen-aligned basis (unitary descriptor), a
different language from the chain's large component. Here the SAME fine data (first-order response of B fine ticks to a
potential plane wave) are read through the large-component descriptor itself:
    coarse amplitude a(k) = Phi(k) c(k),   Phi = (U_L + U_R)/sqrt2  (U: positive band of the fine walk, FW-aligned)
so the coarse generator is the similarity transform  H_LC = Phi H_FW Phi^-1  and the coupling
    Gamma_LC(k', k) = Phi(k') Gamma_FW(k', k) Phi(k)^-1          (NOT Hermitian: a similarity, same eigenvalues).
Phi(k) is a property of the fine walk (its large-component weight), read off, not put in. Gamma_LC is fitted with a
non-Hermitian sandwich stencil (same candidates and same check budgets as phase 2, aggfull_prereg.json).
Physics read off the stencil (small k, k'):  Gamma = g0 + a k'^2 + b k^2 + c k'.k + i gso sigma.(k' x k)
(Dirac, large component: a = -1/4, b = 0, c = +1/4, gso = +1/4; i.e. V_LC = V + (1/4)(pVp - p^2 V) + spin-orbit).

Commands:  timeout 590 python lab/aggfull_p6.py learn     -> lab/out/aggfull_p6_learn.json / .pkl
Units hbar = c = m_e = 1.
"""
import os, sys, json, pickle, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggfull_core as ac
import agg_zz_core as azc

OUT = ac.OUT
TH = ac.TH_A
DEV = os.environ.get("AGGFULL_DEV") == "1"
PFX = "dev_" if DEV else ""


def herm(M):
    return np.conj(np.swapaxes(M, -1, -2))


def large_component(k):
    """Phi(k) (N,2,2): large component of the fine walk's positive band (FW-aligned basis), physical k"""
    W2 = ac.tick_mats(k * TH, TH)
    W2 = W2[1] @ W2[0]
    U = ac.fw_band(W2)
    return (U[:, 0:2, :] + U[:, 2:4, :]) / np.sqrt(2)


def lc_from_fine(k, kp, B):
    """Gamma_LC(kp, k) per unit potential amplitude (coarse units), from the fine walk run B ticks with the field on"""
    Gk, Gp, phi0, GamFW, D = ac.coarse_from_fine(k, kp, TH, B)
    Pk, Pp = large_component(k), large_component(kp)
    return Pp @ GamFW @ np.linalg.inv(Pk), GamFW, Pk, Pp


def fit_nonherm(kc, kpc, Gam, offs, wts=None, rcond=1e-12):
    """least squares for C (2,2,n,n), no Hermiticity: Gamma(kp,k) ~ sum_jl C_jl e^{-i kp.j} e^{i k.l} (cell units)"""
    n = len(offs)
    w = np.ones(len(kc)) if wts is None else wts
    Phi = np.einsum("nj,nl->njl", np.exp(-1j * kpc @ offs.T), np.exp(1j * kc @ offs.T)).reshape(len(kc), n * n)
    sol = np.linalg.lstsq(Phi * w[:, None], Gam.reshape(len(kc), 4) * w[:, None], rcond=rcond)[0]
    return sol.T.reshape(2, 2, n, n)


def lc_physics(C, offs, h):
    """a (k'^2), b (k^2), c (k'.k), gso, from aggfull_core.physics_of_coupling's basis [1, k^2+k'^2, k.k', k^2-k'^2]"""
    p = ac.physics_of_coupling(C, offs, h)
    a = p["gk"] - p["g_asym"]
    b = p["gk"] + p["g_asym"]
    return dict(g0=p["g0"], a=a, b=b, c=p["gx"], gso=p["gso"], so_shape_resid=p["so_shape_resid"])


def n_of_k(h):
    """the learned large-component amplitude N(k) = Phi(k) (scalar part) at small k: N = 1 + n2 k^2 + n4 k^4"""
    rng = np.random.default_rng(5)
    d = rng.normal(size=(60, 3)); d /= np.linalg.norm(d, axis=1, keepdims=True)
    ks = np.linspace(0.02, 0.6, 30) / h
    rows, vals, spin = [], [], 0.0
    for kk in ks:
        P = large_component(kk * d)
        s = 0.5 * np.real(P[:, 0, 0] + P[:, 1, 1])
        spin = max(spin, float(np.max(np.abs(P - s[:, None, None] * np.eye(2)[None])) / np.max(s)))
        rows.append(kk ** 2); vals.append(s.mean())
    x = np.array(rows); y = np.array(vals)
    c = np.polyfit(x, y, 2)
    return dict(n0=float(c[2]), n2=float(c[1]), n4=float(c[0]), nonscalar_max=spin,
                dirac_rel_dev=float(np.max(np.abs(y / np.sqrt((np.sqrt(1 + x) + 1) / (2 * np.sqrt(1 + x))) - 1))))


def check(C, offs, h, k, kp, Gam, nstate, ku, kpu, eps_abs, eps_rel):
    Gf = ac.eval_coupling(C, offs, k[:nstate] * h, kp[:nstate] * h)
    G00 = ac.eval_coupling(C, offs, np.zeros((1, 3)), np.zeros((1, 3)))[0]
    D00 = Gam[:nstate]
    err_abs = float(np.mean(np.linalg.norm(Gf - D00, axis=(1, 2))) / np.linalg.norm(G00) * np.sqrt(2))
    err_rel = float(np.mean(np.linalg.norm(Gf - D00, axis=(1, 2))) / np.mean(np.linalg.norm(D00 - G00[None], axis=(1, 2))))
    Gu = ac.eval_coupling(C, offs, ku * h, kpu * h)
    tame = float(np.max(np.linalg.norm(Gu, axis=(1, 2)) / np.sqrt(2)))
    return dict(err_abs=err_abs, err_rel=err_rel, tame_max=tame, ok=bool(err_abs <= eps_abs and err_rel <= eps_rel and tame <= 1.5))


def learn():
    pre = json.load(open(os.path.join(OUT, "aggfull_prereg.json")))["engine"]["coupling_check"]
    rng = np.random.default_rng(20261014)
    NP = 600 if DEV else 3000
    out, st = dict(levels=[]), []
    for B in (256, 512):
        t0 = time.time()
        h = B * TH
        k_s = 1 / ac.A0 if B <= 256 else 1 / (2 * ac.A0)

        def pairs(n):
            ng = int(0.75 * n)
            kc = azc.sample_state_k(rng, ng, 3, k_s * h); kpc = azc.sample_state_k(rng, ng, 3, k_s * h)
            ok = np.all(np.abs(kc) < np.pi, 1) & np.all(np.abs(kpc) < np.pi, 1)
            kc, kpc = kc[ok], kpc[ok]
            nu = n - len(kc)
            ku = np.pi * (2 * rng.random((nu, 3)) - 1); kpu = np.pi * (2 * rng.random((nu, 3)) - 1)
            w = np.concatenate([np.ones(len(kc)), 0.03 * np.ones(nu)])
            return np.concatenate([kc, ku]) / h, np.concatenate([kpc, kpu]) / h, w, len(kc)
        k, kp, w, ns = pairs(NP)
        Gam, GamFW, Pk, Pp = lc_from_fine(k, kp, B)
        kh, kph, wh, nsh = pairs(NP // 2)
        Gh = lc_from_fine(kh, kph, B)[0]
        ku = np.pi * (2 * rng.random((2000, 3)) - 1) / h
        kpu = np.pi * (2 * rng.random((2000, 3)) - 1) / h
        G0d = lc_from_fine(np.zeros((1, 3)), np.zeros((1, 3)), B)[0][0]
        dep = np.linalg.norm(Gam - G0d[None], axis=(1, 2))
        wr = w / (dep + 1e-3 * np.median(dep[:ns])) * np.sqrt(np.median(dep[:ns]))
        tried, chosen = [], None
        for name in ("point", "star", "cube"):
            offs = ac.stencil(name)
            C = fit_nonherm(k * h, kp * h, Gam, offs, wts=wr)
            ck = check(C, offs, h, kh, kph, Gh, nsh, ku, kpu, pre["eps_abs"], pre["eps_rel_kdep"])
            tried.append(dict(stencil=name, **ck))
            print(f"B={B} {name:5s} abs {ck['err_abs']:.2e} rel {ck['err_rel']:.2e} tame {ck['tame_max']:.4f} {'OK' if ck['ok'] else ''}", flush=True)
            if ck["ok"]:
                chosen = (name, C, offs); break
        rec = dict(B=B, h=h, candidates=tried)
        if chosen is not None:
            name, C, offs = chosen
            ph = lc_physics(C, offs, h)
            nk = n_of_k(h)
            Gp_ = ac.eval_coupling(C, offs, kh[:nsh] * h, kph[:nsh] * h)
            Gm_ = ac.eval_coupling(C, offs, -kh[:nsh] * h, -kph[:nsh] * h)
            par = float(np.mean(np.linalg.norm(Gp_ - Gm_, axis=(1, 2))) / np.mean(np.linalg.norm(Gp_ + Gm_, axis=(1, 2))))
            # the same fine data read in the FW language, fitted with the phase-2 Hermitian fit (for comparison)
            Cfw = ac.fit_coupling(k * h, kp * h, GamFW, offs, wts=wr)
            pfw = ac.physics_of_coupling(Cfw, offs, h)
            rec.update(chosen=name, parity_odd_rel=par, n_params_real=int(8 * len(offs) ** 2),
                       physics={kk: [complex(v).real, complex(v).imag] if not isinstance(v, float) else v for kk, v in ph.items()},
                       N_of_k=nk, fw_same_data={kk: [complex(v).real, complex(v).imag] if not isinstance(v, float) else v for kk, v in pfw.items()},
                       seconds=time.time() - t0)
            st.append(dict(B=B, h=h, C=C, offs=offs, name=name))
            print(f"  chosen {name}: g0 {ph['g0'].real:.7f} a {ph['a'].real:+.5f} (-0.25) b {ph['b'].real:+.5f} (0) "
                  f"c {ph['c'].real:+.5f} (+0.25) gso {ph['gso'].real:+.5f} (+0.25); imag a,b {ph['a'].imag:+.1e},{ph['b'].imag:+.1e}; "
                  f"N: n2 {nk['n2']:+.5f} (-0.125) nonscalar {nk['nonscalar_max']:.1e} ({time.time() - t0:.0f} s)", flush=True)
        out["levels"].append(rec)
    json.dump(out, open(os.path.join(OUT, PFX + "aggfull_p6_learn.json"), "w"), indent=1, default=float)
    pickle.dump(st, open(os.path.join(OUT, PFX + "aggfull_p6_learn.pkl"), "wb"))


if __name__ == "__main__":
    if sys.argv[1] == "learn":
        learn()
