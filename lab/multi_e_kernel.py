"""multi_e_kernel.py - M0: audit del nucleo dei camminatori e permesso del tempo immaginario (gate G0, esami R1, R2, R3).
Sezioni (ognuna < 10 min, risultati a checkpoint in lab/out/multi_e_kernel.json):
  python lab/multi_e_kernel.py g0      parita' della regola fine (G0)
  python lab/multi_e_kernel.py learn   livello 1 con semi nuovi: catena A (seme 2), B (seme 2), A/2 (semi 1 e 2)
  python lab/multi_e_kernel.py mkin    curvatura -> m_kin(theta), fit in theta^2 (R3)
  python lab/multi_e_kernel.py kbeta   K_beta: peso negativo, somma, D, curtosi, beta_min, scelta del livello (R2)
  python lab/multi_e_kernel.py r1      permesso in 1-D: cammino fine contro livello Pauli 1-D imparato (R1)
"""
import os, sys, time, pickle
from multi_e_core import *                                   # noqa: F401,F403
from multi_e_core import zz, OUT, A0, ALPHA, THETA, Band, chain_levels, load_json, save_json, dirs26
import numpy as np
import scipy.fft as sfft
import agg_zz_core as agg

JS = "multi_e_kernel.json"
res = load_json(JS)
t0 = time.time()
sec = sys.argv[1] if len(sys.argv) > 1 else "g0"
PKL = os.path.join(OUT, "multi_e_levels1.pkl")          # livelli 1 ri-imparati (non committati)


def done(msg):
    print(msg, flush=True)
    save_json(JS, res)


# ------------------------------------------------------------------------------------------------ G0
if sec == "g0":
    th = THETA["A"]
    rng = np.random.default_rng(1)
    u = rng.normal(size=(200, 3)); u /= np.linalg.norm(u, axis=1, keepdims=True)
    k = u * (0.5 * rng.random(200) ** (1 / 3))[:, None]
    Wp, Wm = zz.bloch_batch(k, th), zz.bloch_batch(-k, th)
    ep = np.sort(np.angle(np.linalg.eigvals(Wp)), axis=1)
    em = np.sort(np.angle(np.linalg.eigvals(Wm)), axis=1)
    d_eig = float(np.max(np.abs(ep - em)))
    # parita' come matrice: W(-k) = P W(k) P, P = scambio L <-> R
    P = np.zeros((4, 4)); P[0, 2] = P[1, 3] = P[2, 0] = P[3, 1] = 1
    d_par = float(np.max(np.abs(Wm - P @ Wp @ P)))
    # fase della banda positiva vista da pauli2 dopo 16 tick fini (8 doppi)
    Tp = agg.conv_k(np.linalg.matrix_power(Wp, 8), "dirac4", "pauli2")
    Tm = agg.conv_k(np.linalg.matrix_power(Wm, 8), "dirac4", "pauli2")
    d_band = float(np.max(np.abs(agg.band_phase(Tp, "pauli2") - agg.band_phase(Tm, "pauli2"))))
    # controllo incrociato: Bloch = fine3_tick_fast su onde piane (griglia 16^3, k sulla griglia)
    n = 16
    x = np.arange(n)
    X = np.stack(np.meshgrid(x, x, x, indexing="ij"), 0)
    worst_x = 0.0
    for t in range(6):
        m = rng.integers(-3, 4, 3)
        kk = 2 * np.pi * m / n
        sp = rng.normal(size=4) + 1j * rng.normal(size=4); sp /= np.linalg.norm(sp)
        pw = np.exp(1j * np.tensordot(kk, X, axes=(0, 0)))
        psi = sp[None, :, None, None, None] * pw[None, None]
        out = agg.fine3_tick_fast(agg.fine3_tick_fast(psi, th, 0), th, 1)
        pred = (zz.bloch_batch(kk[None], th)[0] @ sp)[None, :, None, None, None] * pw[None, None]
        worst_x = max(worst_x, float(np.max(np.abs(out - pred))))
    # parte dispari dello stencil imparato (livello 1, catena A), come nel brief: (H(k) - H(-k))/H a kh = 0.01 e 1
    L1 = chain_levels("A")[0]
    bu = Band(L1["K"], L1["offs"], L1["B"], th, sym=False)
    odd = {}
    for kh in (0.01, 0.1, 1.0):
        kap = dirs26() * kh
        gp, gm, g0 = bu.gen(kap), bu.gen(-kap), bu.gen(np.zeros(3))
        odd[str(kh)] = float(np.max(np.abs(gp - gm) / np.abs(0.5 * (gp + gm) - g0)))
    ok = d_eig < 1e-12 and d_band < 1e-12 and worst_x < 1e-12
    res["G0"] = dict(max_eigphase_diff=d_eig, max_parity_matrix_diff=d_par, max_pauli_band_phase_diff=d_band,
                     bloch_vs_tick_planewave=worst_x, learned_stencil_odd_rel=odd, PASS=bool(ok),
                     consequence="parte dispari = artefatto dell'apprendimento: stencil simmetrizzato Re(G_j) (DEBT-8)" if ok
                     else "parte dispari tenuta")
    done(f"G0: eig {d_eig:.1e} parity-matrix {d_par:.1e} band {d_band:.1e} tick-vs-bloch {worst_x:.1e} odd {odd} -> "
         f"{'PASS' if ok else 'FAIL'}")

# ------------------------------------------------------------------------------------------------ learn
if sec == "learn":
    st = pickle.load(open(PKL, "rb")) if os.path.exists(PKL) else {}
    jobs = [("A", 2), ("B", 2), ("A2", 1), ("A2", 2)]
    SUB = {"A": 256, "B": 32, "A2": 512}
    for chain, seed in jobs:
        key = f"{chain}_s{seed}"
        if key in st:
            continue
        if time.time() - t0 > 360:
            print("continua alla prossima chiamata"); break
        TH = THETA[chain]
        kfac = 1 / A0
        rng = np.random.default_rng(20261012 + (0 if seed == 1 else 7777))
        prev = agg.Level("fine3", 3, "dirac4", 4, 1, TH, period=2)
        t1 = time.time()
        lv, tried = agg.auto_level(prev, rng, 20, eps_packet=1e-3, eps_energy=3e-4, T_obs=1000, k_obs=kfac * TH,
                                   bs=(16, 8, 4, 3, 2), ws=(1, 2, 3, 4), bmax_total=SUB[chain], log=print,
                                   ks_rule=lambda B: -kfac * TH * B, learner="auto", packets="phys", energy_mode="state")
        st[key] = dict(K=lv.rule.K, offs=lv.rule.offs, B=lv.B, b=lv.b, w=lv.w, phi0=lv.rule.phi0, check=lv.check,
                       n_tried=len(tried), sec=time.time() - t1)
        pickle.dump(st, open(PKL, "wb"))
        res.setdefault("learned_level1", {})[key] = dict(B=lv.B, b=lv.b, w=lv.w, check=lv.check, n_tried=len(tried),
                                                         sec=time.time() - t1)
        done(f"  {key}: B={lv.B} b={lv.b} w={lv.w} in {time.time() - t1:.0f} s")

# ------------------------------------------------------------------------------------------------ m_kin
def mkin_of(band):
    """curvatura: eps = k^2/(2m) + c4 k^4 su |kappa| in [1e-3, 4e-3]/cella, 26 direzioni"""
    D = dirs26()
    kap = np.linspace(1e-3, 4e-3, 7)
    rows, ys = [], []
    for d in D:
        k_au = (kap[:, None] * d[None, :]) / band.h
        e = band.eps(k_au)
        k2 = np.sum(k_au ** 2, 1)
        rows.append(np.stack([k2, k2 ** 2], 1)); ys.append(e)
    Amat = np.concatenate(rows); y = np.concatenate(ys)
    c = np.linalg.lstsq(Amat, y, rcond=None)[0]
    m = 1 / (2 * c[0])
    # isotropia: rapporto eps/(k^2/2m) a k = 1/a0 lungo asse/faccia/corpo
    iso = {}
    for name, d in (("axis", np.array([1., 0, 0])), ("face", np.array([1., 1, 0]) / 2 ** .5), ("body", np.ones(3) / 3 ** .5)):
        for kk in (1.0, 10.0):
            iso[f"{name}_k{kk:g}"] = float(band.eps(kk * d) / (kk ** 2 / 2))
    return float(m), iso


if sec == "mkin":
    st = pickle.load(open(PKL, "rb"))
    out = {}
    for chain in ("A", "B", "A2"):
        TH = THETA[chain]
        for seed in (1, 2):
            key = f"{chain}_s{seed}"
            if chain in ("A", "B") and seed == 1:
                L = chain_levels(chain)[0]
            else:
                L = st[key]
            for sym in (True, False):
                b = Band(L["K"], L["offs"], L["B"], TH, sym=sym)
                m, iso = mkin_of(b)
                out[f"{key}_{'sym' if sym else 'raw'}"] = dict(theta=TH, m_kin=m, tan_law=float(np.tan(TH) / TH), iso=iso,
                                                            D=1 / (2 * m))
                print(f"{key} {'sym' if sym else 'raw'}: theta {TH:.6f} m_kin {m:.7f} (tan/theta {np.tan(TH) / TH:.7f}) "
                      f"iso k=1 axis {iso['axis_k1']:.6f} k=10 {iso['axis_k10']:.6f} {iso['face_k10']:.6f} {iso['body_k10']:.6f}")
    # fit m(theta) = m0 + a th^2 + b th^4 su tutte le combinazioni di semi
    fits = []
    import itertools
    for sa, sb, sc in itertools.product((1, 2), repeat=3):
        th = np.array([THETA["A2"], THETA["A"], THETA["B"]])
        m = np.array([out[f"A2_s{sa}_sym"]["m_kin"], out[f"A_s{sb}_sym"]["m_kin"], out[f"B_s{sc}_sym"]["m_kin"]])
        V = np.stack([np.ones(3), th ** 2, th ** 4], 1)
        fits.append(np.linalg.solve(V, m))
    fits = np.array(fits)
    m0 = float(fits[:, 0].mean()); sm = float(fits[:, 0].std(ddof=1))
    mA = float(np.mean([out["A_s1_sym"]["m_kin"], out["A_s2_sym"]["m_kin"]]))
    res["mkin"] = out
    res["R3"] = dict(m0=m0, sigma_seeds=sm, fits=fits.tolist(), a=float(fits[:, 1].mean()), b=float(fits[:, 2].mean()),
                     tan_law_a=1 / 3, target="1 +- 5e-5", PASS=bool(abs(m0 - 1) <= 5e-5),
                     m_kin_theta_A_used=float(out["A_s1_sym"]["m_kin"]), m_kin_theta_A_seed_mean=mA,
                     D_theta_A=float(1 / (2 * out["A_s1_sym"]["m_kin"])))
    done(f"R3: m_kin(0) = {m0:.7f} +- {sm:.1e} (a = {fits[:, 1].mean():.4f}, tan law 1/3) -> {'PASS' if abs(m0 - 1) <= 5e-5 else 'FAIL'}")

# ------------------------------------------------------------------------------------------------ K_beta
def stencil_grid(hj, offs, n):
    A = np.zeros((n, n, n))
    idx = tuple(np.asarray(offs)[:, i].astype(int) % n for i in range(3))
    np.add.at(A, idx, hj)
    return np.real(sfft.fftn(A, workers=2).conj())       # H(k) = sum_j H_j e^{i k.j}


if sec == "kbeta":
    betas = [1e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 1e-3, 1e-2, 3e-2]
    levels = chain_levels("A")
    TH = THETA["A"]
    out = res.get("kbeta", {})
    for li, L in enumerate(levels):
        b = Band(L["K"], L["offs"], L["B"], TH, sym=True)
        hj = b.hj_eh()
        offs = np.asarray(L["offs"])
        name = f"B{L['B']}"
        rec = out.get(name, {})
        # stoquasticita': tap fuori diagonale col segno sbagliato (H_j > 0 per j != 0)
        off = np.any(offs != 0, axis=1)
        wrong = hj[off] > 0
        rec["wrong_sign_taps"] = int(np.sum(wrong)); rec["largest_wrong_Eh"] = float(hj[off][wrong].max()) if wrong.any() else 0.0
        rec["diag_Eh"] = float(hj[~off][0]); rec["cell_a0"] = b.h
        for n in (64, 128):
            Hk = stencil_grid(hj, offs, n)
            H0 = float(np.sum(hj))
            for be in betas:
                key = f"n{n}_b{be:g}"
                if key in rec:
                    continue
                Kx = np.real(sfft.ifftn(np.exp(-be * (Hk - H0)), workers=2))
                s = Kx.sum()
                neg = float(-Kx[Kx < 0].sum() / s)
                rec[key] = dict(neg=neg, sumK_minus_1=float(s - 1))
                print(f"{name} n={n} beta={be:g}: neg {neg:.2e} sum-1 {s - 1:.1e}", flush=True)
            del Hk
        # D dalla curvatura (k-space) e dallo spazio reale a beta = 1e-3; curtosi dai cumulanti a beta = 1e-2
        m, _ = mkin_of(b)
        rec["D_curv"] = 1 / (2 * m)
        n = 128
        Hk = stencil_grid(hj, offs, n)
        Kx = np.real(sfft.ifftn(np.exp(-1e-3 * (Hk - np.sum(hj))), workers=2))
        xc = np.fft.fftfreq(n, 1 / n) * b.h
        m2 = np.mean([np.sum(Kx * (xc.reshape([-1 if i == a else 1 for i in range(3)])) ** 2) / Kx.sum() for a in range(3)])
        rec["D_real_beta1e-3"] = float(m2 / (2 * 1e-3))
        rec["sigma_cells_beta1e-3"] = float(np.sqrt(m2) / b.h)
        kk = np.linspace(1e-3, 6e-3, 8) / b.h
        e = b.eps(kk[:, None] * np.array([[1., 0, 0]]))
        c = np.linalg.lstsq(np.stack([kk ** 2, kk ** 4], 1), e, rcond=None)[0]
        rec["c2_c4_axis"] = c.tolist()
        rec["excess_kurtosis_beta1e-2"] = float(-6 * c[1] / (1e-2 * c[0] ** 2))
        out[name] = rec
        res["kbeta"] = out
        save_json(JS, res)
    # beta_min e scelta automatica del livello
    choice = None
    summ = {}
    for name, rec in out.items():
        ok_b = [be for be in betas if be <= 1e-2 and all(rec[f"n{n}_b{be:g}"]["neg"] < 1e-10 for n in (64, 128))]
        # beta_min = il piu' piccolo beta della lista tale che tutti i beta >= lui (fino a 1e-2) passano
        bmin = None
        for be in sorted([x for x in betas if x <= 1e-2]):
            if all(rec[f"n{n}_b{x:g}"]["neg"] < 1e-10 for n in (64, 128) for x in betas if be <= x <= 1e-2):
                bmin = be; break
        maxneg_small = max(rec[f"n{n}_b{be:g}"]["neg"] for n in (64, 128) for be in betas if be <= 0.03)
        sums = max(abs(rec[f"n{n}_b{be:g}"]["sumK_minus_1"]) for n in (64, 128) for be in betas)
        summ[name] = dict(beta_min=bmin, max_neg_beta_le_0p03=maxneg_small, max_abs_sum_minus_1=sums,
                          walker_level=bool(bmin is not None and bmin <= 1e-4))
    walker = [n for n in summ if summ[n]["walker_level"]]
    choice = max(walker, key=lambda n: int(n[1:])) if walker else None
    r1 = out["B16"]
    p_l1 = (summ["B16"]["walker_level"] and summ["B16"]["max_abs_sum_minus_1"] < 1e-9)
    p_D = abs(r1["D_real_beta1e-3"] - r1["D_curv"]) < 2e-5
    p_k = abs(r1["excess_kurtosis_beta1e-2"]) < 0.02
    p_coarse = all(summ[n]["max_neg_beta_le_0p03"] > 1e-3 for n in ("B256", "B512"))
    res["R2"] = dict(summary=summ, chosen_walker_level=choice, D=r1["D_curv"], D_real=r1["D_real_beta1e-3"],
                     kurtosis=r1["excess_kurtosis_beta1e-2"],
                     PASS_level1=bool(p_l1), PASS_D=bool(p_D), PASS_kurtosis=bool(p_k), coarse_levels_non_walker=bool(p_coarse),
                     PASS=bool(p_l1 and p_D and p_k and p_coarse))
    done(f"R2: {summ}\n   scelto {choice}; D curv {r1['D_curv']:.6f} real {r1['D_real_beta1e-3']:.6f}; kurt "
         f"{r1['excess_kurtosis_beta1e-2']:.2e} -> {'PASS' if res['R2']['PASS'] else 'FAIL'}")

# ------------------------------------------------------------------------------------------------ R1
if sec == "r1":
    th, N, V0, w = 0.1, 1024, 0.05, 20.0
    x = np.arange(N) - N // 2
    V = -th * V0 / np.sqrt(1 + (x / w) ** 2)                 # fase per tick = (energia in mc^2) x theta
    # cammino fine 1-D con fase di potenziale a ogni tick: W = e^{-iV} C S  (dt = 1 tick = 1 cella / c)
    Wm = np.zeros((2 * N, 2 * N), complex)
    c, s = np.cos(th), np.sin(th)
    for i in range(N):
        ip, im = (i + 1) % N, (i - 1) % N
        # R in i -> R in i+1 ; L in i -> L in i-1 ; poi giro
        Wm[ip, i] += c; Wm[N + ip, i] += -1j * s           # R(i) va in i+1, contribuisce a R' e L'
        Wm[N + im, N + i] += c; Wm[im, N + i] += -1j * s
    Wm = np.exp(-1j * np.concatenate([V, V]))[:, None] * Wm
    lam, U = np.linalg.eig(Wm)
    Heff = -np.angle(lam)                                    # autofasi in (-pi, pi] / dt
    Uinv = np.linalg.inv(U)
    # stato iniziale: gaussiana a banda positiva (phi = (R+L)/sqrt2, chi piccolo) larga 15 celle
    g = np.exp(-x ** 2 / (2 * 15.0 ** 2))
    psi0 = np.concatenate([g, g]) / np.sqrt(2 * np.sum(g ** 2))
    a0c = Uinv @ psi0
    traj = {}
    for be in (0, 0.5, 1, 2, 5, 10, 20, 50):                 # beta in hbar/mc^2 = beta/theta tick
        wgt = np.exp(-(be / th) * (Heff - Heff.min()))
        a = a0c * wgt
        E = float(np.real(np.sum(np.abs(a) ** 2 * Heff) / np.sum(np.abs(a) ** 2)))   # U quasi unitaria
        traj[str(be)] = E / th                               # energia in mc^2 (m = th per cella, c = 1)
    E_bottom = float(Heff.min() / th)
    # livello Pauli 1-D imparato: banda positiva del cammino fine, b = 16 tick, generatore a stencil w = 4
    b, wst = 16, 4
    rng = np.random.default_rng(3)
    kk = np.concatenate([np.clip(rng.normal(size=3000) * 0.3, -np.pi, np.pi), np.pi * (2 * rng.random(1000) - 1), [0.0]])
    Tb = np.linalg.matrix_power(agg.fine1_bloch(kk / b, th), b)
    lamb, Vb = np.linalg.eig(Tb)
    wphi = agg.phi_weight(Vb)
    o = np.argmax(wphi, axis=-1)
    lp = np.take_along_axis(lamb, o[:, None], 1)[:, 0]
    phi0 = float(np.angle(lp[-1]))
    lg = -np.angle(lp * np.exp(-1j * phi0))
    lg = np.where(lg < -0.5, lg + 2 * np.pi, lg)
    offs = np.arange(-wst, wst + 1)
    sw = np.ones(len(kk)); sw[3000:-1] = 0.03
    Phi = np.cos(np.outer(kk, offs)) * sw[:, None]
    gj = np.linalg.lstsq(Phi, lg * sw, rcond=None)[0]
    # griglia grossa: n = N/b celle grosse, potenziale mediato sulla cella grossa, tau = b tick
    n = N // b
    Vc = V.reshape(n, b).mean(1)
    kc = 2 * np.pi * np.fft.fftfreq(n)
    Gk = np.cos(np.outer(kc, offs)) @ gj
    Gk = Gk - (np.cos(0 * offs) @ gj)
    tau = float(b)
    F = np.fft.ifft(np.eye(n), axis=0) * np.sqrt(n)
    Hc = np.real((F @ np.diag(Gk / tau) @ F.conj().T)) + np.diag(Vc)
    Ee = np.linalg.eigvalsh(Hc)[0]
    # proiezione in tempo immaginario con Strang: e^{-beta V/2} e^{-beta G/tau} e^{-beta V/2}, beta = 0.5 tau
    beta = 0.5 * tau
    f = np.exp(-(np.arange(n) - n // 2) ** 2 / 50.0).astype(complex)
    Eps = []
    for it in range(20000):
        f = np.exp(-beta * Vc / 2) * np.fft.ifft(np.exp(-beta * Gk / tau) * np.fft.fft(np.exp(-beta * Vc / 2) * f))
        nrm = np.linalg.norm(f); f /= nrm
        Eps.append(-np.log(nrm) / beta)
        if it > 50 and abs(Eps[-1] - Eps[-2]) < 1e-15:
            break
    # estrapolazione beta -> 0: ripeto a beta/2, beta/4 e prendo l'energia del generatore sull'autovettore
    Ev = float(np.real(np.vdot(f, Hc @ f)))
    Eproj = {}
    for bb in (beta, beta / 2, beta / 4):
        g2 = f.copy()
        for it in range(40000):
            g2n = np.exp(-bb * Vc / 2) * np.fft.ifft(np.exp(-bb * Gk / tau) * np.fft.fft(np.exp(-bb * Vc / 2) * g2))
            nr = np.linalg.norm(g2n); g2n /= nr
            if np.linalg.norm(g2n - g2) < 1e-14:
                g2 = g2n; break
            g2 = g2n
        Eproj[str(bb)] = float(-np.log(nr) / bb)
    bs_ = np.array([beta, beta / 2, beta / 4]); es = np.array([Eproj[str(x)] for x in bs_])
    cpoly = np.polyfit(bs_, es, 2)
    E0 = float(cpoly[-1])
    relp = abs(E0 - Ee) / abs(Ee)
    p1 = traj["50"] < -1.0
    res["R1"] = dict(fine_energy_vs_beta_mc2=traj, fine_band_bottom_mc2=E_bottom, minus_pi_over_dt_mc2=-np.pi / th,
                     pauli1d_stencil=gj.tolist(), pauli1d_E_eigh=float(Ee), pauli1d_E_proj_by_beta=Eproj,
                     pauli1d_E_proj_beta0=E0, pauli1d_E_mixed_generator=Ev, rel_err=float(relp),
                     PASS_fine_sinks=bool(p1), PASS_pauli=bool(relp < 1e-6), PASS=bool(p1 and relp < 1e-6))
    done(f"R1: fine E(beta) {traj} (fondo {E_bottom:.2f} mc^2); Pauli 1-D proj {E0:.10f} vs eigh {Ee:.10f} rel {relp:.1e} "
         f"-> {'PASS' if res['R1']['PASS'] else 'FAIL'}")
