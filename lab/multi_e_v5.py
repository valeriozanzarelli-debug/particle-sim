"""multi_e_v5.py - corse della pre-registrazione v5 (lab/out/multi_e_prereg_v5.json): livello 1 T6 con la scelta di fase 5 / M6 (AGGLOMERATO_COMPLETO par. 12-13).
  learn   impara T6 (A, Z = 10, semi 1, 3, 4; w = il piu' piccolo fra 4..7 che passa TUTTO il controllo, esame su un campione
          indipendente, come aggfull_p5.py) / M6 (A, Z = 2, w = 8, come aggfull_p4.py) / T6 (A2, A, B, Z = 1, w = 5)
          (funzioni di aggfull_l1 e agg_zz_core importate, non modificate) -> lab/out/multi_e_v5_<key>.pkl, multi_e_v5.json
  checks  R3v5, L1c_v5, G0v5, R2v5
  tables  tabelle F(s) delle bande T6 e M6
  he      M1v5, E1v5, Q1v5
Ogni sezione riprende da sola (checkpoint). Run: python lab/multi_e_v5.py learn | checks | tables | he"""
import os, sys, time, pickle, itertools
import numpy as np
import scipy.fft as sfft
from multi_e_core import OUT, ALPHA, A0, THETA, EH_EV, NUC_MASS, load_json, save_json, Band, fib_dirs, dirs26

JS = "multi_e_v5.json"
res = load_json(JS)
what = sys.argv[1] if len(sys.argv) > 1 else "learn"
TH = {"A": THETA["A"], "B": THETA["B"], "A2": THETA["A2"]}
KT, KM = "A_Z10_T6_s1", "A_Z2_w8_M6_s1"
JOBS = [("A", 10, None, "T", s) for s in (1, 3, 4)] + [("A", 2, 8, "M", 1)] + [(c, 1, 5, "T", 1) for c in ("A2", "A", "B")]


def pkl(key):
    return os.path.join(OUT, f"multi_e_v5_{key}.pkl")


def level(key):
    return pickle.load(open(pkl(key), "rb"))


def check(prev, th, kf, r, w, seed, nk):
    """controllo di aggfull_p5.check (stesse soglie)"""
    import agg_zz_core as azc
    cand = azc.Level("learned", 3, "pauli2", 2, 16, th, rule=r, b=16, w=w)
    ck = azc.kcheck(prev, cand, 16, np.random.default_rng(seed), 1000, 10 / A0 * th, sig=max(1.0, 2.0 / kf), k0_cells=0.5 * kf,
                    nk=nk, energy_mode="state")
    ok = ck["packet"] <= 1e-3 and ck["bz_rel90"] <= 0.2 and ck["neg_gen"] <= 0.01 and ck["energy_rel"] <= 3e-5 and ck["rest_rel"] <= 3e-5
    return ck, bool(ok)


if what == "learn":
    import aggfull_l1 as L
    import agg_zz_core as azc
    lr = res.get("learn", {}); c246s = res.get("c246", {})
    t0 = time.time()
    B = b = 16
    D = L.dirs13()
    for chain, obs, w0, m, seed in JOBS:
        key = f"{chain}_Z{obs}_T6_s{seed}" if w0 is None else f"{chain}_Z{obs}_w{w0}_{m}6_s{seed}"
        if key in lr:
            continue
        t1 = time.time()
        th = TH[chain]
        prev = azc.Level("fine3", 3, "dirac4", 4, 1, th, period=2)
        if chain not in c246s:
            c246s[chain] = L.taylor6(prev, b, D).tolist()
            res = load_json(JS); res["c246"] = c246s; save_json(JS, res)
        c246 = np.array(c246s[chain])
        kf = obs / A0 * th * B
        kap, Gs, sw, phi0 = L.samples(prev, b, kf, np.random.default_rng(seed), 8000, 0.25, 1e-3)
        out = dict(tried=[])
        chosen = None
        for w in ((4, 5, 6, 7) if w0 is None else (w0,)):
            r = L.fit(kap, Gs, sw, phi0, w, m, cons=c246, D=D)
            if obs != 10:
                chosen = (w, r)
                break
            ck, ok = check(prev, th, kf, r, w, seed + 100, 3000)
            out["tried"].append(dict(w=w, energy=ck["energy_rel"], rest=ck["rest_rel"], bz=ck["bz_rel90"], neg=ck["neg_gen"], ok=ok))
            print(key, f"w={w}: scelta energy {ck['energy_rel']:.3e} ok {ok} ({time.time() - t1:.0f}s)", flush=True)
            if ok:
                chosen = (w, r)
                break
        if chosen is None:
            out["w"] = None; lr[key] = out
            res = load_json(JS); res["learn"] = lr; save_json(JS, res)
            print(key, "nessun w passa: FAIL dichiarato", flush=True)
            continue
        w, r = chosen
        if obs == 10:
            ck, ok = check(prev, th, kf, r, w, seed + 200, 10000)
            out["exam"] = dict(check=ck, ok=ok)
            print(key, f"w={w} ESAME indipendente energy {ck['energy_rel']:.3e} rest {ck['rest_rel']:.2e} ok {ok}", flush=True)
        g, o = L.band_of(r)
        m_rel, m_abs, fr = L.mass_k4(g, o, prev, b, B * th, c246[:, :2])
        off = np.any(o != 0, axis=1)
        dev = [float(np.max(np.abs(np.array([s_ * np.sum(g * (o @ d) ** p) / f for d in D]) / c246[:, i] - 1)))
               for i, (p, f, s_) in enumerate(((2, 2, -1), (4, 24, 1), (6, 720, -1)))]
        out.update(w=w, m_over_mfine_minus_1=float(m_rel - 1), k4frac=fr, max_dev_c2_c4_c6=dev,
                   wrong_sign_taps=int(np.sum(g[off] > 1e-12 * abs(g[~off][0]))), sec=time.time() - t1)
        pickle.dump(dict(K=r.K, offs=r.offs, phi0=r.phi0, B=16, b=16, w=w, C=2, desc="pauli2", chain=chain, learner=m + "6"),
                    open(pkl(key), "wb"))
        lr[key] = out
        res = load_json(JS); res["learn"] = lr; res["c246"] = c246s; save_json(JS, res)
        print(key, f"{out['sec']:.0f}s w={w} m-1={m_rel - 1:.1e} dev={dev}", flush=True)
        if time.time() - t0 > 400:
            print("checkpoint: run again", flush=True); sys.exit(0)
    print("learn done", flush=True)

if what == "checks":
    from multi_e_kernel import mkin_of, stencil_grid
    out = {}
    # R3v5
    mk = {}
    for c in ("A2", "A", "B"):
        for s in (1,):
            L = level(f"{c}_Z1_w5_T6_s1")
            mk[f"{c}_s{s}"] = mkin_of(Band(L["K"], L["offs"], 16, TH[c], sym=True))[0]
    fits = []
    for sa, sb, sc in ((1, 1, 1),):
        th = np.array([TH["A2"], TH["A"], TH["B"]])
        mm = np.array([mk[f"A2_s{sa}"], mk[f"A_s{sb}"], mk[f"B_s{sc}"]])
        fits.append(np.linalg.solve(np.stack([np.ones(3), th ** 2, th ** 4], 1), mm))
    fits = np.array(fits); m0 = float(fits[:, 0].mean())
    # un solo seme (massa fissata): la barra e' la deriva di m dalla regola fine, |m/m_fine - 1| (dichiarato in v5)
    sm = float(max(abs(v["m_over_mfine_minus_1"]) for v in res["learn"].values() if v.get("w")) + 1e-12)
    LT = level(KT); LM = level(KM)
    bT = Band(LT["K"], LT["offs"], 16, TH["A"]); bM = Band(LM["K"], LM["offs"], 16, TH["A"])
    mT, isoT = mkin_of(bT); mM, isoM = mkin_of(bM)
    out["R3v5"] = dict(m_kin=mk, m0=m0, sigma=sm, a=float(fits[:, 1].mean()), PASS=bool(abs(m0 - 1) <= 5e-5),
                       within_1e9=bool(abs(m0 - 1) <= 1e-9), m_T_thetaA=mT, m_M_thetaA=mM, tan_law_A=float(np.tan(TH["A"]) / TH["A"]))
    print("R3v5", out["R3v5"], flush=True)
    # L1c_v5
    lr = res["learn"]
    P5 = load_json("aggfull_p5.json")
    theirs = {s: P5[f"s{s}"] for s in range(1, 7)}
    ours = {s: lr[f"A_Z10_T6_s{s}"] for s in (1, 3, 4)}
    out["L1c_v5"] = dict(theirs=dict(w={s: v["w"] for s, v in theirs.items()},
                                     energy_rel={s: v["exam_check"]["energy_rel"] for s, v in theirs.items()},
                                     n_pass=int(sum(v["exam_ok"] for v in theirs.values()))),
                         ours=dict(w={s: v["w"] for s, v in ours.items()},
                                   energy_rel={s: v["exam"]["check"]["energy_rel"] for s, v in ours.items()},
                                   rest_rel={s: v["exam"]["check"]["rest_rel"] for s, v in ours.items()},
                                   ok={s: v["exam"]["ok"] for s, v in ours.items()}),
                         reproduces=bool(all(ours[s]["w"] == theirs[s]["w"] and abs(ours[s]["exam"]["check"]["energy_rel"] -
                                             theirs[s]["exam_check"]["energy_rel"]) <= 1e-7 for s in ours)))
    out["L1c_v5"]["PASS"] = bool(all(ours[s]["exam"]["ok"] for s in ours) and all(v["exam_ok"] for v in theirs.values()))
    print("L1c_v5", out["L1c_v5"], flush=True)
    # G0v5: parita'
    g0 = {}
    rng = np.random.default_rng(3)
    D = fib_dirs(200)
    for nm, L in (("T6", LT), ("M6", LM)):
        raw = Band(L["K"], L["offs"], 16, TH["A"], sym=False)
        for kh in (0.01, 1.0):
            kap = kh * D
            gp = raw.gen(kap) - raw.gen(np.zeros(3)); gm = raw.gen(-kap) - raw.gen(np.zeros(3))
            g0[f"{nm}_kh{kh:g}"] = float(np.max(np.abs(gp - gm)) / np.max(np.abs(gp + gm)))
        # parte immaginaria dello stencil (parte dispari in k)
        g0[f"{nm}_imag_over_real"] = float(np.max(np.abs(np.imag(0.5 * (L["K"][0, 0] + L["K"][1, 1])))) /
                                           np.max(np.abs(np.real(0.5 * (L["K"][0, 0] + L["K"][1, 1])))))
    out["G0v5"] = dict(odd_rel=g0, PASS=bool(all(v < 1e-10 for k, v in g0.items() if "_kh" in k)))
    print("G0v5", out["G0v5"], flush=True)
    # R2v5: K_beta
    betas = [1e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 1e-3, 1e-2]
    r2 = {}
    for nm, b, L in (("M6", bM, LM), ("T6", bT, LT)):
        hj = b.hj_eh(); offs = np.asarray(L["offs"])
        off = np.any(offs != 0, axis=1)
        rec = dict(wrong_sign_taps=int(np.sum(hj[off] > 1e-12 * abs(hj[~off][0]))))
        for n in (64, 128):
            Hk = stencil_grid(hj, offs, n); H0 = float(np.sum(hj))
            for be in betas:
                Kx = np.real(sfft.ifftn(np.exp(-be * (Hk - H0)), workers=2)); s = Kx.sum()
                rec[f"n{n}_b{be:g}"] = dict(neg=float(-Kx[Kx < 0].sum() / s), sum_minus_1=float(s - 1))
        m, _ = mkin_of(b); rec["D_curv"] = 1 / (2 * m)
        n = 128; Hk = stencil_grid(hj, offs, n)
        Kx = np.real(sfft.ifftn(np.exp(-1e-3 * (Hk - np.sum(hj))), workers=2))
        xc = np.fft.fftfreq(n, 1 / n) * b.h
        m2 = np.mean([np.sum(Kx * (xc.reshape([-1 if i == a else 1 for i in range(3)])) ** 2) / Kx.sum() for a in range(3)])
        rec["D_real_beta1e-3"] = float(m2 / 2e-3)
        kk = np.linspace(1e-3, 6e-3, 8) / b.h
        e = b.eps(kk[:, None] * np.array([[1., 0, 0]]))
        c = np.linalg.lstsq(np.stack([kk ** 2, kk ** 4], 1), e, rcond=None)[0]
        rec["excess_kurtosis_beta1e-2"] = float(-6 * c[1] / (1e-2 * c[0] ** 2))
        rec["max_neg"] = max(rec[f"n{n}_b{be:g}"]["neg"] for n in (64, 128) for be in betas)
        rec["max_sum_minus_1"] = max(abs(rec[f"n{n}_b{be:g}"]["sum_minus_1"]) for n in (64, 128) for be in betas)
        rec["PASS"] = bool(rec["wrong_sign_taps"] == 0 and rec["max_neg"] < 1e-10 and rec["max_sum_minus_1"] < 1e-9 and
                           abs(rec["D_real_beta1e-3"] - rec["D_curv"]) < 2e-5 and abs(rec["excess_kurtosis_beta1e-2"]) < 0.02)
        r2[nm] = rec
        print("R2v5", nm, {k: v for k, v in rec.items() if not k.startswith("n")}, flush=True)
    out["R2v5"] = dict(levels=r2, PASS=r2["M6"]["PASS"])
    out["use_for_He"] = bool(out["R3v5"]["PASS"] and out["G0v5"]["PASS"] and out["R2v5"]["PASS"])
    res = load_json(JS); res["checks"] = out; save_json(JS, res)
    print("use_for_He:", out["use_for_He"])

if what == "tables":
    # F(s) = media di eps su p ~ N(0, s I3), con la media sferica della banda (come multi_e_disp.py, solo dentro la zona:
    # la parte oltre la zona si controlla con G4 ed e' trascurabile per l'elio)
    tabs = {}
    D = fib_dirs(400)
    u = np.linspace(0, 14, 7001)[1:]
    wu = np.sqrt(2 / np.pi) * u ** 2 * np.exp(-u ** 2 / 2) * (u[1] - u[0])
    svals = np.geomspace(1e-6, 1e8, 1401)
    for nm, key in (("T6", KT), ("M6", KM)):
        L = level(key); b = Band(L["K"], L["offs"], 16, TH["A"])
        from multi_e_kernel import mkin_of
        m, _ = mkin_of(b)
        kz = np.pi / b.h
        kin = np.unique(np.concatenate([np.geomspace(1e-3, min(200.0, kz), 2001), np.geomspace(min(200.0, kz), kz, 600)]))
        eb = np.concatenate([np.mean(b.eps(kin[i:i + 20, None, None] * D[None]), axis=1) for i in range(0, len(kin), 20)])
        kk = np.linspace(1e-3, 4e-3, 7) / b.h
        c = np.linalg.lstsq(np.stack([kk ** 2, kk ** 4], 1), np.mean(b.eps(kk[:, None, None] * D[None]), axis=1), rcond=None)[0]
        rtab = eb / (kin ** 2 / 2)
        lo = kin < 1.0
        rtab[lo] = 2 * c[0] + 2 * c[1] * kin[lo] ** 2
        F = np.empty_like(svals)
        for i, s in enumerate(svals):
            p = u * np.sqrt(s)
            F[i] = np.sum(wu * np.interp(p, kin, rtab, right=rtab[-1]) * p ** 2 / 2)
        tabs[nm] = dict(s=svals, R=F / (1.5 * svals), m_kin=m, kz=kz)
        print(nm, "m_kin", m, "R(s=1)*m", float(np.interp(0.0, np.log(svals), F / (1.5 * svals))) * m, flush=True)
    np.savez(os.path.join(OUT, "multi_e_v5_Ftab.npz"), **{f"{n}__{k}": v for n, t in tabs.items() for k, v in t.items()})
    res = load_json(JS); res["tables"] = {n: dict(m_kin=t["m_kin"], kz=t["kz"]) for n, t in tabs.items()}; save_json(JS, res)

if what == "he":
    from multi_e_ecg import System, SVM
    ck = res["checks"]
    E = load_json("multi_e_ecg.json"); nr = E["prod"]["nr"]
    K1 = load_json("multi_e_kernel.json"); m_A = K1["R3"]["m_kin_theta_A_used"]
    m0, sm = ck["R3v5"]["m0"], ck["R3v5"]["sigma"]
    out = dict(use_for_He=ck["use_for_He"])
    # M1v5 / E1v5: stessa base, massa nuova (scala esatta E proporzionale a m)
    eHe = nr["he1S_Z2"]["E"] * m0 / m_A; eHp = nr["one_Z2"]["E"] * m0 / m_A
    bar = float(np.hypot(nr["he1S_Z2"]["conv_err"], abs(eHe) * sm))
    out["M1v5"] = dict(E_theta0=eHe, error_bar=bar, diff_mEh=(eHe + 2.903724377) * 1e3, PASS=bool(abs(eHe + 2.903724377) <= 2.9e-4))
    mu = NUC_MASS[2] / (NUC_MASS[2] + 1)
    ie = (eHp - eHe) * mu * EH_EV
    out["E1v5"] = dict(IE_eV=ie, diff_meV=(ie - 24.587389) * 1e3, PASS=bool(abs(ie - 24.587389) <= 0.01),
                       error_bar_eV=float(np.hypot(nr["one_Z2"]["conv_err"], bar) * mu * EH_EV))
    # Q1v5
    FT = np.load(os.path.join(OUT, "multi_e_v5_Ftab.npz"))
    z = np.load(os.path.join(OUT, "multi_e_ecg_he1S_Z2.npz"))
    s_nr = System(2, 1, 2.0, m_A)
    sv = SVM(s_nr, np.random.default_rng(1), 2.0, 20.0, 30); sv.A = z["A"]; E_A = sv.rebuild()
    from multi_e_ecg_prod import expect
    a = ALPHA; c = 1 / a
    u = np.linspace(0, 14, 7001)[1:]
    wu = np.sqrt(2 / np.pi) * u ** 2 * np.exp(-u ** 2 / 2) * (u[1] - u[0])
    q = {}
    for nm in ("T6", "M6"):
        sv_t = FT[f"{nm}__s"]; R = FT[f"{nm}__R"]; m_b = float(FT[f"{nm}__m_kin"])
        s_f = System(2, 1, 2.0, m_A, kin="full", Ftab=(sv_t, R))
        Fd = np.array([np.sum(wu * (np.sqrt((u * np.sqrt(s)) ** 2 * c * c + m_b ** 2 * c ** 4) - m_b * c * c)) for s in sv_t])
        s_d = System(2, 1, 2.0, m_A, kin="full", Ftab=(sv_t, Fd / (1.5 * sv_t)))
        E_nr_b = E_A * m_b / m_A
        d_band = expect(sv, s_f)["E"] - E_nr_b; d_dir = expect(sv, s_d)["E"] - E_nr_b
        sv2 = SVM(s_nr, np.random.default_rng(1), 2.0, 20.0, 30); sv2.A = z["A"]; sv2.rebuild()
        E_var = sv2.rebuild(s_f)
        q[nm] = dict(m_band=m_b, delta_band_first_order_mEh=d_band * 1e3, delta_dirac_first_order_mEh=d_dir * 1e3,
                     band_over_dirac=d_band / d_dir, delta_band_variational_same_basis_mEh=(E_var - E_nr_b) * 1e3,
                     non_dirac_part_mEh=(d_band - d_dir) * 1e3)
        print("Q1v5", nm, q[nm], flush=True)
    old = load_json("multi_e_ecg.json")["prod"]["dirac_ref_posthoc"]["existing_A:he1S_Z2"]
    q["old_level"] = dict(delta_band_mEh=old["delta_band_first_order"] * 1e3, delta_dirac_mEh=old["delta_dirac_first_order"] * 1e3,
                          non_dirac_part_mEh=(old["delta_band_first_order"] - old["delta_dirac_first_order"]) * 1e3)
    v3 = load_json("multi_e_v3.json")["he"]["Q1v3"]
    q["v3_levels"] = {k: dict(non_dirac_part_mEh=v3[k]["non_dirac_part_mEh"], band_over_dirac=v3[k]["band_over_dirac"]) for k in ("T", "M")}
    q["PASS"] = bool(abs(q["T6"]["band_over_dirac"] - 1) <= 0.01 and abs(q["M6"]["band_over_dirac"] - 1) <= 0.02)
    out["Q1v5"] = q
    res = load_json(JS); res["he"] = out; save_json(JS, res)
    print({k: v for k, v in out.items() if k != "Q1v5"})
