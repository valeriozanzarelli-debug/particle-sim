"""multi_e_v3.py - corse della pre-registrazione v3 (lab/out/multi_e_prereg_v3.json): livello 1 T/M di AGGLOMERATO_COMPLETO.
  learn   impara i livelli con aggfull_l1.prod_one (importato, non modificato) -> lab/out/multi_e_v3_<key>.pkl, multi_e_v3.json
  checks  R3v3 (m0 dal nostro mkin_of), L1c_v3 (controllo Z = 10 di auto_level), G0v3 (parita'), R2v3 (K_beta del livello M)
  tables  tabelle F(s) delle bande T e M (stesso algoritmo di multi_e_disp.py)
  he      M1v3, E1v3, Q1v3 sull'elio (stato ECG NR della v2)
Run: python lab/multi_e_v3.py learn | checks | tables | he"""
import os, sys, time, pickle, itertools
import numpy as np
import scipy.fft as sfft
from multi_e_core import OUT, ALPHA, A0, THETA, EH_EV, NUC_MASS, load_json, save_json, Band, fib_dirs, dirs26

JS = "multi_e_v3.json"
res = load_json(JS)
what = sys.argv[1] if len(sys.argv) > 1 else "learn"
TH = {"A": THETA["A"], "B": THETA["B"], "A2": THETA["A2"]}
JOBS = [("A", 10, 4, "T", s) for s in (1, 2, 3)] + [("A", 2, 4, "M", 1)] + \
       [(c, 1, 2, "T", s) for c in ("A2", "A", "B") for s in (1, 2)]


def pkl(key):
    return os.path.join(OUT, f"multi_e_v3_{key}.pkl")


def level(key):
    return pickle.load(open(pkl(key), "rb"))


if what == "learn":
    import aggfull_l1 as al
    lr = res.get("learn", {}); fine = res.get("fine", {})
    t0 = time.time()
    for chain, obs, w, m, seed in JOBS:
        key = f"{chain}_Z{obs}_w{w}_{m}_s{seed}"
        if key in lr:
            continue
        t1 = time.time()
        out, r = al.prod_one(chain, obs, w, m, seed, fine)
        pickle.dump(dict(K=r.K, offs=r.offs, phi0=r.phi0, B=16, b=16, w=w, C=2, desc="pauli2", chain=chain, learner=m),
                    open(pkl(key), "wb"))
        out["sec"] = time.time() - t1
        lr[key] = out
        res = load_json(JS); res["learn"] = lr; res["fine"] = fine; save_json(JS, res)
        print(key, f"{out['sec']:.0f}s", "check_Z10", out.get("check_Z10", {}).get("energy_rel"), flush=True)
        if time.time() - t0 > 480:
            print("checkpoint: run again", flush=True); sys.exit(0)
    print("learn done", flush=True)

if what == "checks":
    from multi_e_kernel import mkin_of, stencil_grid
    out = {}
    # R3v3
    mk = {}
    for c in ("A2", "A", "B"):
        for s in (1, 2):
            L = level(f"{c}_Z1_w2_T_s{s}")
            mk[f"{c}_s{s}"] = mkin_of(Band(L["K"], L["offs"], 16, TH[c], sym=True))[0]
    fits = []
    for sa, sb, sc in itertools.product((1, 2), repeat=3):
        th = np.array([TH["A2"], TH["A"], TH["B"]])
        mm = np.array([mk[f"A2_s{sa}"], mk[f"A_s{sb}"], mk[f"B_s{sc}"]])
        fits.append(np.linalg.solve(np.stack([np.ones(3), th ** 2, th ** 4], 1), mm))
    fits = np.array(fits); m0 = float(fits[:, 0].mean()); sm = float(fits[:, 0].std(ddof=1))
    LT = level("A_Z10_w4_T_s1"); LM = level("A_Z2_w4_M_s1")
    bT = Band(LT["K"], LT["offs"], 16, TH["A"]); bM = Band(LM["K"], LM["offs"], 16, TH["A"])
    mT, isoT = mkin_of(bT); mM, isoM = mkin_of(bM)
    out["R3v3"] = dict(m_kin=mk, m0=m0, sigma=sm, a=float(fits[:, 1].mean()), PASS=bool(abs(m0 - 1) <= 5e-5),
                       within_1e9=bool(abs(m0 - 1) <= 1e-9), m_T_thetaA=mT, m_M_thetaA=mM, tan_law_A=float(np.tan(TH["A"]) / TH["A"]))
    print("R3v3", out["R3v3"], flush=True)
    # L1c_v3
    lr = res["learn"]
    ch = {s: lr[f"A_Z10_w4_T_s{s}"].get("check_Z10", {}) for s in (1, 2, 3)}
    out["L1c_v3"] = dict(energy_rel={s: v.get("energy_rel") for s, v in ch.items()},
                         PASS=bool(all(v.get("ok_budget_3e-5") for v in ch.values())))
    print("L1c_v3", out["L1c_v3"], flush=True)
    # G0v3: parita'
    g0 = {}
    rng = np.random.default_rng(3)
    D = fib_dirs(200)
    for nm, L in (("T", LT), ("M", LM)):
        raw = Band(L["K"], L["offs"], 16, TH["A"], sym=False)
        for kh in (0.01, 1.0):
            kap = kh * D
            gp = raw.gen(kap) - raw.gen(np.zeros(3)); gm = raw.gen(-kap) - raw.gen(np.zeros(3))
            g0[f"{nm}_kh{kh:g}"] = float(np.max(np.abs(gp - gm)) / np.max(np.abs(gp + gm)))
        # parte immaginaria dello stencil (parte dispari in k)
        g0[f"{nm}_imag_over_real"] = float(np.max(np.abs(np.imag(0.5 * (L["K"][0, 0] + L["K"][1, 1])))) /
                                           np.max(np.abs(np.real(0.5 * (L["K"][0, 0] + L["K"][1, 1])))))
    out["G0v3"] = dict(odd_rel=g0, PASS=bool(all(v < 1e-10 for k, v in g0.items() if "_kh" in k)))
    print("G0v3", out["G0v3"], flush=True)
    # R2v3: K_beta
    betas = [1e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 1e-3, 1e-2]
    r2 = {}
    for nm, b, L in (("M", bM, LM), ("T", bT, LT)):
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
        rec["PASS"] = bool(rec["max_neg"] < 1e-10 and rec["max_sum_minus_1"] < 1e-9 and
                           abs(rec["D_real_beta1e-3"] - rec["D_curv"]) < 2e-5 and abs(rec["excess_kurtosis_beta1e-2"]) < 0.02)
        r2[nm] = rec
        print("R2v3", nm, {k: v for k, v in rec.items() if not k.startswith("n")}, flush=True)
    out["R2v3"] = dict(levels=r2, PASS=r2["M"]["PASS"])
    out["use_for_He"] = bool(out["R3v3"]["PASS"] and out["G0v3"]["PASS"] and out["R2v3"]["PASS"])
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
    for nm, key in (("T", "A_Z10_w4_T_s1"), ("M", "A_Z2_w4_M_s1")):
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
    np.savez(os.path.join(OUT, "multi_e_v3_Ftab.npz"), **{f"{n}__{k}": v for n, t in tabs.items() for k, v in t.items()})
    res = load_json(JS); res["tables"] = {n: dict(m_kin=t["m_kin"], kz=t["kz"]) for n, t in tabs.items()}; save_json(JS, res)

if what == "he":
    from multi_e_ecg import System, SVM
    ck = res["checks"]
    E = load_json("multi_e_ecg.json"); nr = E["prod"]["nr"]
    K1 = load_json("multi_e_kernel.json"); m_A = K1["R3"]["m_kin_theta_A_used"]
    m0, sm = ck["R3v3"]["m0"], ck["R3v3"]["sigma"]
    out = dict(use_for_He=ck["use_for_He"])
    # M1v3 / E1v3: stessa base, massa nuova (scala esatta E proporzionale a m)
    eHe = nr["he1S_Z2"]["E"] * m0 / m_A; eHp = nr["one_Z2"]["E"] * m0 / m_A
    bar = float(np.hypot(nr["he1S_Z2"]["conv_err"], abs(eHe) * sm))
    out["M1v3"] = dict(E_theta0=eHe, error_bar=bar, diff_mEh=(eHe + 2.903724377) * 1e3, PASS=bool(abs(eHe + 2.903724377) <= 2.9e-4))
    mu = NUC_MASS[2] / (NUC_MASS[2] + 1)
    ie = (eHp - eHe) * mu * EH_EV
    out["E1v3"] = dict(IE_eV=ie, diff_meV=(ie - 24.587389) * 1e3, PASS=bool(abs(ie - 24.587389) <= 0.01),
                       error_bar_eV=float(np.hypot(nr["one_Z2"]["conv_err"], bar) * mu * EH_EV))
    # Q1v3
    FT = np.load(os.path.join(OUT, "multi_e_v3_Ftab.npz"))
    z = np.load(os.path.join(OUT, "multi_e_ecg_he1S_Z2.npz"))
    s_nr = System(2, 1, 2.0, m_A)
    sv = SVM(s_nr, np.random.default_rng(1), 2.0, 20.0, 30); sv.A = z["A"]; E_A = sv.rebuild()
    from multi_e_ecg_prod import expect
    a = ALPHA; c = 1 / a
    u = np.linspace(0, 14, 7001)[1:]
    wu = np.sqrt(2 / np.pi) * u ** 2 * np.exp(-u ** 2 / 2) * (u[1] - u[0])
    q = {}
    for nm in ("T", "M"):
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
        print("Q1v3", nm, q[nm], flush=True)
    old = load_json("multi_e_ecg.json")["prod"]["dirac_ref_posthoc"]["existing_A:he1S_Z2"]
    q["old_level"] = dict(delta_band_mEh=old["delta_band_first_order"] * 1e3, delta_dirac_mEh=old["delta_dirac_first_order"] * 1e3,
                          non_dirac_part_mEh=(old["delta_band_first_order"] - old["delta_dirac_first_order"]) * 1e3)
    q["PASS"] = bool(abs(q["T"]["band_over_dirac"] - 1) <= 0.01 and abs(q["M"]["band_over_dirac"] - 1) <= 0.02)
    out["Q1v3"] = q
    res = load_json(JS); res["he"] = out; save_json(JS, res)
    print({k: v for k, v in out.items() if k != "Q1v3"})
