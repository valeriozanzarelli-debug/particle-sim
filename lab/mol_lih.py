"""mol_lih.py - study mol_ (docs/MOLECOLE.md), tier "esatto": LiH with all FOUR electrons (ECG, mol_ecg.py).
Pre-registration: lab/out/mol_prereg.json (lab/mol_prereg.py).

Jobs (each <= ~9 minutes, checkpoints in lab/out/mol_lih_*.npz, not committed; results in lab/out/mol_lih.json):
  li       Li atom (3 e, doublet) with the same machine, grown to K_LI and refined
  re       LiH at R0 = 3.015 a0 (the starting geometry; the bond length is then FOUND from the curve), grown to K_RE
  curve    every R of RGRID: start from the R0 basis with centres stretched (z -> z R/R0, Li stays at 0, H at R),
           NSWEEP refinement sweeps; energy and dipole moment
  analyze  R_e (minimum of the curve), D_e = E(Li) + E(H) - E(R_e), omega_e, vibrational levels (1-D, nuclear masses),
           dipole moment at R_e and averaged over v = 0; verdicts against the pre-registration
H atom: E(H) = -m_kin/2 exactly for this operator (closed form, not computed)."""
import os, sys, json, time
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mol_ecg as M

OUT = M.OUT
JS = os.path.join(OUT, "mol_lih.json")
PRE = json.load(open(os.path.join(OUT, "mol_prereg.json")))["LiH"]
CFG = PRE["config"]
R0 = CFG["R0"]
BUDGET = 530.0


def load():
    return json.load(open(JS)) if os.path.exists(JS) else {}


def save(res):
    json.dump(res, open(JS + ".tmp", "w"), indent=1, default=float)
    os.replace(JS + ".tmp", JS)


def groups_lih(R):
    kc = (0, 0.0, -0.3, 2.6)            # Li 1s pair: tight, on Li
    kb = (0.5 * R, 0.8, -1.8, 0.2)      # bond region: diffuse, between the nuclei
    kh = (R, 0.2, -1.3, 1.0)            # on H
    kl = (0, 0.5, -1.8, 0.3)            # Li valence
    return [(1, [kc, kc, kb, kb]), (1, [kc, kc, kh, kb]), (1, [kc, kc, kh, kh]), (1, [kc, kc, kl, kh]),
            (0.3, [kc, kc, kl, kl]), (0.3, [kc, kh, kc, kh])]


def groups_li():
    kc = (0, 0.0, -0.3, 2.6); kv = (0, 0.0, -1.8, 0.6)
    return [(1, [kc, kc, kv]), (0.3, [kc, kv, kv]), (0.3, [kc, kc, kc]), (0.2, [kv, kc, kc])]


def sys_lih(R):
    return M.System(4, [(3, 0.0), (1, R)], M.spin_singlet_pairs(2))


def sys_li():
    return M.System(3, [(3, 0.0)], M.spin_singlet_pairs(1, 1))


def ck_path(name):
    return os.path.join(OUT, f"mol_lih_{name}.npz")


def ck_save(name, sv, stage):
    np.savez(ck_path(name), A=sv.A, s=sv.s, stage=json.dumps(stage),
             rng=json.dumps(sv.rng.bit_generator.state))


def ck_load(name, sv):
    p = ck_path(name)
    if not os.path.exists(p):
        return None
    z = np.load(p, allow_pickle=False)
    sv.set(z["A"], z["s"])
    sv.rng.bit_generator.state = json.loads(str(z["rng"]))
    return json.loads(str(z["stage"]))


def grow_refine(name, sv, K, nsweep, ntg, ntr, t0, log):
    """resumable: grow to K (ntg trials), then nsweep sweeps (ntr trials); returns True when finished"""
    st = ck_load(name, sv) or dict(sweep=0, k=0)
    while len(sv.A) < K:
        for _ in range(10):
            if len(sv.A) < K:
                sv.grow(ntg)
        if time.time() - t0 > BUDGET:
            ck_save(name, sv, st); log(f"{name}: K={len(sv.A)} E={sv.E:.7f} checkpoint"); return False
    while st["sweep"] < nsweep:
        while st["k"] < len(sv.A):
            sv.refine_one(st["k"], ntr); st["k"] += 1
            if time.time() - t0 > BUDGET:
                ck_save(name, sv, st); log(f"{name}: sweep {st['sweep']} k={st['k']} E={sv.E:.7f} checkpoint"); return False
        st["sweep"] += 1; st["k"] = 0
        log(f"{name}: sweep {st['sweep']} done E={sv.E:.8f}")
    ck_save(name, sv, st)
    return True


def dipole_and_energy(sy, A, s):
    """<sum z_i>, E with the eigenvector of the full (numpy) matrices; also returns the coefficients"""
    import scipy.linalg as sla
    S, H, D = sy.mats(A, s, dip=True)
    w, U = sla.eigh(H, S, subset_by_index=[0, 0])
    c = U[:, 0] / np.sqrt(U[:, 0] @ S @ U[:, 0])
    return float(w[0]), float(c @ D @ c), c


def job_li(log=print):
    t0 = time.time()
    sv = M.SVM(sys_li(), M.make_gen(3, groups_li()), np.random.default_rng(CFG["seed_li"]))
    done = grow_refine("li", sv, CFG["K_LI"], CFG["nsweep_li"], CFG["ntrial_grow"], CFG["ntrial_ref"], t0, log)
    if done:
        res = load(); res["Li"] = dict(E=sv.E, K=len(sv.A)); save(res); log(f"LI DONE E={sv.E:.8f}")


def job_re(log=print):
    t0 = time.time()
    sv = M.SVM(sys_lih(R0), M.make_gen(4, groups_lih(R0)), np.random.default_rng(CFG["seed_re"]))
    done = grow_refine("re", sv, CFG["K_RE"], CFG["nsweep_re"], CFG["ntrial_grow"], CFG["ntrial_ref"], t0, log)
    if done:
        res = load(); res["R0"] = dict(E=sv.E, K=len(sv.A)); save(res); log(f"RE DONE E={sv.E:.8f}")


def job_curve(log=print):
    t0 = time.time()
    z = np.load(ck_path("re"))
    Are, sre = z["A"], z["s"]
    res = load(); cur = res.setdefault("curve", {})
    for R in CFG["RGRID"]:
        key = f"{R:.3f}"
        if key in cur:
            continue
        name = f"R{key}"
        sy = sys_lih(R)
        sv = M.SVM(sy, M.make_gen(4, groups_lih(R)), np.random.default_rng(CFG["seed_curve"] + int(1000 * R)))
        if not os.path.exists(ck_path(name)):
            s = sre.copy(); s[..., 2] *= R / R0
            sv.set(Are, s)
            ck_save(name, sv, dict(sweep=0, k=0))
            log(f"{name}: start E={sv.E:.7f}")
        done = grow_refine(name, sv, len(Are), CFG["nsweep_curve"], CFG["ntrial_grow"], CFG["ntrial_ref"], t0, log)
        if not done:
            return
        E, zsum, _ = dipole_and_energy(sy, sv.A, sv.s)
        mu = 1.0 * R - zsum                     # nuclear charges (Li at 0, H at R) minus electrons, a.u.; neutral
        res = load(); res.setdefault("curve", {})[key] = dict(R=R, E=sv.E, E_numpy=E, mu_au=mu, K=len(sv.A))
        save(res); cur = res["curve"]
        log(f"{name}: E={sv.E:.8f} mu={mu:.5f} au")
        if time.time() - t0 > BUDGET - 60:
            log("checkpoint: run again"); return
    log("CURVE DONE")


def job_analyze(log=print):
    from scipy.interpolate import CubicSpline
    res = load()
    C = res["curve"]
    Rs = np.array(sorted(v["R"] for v in C.values()))
    Es = np.array([C[f"{R:.3f}"]["E"] for R in Rs]); mus = np.array([C[f"{R:.3f}"]["mu_au"] for R in Rs])
    cs = CubicSpline(Rs, Es)
    x = np.linspace(Rs[0], Rs[-1], 200001); i = int(np.argmin(cs(x)))
    Re_spl, Emin_spl = float(x[i]), float(cs(x[i]))
    # pre-registered omega_e: polynomial fit of degree 4 in the window
    lo, hi = CFG["fit_window"]
    sel = (Rs >= lo - 1e-9) & (Rs <= hi + 1e-9)
    p = np.polyfit(Rs[sel] - R0, Es[sel], 4)
    xx = np.linspace(lo, hi, 20001); yy = np.polyval(p, xx - R0); j = int(np.argmin(yy))
    Re, Emin = float(xx[j]), float(yy[j])
    k2 = float(np.polyval(np.polyder(p, 2), Re - R0))
    u = 1822.888486209                               # m_e per dalton (constant)
    mLi = 7.016003434 * u - 3.0                      # 7Li nucleus (atomic mass minus 3 electrons; binding < 1e-7)
    mH = 1836.15267343                               # proton
    mu_red = mLi * mH / (mLi + mH)
    EH_EV = 27.211386245988; HART_CM = 219474.6313632
    we = float(np.sqrt(k2 / mu_red) * HART_CM)
    EH = -0.5 * M.M_KIN
    ELi = res["Li"]["E"]
    De = (ELi + EH - Emin) * EH_EV
    # vibrational levels: sinc-DVR on the spline (R in [Rs0, Rsmax])
    n = 1200; g = np.linspace(Rs[0], Rs[-1], n); h = g[1] - g[0]
    ii = np.arange(n)
    Tm = np.where(ii[:, None] == ii[None], np.pi ** 2 / 3, 2.0 / np.where(ii[:, None] == ii[None], 1, (ii[:, None] - ii[None]) ** 2))
    Tm = Tm * (-1.0) ** (ii[:, None] - ii[None]) / (2 * mu_red * h * h)
    Hm = Tm + np.diag(cs(g))
    w, U = np.linalg.eigh(Hm)
    nu01 = float((w[1] - w[0]) * HART_CM)
    D0 = float((ELi + EH - w[0]) * EH_EV)
    csm = CubicSpline(Rs, mus)
    mu0 = float(np.sum(U[:, 0] ** 2 * csm(g)))
    AU_D = 2.541746473
    mu_e = float(abs(csm(Re)) * AU_D); mu_v0 = float(abs(mu0) * AU_D)
    T = PRE["exams"]
    out = dict(R_e=Re, R_e_spline=Re_spl, E_min=Emin, E_min_spline=Emin_spl, E_Li=ELi, E_H=EH, D_e_eV=De,
               omega_e_cm=we, nu01_cm=nu01, D0_eV=D0, mu_e_D=mu_e, mu_v0_D=mu_v0, mu_red_me=mu_red,
               levels_cm=((w[:6] - w[0]) * HART_CM).tolist())
    E_R0 = C.get(f"{R0:.3f}", {}).get("E")
    v = {}
    v["L0"] = dict(E_R0=E_R0, ref=T["L0"]["ref"], err_mEh=1e3 * (E_R0 - T["L0"]["ref"]),
                   PASS=bool(abs(E_R0 - T["L0"]["ref"]) <= T["L0"]["tol_Eh"]))
    v["L1"] = dict(val=Re, ref=T["L1"]["ref"], rel=Re / T["L1"]["ref"] - 1, PASS=bool(abs(Re / T["L1"]["ref"] - 1) <= T["L1"]["tol_rel"]))
    v["L2"] = dict(val=De, ref=T["L2"]["ref"], diff=De - T["L2"]["ref"], PASS=bool(abs(De - T["L2"]["ref"]) <= T["L2"]["tol_eV"]))
    v["L3"] = dict(val=we, ref=T["L3"]["ref"], rel=we / T["L3"]["ref"] - 1, PASS=bool(abs(we / T["L3"]["ref"] - 1) <= T["L3"]["tol_rel"]))
    v["L4"] = dict(val=mu_v0, ref=T["L4"]["ref"], rel=mu_v0 / T["L4"]["ref"] - 1, PASS=bool(abs(mu_v0 / T["L4"]["ref"] - 1) <= T["L4"]["tol_rel"]))
    v["L5"] = dict(val=ELi, ref=T["L5"]["ref"], err_mEh=1e3 * (ELi - T["L5"]["ref"]), PASS=bool(abs(ELi - T["L5"]["ref"]) <= T["L5"]["tol_Eh"]))
    res["analysis"] = out; res["verdict"] = v
    save(res)
    log(json.dumps(out, indent=1)); log(json.dumps(v, indent=1, default=float))


if __name__ == "__main__":
    job = sys.argv[1]
    {"li": job_li, "re": job_re, "curve": job_curve, "analyze": job_analyze}[job]()
