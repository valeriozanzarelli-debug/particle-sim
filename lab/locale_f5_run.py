"""locale_f5_run.py - phase 5: the breathing pair (electron breathing coin + light breathing between F+ and F-),
exams E1-E6 of lab/out/locale_f5_prereg.json. Writes lab/out/locale_f5_run.json (exit 3 = run again)."""
import os, sys, json, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locale_core import SplitWalk, sdot, zigzag_ticks, n_for_tol, dirs_for, RealWalk, PI
from locale_f4_core import walk as ewalk, measure, E_low, DIRS

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "locale_f5_run.json")
T0 = time.time()
R = json.load(open(OUT)) if os.path.exists(OUT) else {}
X6 = np.kron(np.array([[0, 1], [1, 0]], complex), np.eye(3))


def lwalk(mu, T_phys=2.0):
    fac = []
    for it, tk in enumerate(zigzag_ticks(3)):
        s = -1 if it == 0 else +1
        for d in tk:
            S = sdot(d, 1)
            M = np.zeros((6, 6), complex); M[:3, :3] = S; M[3:, 3:] = -S
            fac.append(("shift", d, M))
        if mu != 0:
            fac.append(("coin", np.cos(s * mu) * np.eye(6) - 1j * np.sin(s * mu) * X6))
    return SplitWalk(fac, 6, np.eye(3), T_phys=T_phys, dim=3)


def c_light(w, k=1e-3):
    out = []
    for d in DIRS:
        eps = np.sort(np.abs(w.phases((k * d)[None])[0][0]))
        out.append(eps[-1] / (k * w.T))
    return float(np.mean(out)), float(np.max(out) - np.min(out))


def c_electron(w, th):
    cs = []
    for d in DIRS:
        ps = np.linspace(0.02, 0.1, 9) * th
        E = E_low(w, ps[:, None] * d[None]); E0 = E_low(w, np.zeros((1, 3)))[0]
        cc = np.linalg.lstsq(np.stack([ps ** 2, ps ** 4], 1), E ** 2 - E0 ** 2, rcond=None)[0]
        cs.append(np.sqrt(cc[0]))
    return float(np.mean(cs)), float(np.max(cs) - np.min(cs))


def helicity(w, k=2 * PI / 20):
    K = np.array([[0, 0, k]])
    eps, v = w.phases(K)
    j = np.argmax(eps[0])                         # the + band (largest phase, = +c k T)
    vec = v[0, :, j]
    return float(np.sum(np.abs(vec[:3]) ** 2) / np.sum(np.abs(vec) ** 2))


def norm_check(w, nint, n=24, per=50):
    rng = np.random.default_rng(1)
    psi = rng.normal(size=(nint, n, n, n)) + 1j * rng.normal(size=(nint, n, n, n))
    psi /= np.linalg.norm(psi)
    rw = RealWalk(w, n, 3)
    for _ in range(per):
        psi = rw.period(psi)
    return float(abs(np.linalg.norm(psi) - 1))


def todo(k):
    if k in R:
        return False
    if time.time() - T0 > 480:
        json.dump(R, open(OUT, "w"), indent=1); print("rerun"); sys.exit(3)
    print(f"[{time.time() - T0:.0f}s] {k}", flush=True)
    return True


# light: speed, N_1%, helicity, unitarity, norm
for mu in (0.0, 0.3, 0.5):
    key = f"light|{mu}"
    if not todo(key):
        continue
    w = lwalk(mu)
    cl, cl_spread = c_light(w)
    wr = lwalk(mu, T_phys=2.0 * cl)              # time unit rescaled so that the measured light speed is 1
    N1 = n_for_tol(wr, dirs_for(3, 80))
    U = w.bloch(np.random.default_rng(0).uniform(-PI, PI, (50, 3)))
    R[key] = dict(c_l=cl, c_l_spread=cl_spread, cos_mu=float(np.cos(mu)), N1=float(N1),
                  omega0=float(np.max(np.abs(w.phases(np.zeros((1, 3)))[0]))),
                  unitarity=float(np.max(np.abs(np.einsum("nba,nbc->nac", U.conj(), U) - np.eye(6)))),
                  norm_drift=norm_check(w, 6), helicity_plus=helicity(w))
    print(key, R[key], flush=True)
    json.dump(R, open(OUT, "w"), indent=1)

# electron: speed, gaps, E0, unitarity, norm
for th in (0.1, 0.02):
    for mu in (0.0, 0.3, 0.5):
        for hu in (0.0, 0.1):
            key = f"electron|{th}|{mu}|{hu}"
            if not todo(key):
                continue
            w = ewalk(th, mu_x=mu, hu=hu)
            r = measure(w, th)
            ce, ce_spread = c_electron(w, th)
            r.update(c_e=ce, c_e_spread=ce_spread, norm_drift=norm_check(w, 4) if th == 0.1 else None)
            R[key] = r
            print(key, f"c_e {ce:.6f} gap {r['min_doubler_gap']:.3f} pi {min(r['X_x'], r['R_xyz'], r['zone_gpi']):.3f} "
                  f"E0 {r['E0_rel']:.0e}", flush=True)
            json.dump(R, open(OUT, "w"), indent=1)

# verdicts
V = {}
for th in (0.1, 0.02):
    for mu in (0.3, 0.5):
        for hu in (0.0, 0.1):
            e = R[f"electron|{th}|{mu}|{hu}"]; e0 = R[f"electron|{th}|0.0|0.0"]
            l = R[f"light|{mu}"]; l0 = R["light|0.0"]
            mis = e["c_e"] / l["c_l"] - 1; mis0 = e0["c_e"] / l0["c_l"] - 1
            E1 = abs(mis) <= 1e-4 if th == 0.02 else abs(mis - mis0) <= 1e-4
            pi_gap = min(e["X_x"], e["X_y"], e["X_z"], e["R_xyz"], e["zone_gpi"])
            E2 = pi_gap >= 0.3 and (hu == 0 or min(e["N_xy"], e["zone_g0"]) >= 0.3)
            mkin = (l["c_l"] / e["c_e"]) ** 2 - 1; mkin0 = (l0["c_l"] / e0["c_e"]) ** 2 - 1
            E3 = e["E0_rel"] <= 1e-12 and abs(mkin - mkin0) <= 1e-4
            E4 = l["N1"] <= 1.1 * l0["N1"] and l["omega0"] <= 1e-12
            E5 = l["unitarity"] <= 1e-12 and e["unitarity"] <= 1e-12 and l["norm_drift"] <= 1e-12 and \
                (e["norm_drift"] is None or e["norm_drift"] <= 1e-12)
            V[f"{th}|{mu}|{hu}"] = dict(mismatch=mis, mismatch_frozen=mis0, pi_gap=pi_gap, all_doublers=e["min_doubler_gap"],
                                        mkin_minus_frozen=mkin - mkin0, E1=E1, E2=E2, E3=E3, E4=E4, E5=E5,
                                        ADOPTABLE=bool(E1 and E2 and E3 and E4 and E5))
            print(th, mu, hu, V[f"{th}|{mu}|{hu}"], flush=True)
R["verdicts"] = V
json.dump(R, open(OUT, "w"), indent=1)
print("ALL DONE", time.time() - T0)
