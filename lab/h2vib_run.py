"""h2vib_run.py - why is the H2 vibration of the balls -1.4 % low?  Checklist B0-B7 of lab/out/h2vib_prereg.json.
Reads the DMC points from assets/molecole/risultati.txt (no new DMC), re-implements in Python the fit and the
vibrational solver of src/molecole/misura.rs, compares with the exact Born-Oppenheimer curve (Pachucki 2010).
Writes lab/out/h2vib_run.json.  (~2-4 min, 1 thread)"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, sys, time
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from h2vib_pachucki import TABLE

# --- constants exactly as src/molecole/costanti.rs (CODATA 2018)
HBAR, ECH, ME, EPS0, CL = 1.054571817e-34, 1.602176634e-19, 9.1093837015e-31, 8.8541878128e-12, 299792458.0
MP = 1.67262192369e-27
A0 = 4 * np.pi * EPS0 * HBAR**2 / (ME * ECH**2)
EH = HBAR**2 / (ME * A0**2)
H2CM = EH / (2 * np.pi * HBAR * CL) / 100.0
BOHR_A = A0 * 1e10
MP_ME = MP / ME
MU_ATOM = (MP_ME + 1.0) / 2.0     # what misura.rs uses
MU_NUC = MP_ME / 2.0
EXP = dict(we=4401.21, wexe=121.33, nu01=4161.17, re=0.74144)
GRID = [0.9, 1.1, 1.2, 1.25, 1.3, 1.33, 1.4, 1.47, 1.5, 1.55, 1.6, 1.7, 1.9, 2.2, 2.6, 3.2, 4.0, 5.5]


def load_points():
    kv = {}
    for line in open(os.path.join(ROOT, "assets", "molecole", "risultati.txt"), encoding="utf8"):
        if "=" in line:
            k, v = line.split("=", 1)
            try:
                kv[k.strip()] = float(v)
            except ValueError:
                pass
    pts, runs = [], {}
    for r in GRID:
        rr = []
        for sfx in ["", "_run1", "_run2"]:
            k = f"h2_mol_E_R{r:g}{sfx}"
            if k in kv:
                rr.append((kv[k], kv[f"h2_mol_Eerr_R{r:g}{sfx}"]))
        w = sum(1 / d**2 for _, d in rr)
        e = sum(x / d**2 for x, d in rr) / w
        pts.append((r, e, 1 / np.sqrt(w)))
        runs[r] = rr
    return np.array(pts), runs, kv


EXACT = CubicSpline([r for r, _ in TABLE], [e for _, e in TABLE])


# --- the fit of misura.rs
class Curva:
    def __init__(self, e_inf, beta, c):
        self.e_inf, self.beta, self.c = e_inf, beta, np.asarray(c)

    def v(self, r):
        r = np.asarray(r, float)
        k = np.arange(1, len(self.c) + 1)
        return self.e_inf + (self.c[None, :] * np.exp(-np.outer(r, k) * self.beta)).sum(1).reshape(r.shape)

    @staticmethod
    def fit(pts, e_inf, nb=5, free_einf=False):
        r, e, de = pts[:, 0], pts[:, 1], pts[:, 2]
        best = None
        for ib in range(30):
            beta = 0.2 + 0.04 * ib
            A = np.exp(-np.outer(r, np.arange(1, nb + 1)) * beta)
            y = e - e_inf
            if free_einf:
                A = np.hstack([np.ones((len(r), 1)), A]); y = e
            w = 1 / de
            sol = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)[0]
            if free_einf:
                cv = Curva(sol[0], beta, sol[1:])
            else:
                cv = Curva(e_inf, beta, sol)
            chi = float((((cv.v(r) - e) / de) ** 2).sum())
            if best is None or chi < best[1]:
                best = (cv, chi)
        return best

    def minimo(self):
        rr = np.arange(0.5, 8.0, 1e-4)
        vv = self.v(rr)
        i = int(np.argmin(vv))
        re = rr[i]; h = 1e-3
        kk = (self.v(re + h) - 2 * self.v(re) + self.v(re - h)) / h**2
        return float(re), float(vv[i]), float(kk)


def levels(vf, mu, r0=0.4, r1=14.0, npt=6000, nlev=2):
    h = (r1 - r0) / (npt + 1)
    x = r0 + h * np.arange(1, npt + 1)
    d = vf(x) + 1.0 / (mu * h * h)
    off = np.full(npt - 1, -0.5 / (mu * h * h))
    return eigh_tridiagonal(d, off, select="i", select_range=(0, nlev - 1), eigvals_only=True)


def props(vf, mu, minimo=None, npt=6000, r0=0.4, r1=14.0, nlev=2):
    if minimo is None:
        rr = np.arange(0.9, 2.2, 1e-4); vv = vf(rr); i = int(np.argmin(vv)); re = rr[i]; h = 1e-3
        kk = (vf(np.array([re + h]))[0] - 2 * vf(np.array([re]))[0] + vf(np.array([re - h]))[0]) / h**2
        vmin = vv[i]
    else:
        re, vmin, kk = minimo
    we = np.sqrt(kk / mu) * H2CM
    lv = levels(vf, mu, r0, r1, npt, nlev)
    return dict(re_A=re * BOHR_A, we=float(we), nu01=float((lv[1] - lv[0]) * H2CM), zpe=float((lv[0] - vmin) * H2CM),
                levels_cm=((lv - vmin) * H2CM).tolist(), k=float(kk))


def curve_props(cv, mu=MU_ATOM, **kw):
    return props(lambda r: cv.v(r), mu, cv.minimo(), **kw)


def dunham(levels_cm, nv):
    """fit G(v) = we(v+1/2) - wexe(v+1/2)^2 + weye(v+1/2)^3 (+const) to the first nv levels"""
    v = np.arange(nv) + 0.5
    G = np.asarray(levels_cm[:nv])
    A = np.stack([v, -v**2, v**3, np.ones(nv)], 1)
    s = np.linalg.lstsq(A, G, rcond=None)[0]
    return dict(we=float(s[0]), wexe=float(s[1]), weye=float(s[2]))


def local_poly(pts, re0, win, deg):
    m = np.abs(pts[:, 0] - re0) <= win + 1e-9
    p = pts[m]
    if len(p) < deg + 2:
        return None
    x = p[:, 0] - re0
    A = np.stack([x**i for i in range(deg + 1)], 1)
    w = 1 / p[:, 2]
    c = np.linalg.lstsq(A * w[:, None], p[:, 1] * w, rcond=None)[0]
    # minimum of the polynomial near 0 and curvature there
    P = np.polynomial.Polynomial(c)
    dP, d2P = P.deriv(), P.deriv(2)
    roots = [r.real for r in dP.roots() if abs(r.imag) < 1e-9 and abs(r.real) < win]
    x0 = min(roots, key=abs) if roots else 0.0
    kk = d2P(x0)
    return dict(n=int(m.sum()), we=float(np.sqrt(max(kk, 0) / MU_ATOM) * H2CM), re_A=float((re0 + x0) * BOHR_A))


if __name__ == "__main__":
    t0 = time.time()
    res = {}
    pts, runs, kv = load_points()
    e_inf = 2 * kv["h_E_guided"]
    # ---------------- B0 reproduce
    cv, chi = Curva.fit(pts, e_inf)
    p0 = curve_props(cv)
    res["B0"] = dict(python=p0, rust=dict(we=kv["h2_mol_we"], nu01=kv["h2_mol_nu01"], re=kv["h2_mol_re"]),
                     beta=cv.beta, chi2_per_point=chi / len(pts))
    res["B0"]["PASS"] = bool(abs(p0["we"] - kv["h2_mol_we"]) <= 1 and abs(p0["nu01"] - kv["h2_mol_nu01"]) <= 1
                             and abs(p0["re_A"] - kv["h2_mol_re"]) <= 2e-4)
    print("B0", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in p0.items() if k != "levels_cm"},
          "rust", res["B0"]["rust"], "beta", cv.beta, "PASS", res["B0"]["PASS"], flush=True)

    # ---------------- exact BO curve through the same solver
    pe = props(EXACT, MU_ATOM, nlev=8)
    pe_n = props(EXACT, MU_NUC, nlev=8)
    res["exact_BO"] = dict(mu_atom=pe, mu_nuc=pe_n)
    print("exact BO (mu atom):", round(pe["we"], 2), round(pe["nu01"], 2), round(pe["re_A"], 5), flush=True)
    print("exact BO (mu nuc): ", round(pe_n["we"], 2), round(pe_n["nu01"], 2), flush=True)
    # ---------------- B1 harmonic vs anharmonic
    dn = {nv: dunham(pe_n["levels_cm"], nv) for nv in (4, 5, 6)}
    exp_check = EXP["we"] - 2 * EXP["wexe"] + 13 / 4 * 0.8129  # H&H omega_e y_e = 0.8129 (not a free choice: tabulated)
    res["B1"] = dict(harmonic_we_exactBO_mu_nuc=pe_n["we"], harmonic_we_exactBO_mu_atom=pe["we"],
                     dunham_fit_exactBO_mu_nuc=dn, exp_selfconsistency_nu01=exp_check,
                     dev_harmonic_vs_HH=pe_n["we"] / EXP["we"] - 1)
    res["B1"]["CAUSE"] = bool(abs(res["B1"]["dev_harmonic_vs_HH"]) > 0.01)
    # our own curve, Dunham constants
    pc = curve_props(cv, nlev=6)
    res["B1"]["dunham_fit_DMC"] = dunham(pc["levels_cm"], 6)
    print("B1", json.dumps(res["B1"], default=float)[:600], flush=True)
    # ---------------- B2 reduced mass
    pn = curve_props(cv, mu=MU_NUC)
    res["B2"] = dict(we_mu_atom=p0["we"], we_mu_nuc=pn["we"], nu01_mu_atom=p0["nu01"], nu01_mu_nuc=pn["nu01"],
                     rel_change=pn["we"] / p0["we"] - 1)
    res["B2"]["CAUSE"] = bool(abs(res["B2"]["rel_change"]) >= 0.005)
    print("B2", res["B2"], flush=True)
    # ---------------- B3 BO corrections
    res["B3"] = dict(nu01_exactBO_mu_nuc=pe_n["nu01"], nu01_exactBO_mu_atom=pe["nu01"], exp=EXP["nu01"],
                     BO_minus_exp_cm=pe_n["nu01"] - EXP["nu01"], rel=pe_n["nu01"] / EXP["nu01"] - 1,
                     literature="Komasa, Piszczatowski, Lach, Przybytek, Jeziorski, Pachucki, JCTC 7, 3105 (2011): "
                                "adiabatic + nonadiabatic + relativistic + QED shift nu(0->1) by about -2 cm-1 in total")
    res["B3"]["CAUSE"] = bool(abs(res["B3"]["rel"]) >= 0.005)
    print("B3", res["B3"], flush=True)
    # ---------------- B5 model bias on exact data (no noise), and variants on DMC data
    ex_pts = np.array([(r, float(EXACT(r)), de) for r, _, de in pts])
    cve, _ = Curva.fit(ex_pts, -1.0)            # exact E_inf = -1
    pfe = curve_props(cve)
    cve2, _ = Curva.fit(ex_pts, e_inf)           # with our E_inf
    pfe2 = curve_props(cve2)
    res["B5"] = dict(exact_points_fit=dict(we=pfe["we"], nu01=pfe["nu01"], beta=cve.beta),
                     exact_points_fit_ourEinf=dict(we=pfe2["we"], nu01=pfe2["nu01"], beta=cve2.beta),
                     exact_curve=dict(we=pe["we"], nu01=pe["nu01"]),
                     model_bias_we=pfe["we"] / pe["we"] - 1, model_bias_nu01=pfe["nu01"] / pe["nu01"] - 1)
    res["B5"]["CAUSE"] = bool(abs(res["B5"]["model_bias_we"]) >= 0.005)
    var = {}
    cvf, _ = Curva.fit(pts, e_inf, free_einf=True); var["Einf_free"] = curve_props(cvf)
    m = (pts[:, 0] > 0.95) & (pts[:, 0] < 5.0)
    cvc, _ = Curva.fit(pts[m], e_inf); var["no_0.9_5.5"] = curve_props(cvc)
    for nb in (3, 4, 6):
        c_, _ = Curva.fit(pts, e_inf, nb=nb); var[f"nb{nb}"] = curve_props(c_)
    for k in var:
        var[k] = dict(we=var[k]["we"], nu01=var[k]["nu01"], re_A=var[k]["re_A"])
    re0 = cv.minimo()[0]
    loc = {}
    for win in (0.2, 0.3, 0.35, 0.45):
        for deg in (2, 3, 4):
            q = local_poly(pts, re0, win, deg)
            qe = local_poly(ex_pts, re0, win, deg)
            if q:
                loc[f"w{win}_d{deg}"] = dict(dmc=q, exact_points=qe)
    res["B5"]["variants_DMC"] = var
    res["B5"]["local_poly"] = loc
    print("B5", json.dumps({k: v for k, v in res["B5"].items() if k != "local_poly"}, default=float)[:900], flush=True)
    for k, v in loc.items():
        print("  local", k, v["dmc"]["n"], round(v["dmc"]["we"], 1), "exact pts:", round(v["exact_points"]["we"], 1))
    # ---------------- B6 grid
    g = {}
    for npt in (1500, 3000, 6000, 12000):
        g[f"n{npt}"] = curve_props(cv, npt=npt)["nu01"]
    g["range_0.6_8"] = curve_props(cv, r0=0.6, r1=8.0, npt=3000)["nu01"]
    res["B6"] = dict(nu01=g, max_rel_change=max(abs(v / g["n6000"] - 1) for v in g.values()))
    res["B6"]["CAUSE"] = bool(res["B6"]["max_rel_change"] >= 0.001)
    print("B6", res["B6"], flush=True)
    # ---------------- B7 residuals and their curvature
    resid = pts[:, 1] - EXACT(pts[:, 0])
    res["residuals"] = dict(R=pts[:, 0].tolist(), dE=resid.tolist(), err=pts[:, 2].tolist(),
                            chi2_per_point_vs_exact=float(((resid / pts[:, 2]) ** 2).mean()),
                            chi2_per_point_vs_exact_after_mean_shift=float((((resid - np.average(resid, weights=1 / pts[:, 2]**2)) / pts[:, 2]) ** 2).mean()),
                            per_run={str(r): [(e - float(EXACT(r)), d) for e, d in rr] for r, rr in runs.items()})
    m = (pts[:, 0] >= 1.1) & (pts[:, 0] <= 1.9)
    x = pts[m, 0] - 1.4011; w = 1 / pts[m, 2]
    A = np.stack([np.ones_like(x), x, x**2], 1)
    cov = np.linalg.inv((A * w[:, None] ** 2).T @ A)
    c = cov @ ((A * w[:, None] ** 2).T @ resid[m])
    chi_r = float((((A @ c - resid[m]) * w) ** 2).sum() / (m.sum() - 3))
    sc = max(1.0, np.sqrt(chi_r))
    k_res, dk_res = 2 * c[2], 2 * np.sqrt(cov[2, 2]) * sc
    kk_ex = pe["k"]
    shift = np.sqrt((kk_ex + k_res) / kk_ex) - 1
    dshift = 0.5 * dk_res / kk_ex
    res["B7"] = dict(residual_quadratic=c.tolist(), k_residual=float(k_res), k_residual_err=float(dk_res),
                     k_exact=kk_ex, we_shift_rel=float(shift), we_shift_rel_err=float(dshift), chi2_quad=chi_r)
    res["B7"]["CAUSE"] = bool(abs(shift) >= 0.005 and abs(shift) > 2 * dshift)
    print("B7", res["B7"], flush=True)
    print("residuals:", [(r, round(d * 1e4, 2), round(e * 1e4, 2)) for r, d, e in zip(pts[:, 0], resid, pts[:, 2])])
    # ---------------- B4 synthetic noise on exact energies
    rng = np.random.default_rng(20261012)
    NREP = 2000
    syn = []
    for i in range(NREP):
        q = ex_pts.copy()
        q[:, 1] = q[:, 1] + 2.3 * q[:, 2] * rng.standard_normal(len(q))
        c_, _ = Curva.fit(q, -1.0)
        mn = c_.minimo()
        lv = levels(lambda r: c_.v(r), MU_ATOM, npt=3000)
        syn.append((np.sqrt(max(mn[2], 0) / MU_ATOM) * H2CM, (lv[1] - lv[0]) * H2CM))
        if time.time() - t0 > 300:
            break
    syn = np.array(syn)
    def rob(v):
        q16, q50, q84 = np.percentile(v, [16, 50, 84]); return float(q50), float(0.5 * (q84 - q16))
    we_med, we_sd = rob(syn[:, 0]); nu_med, nu_sd = rob(syn[:, 1])
    res["B4"] = dict(nrep=len(syn), we_median=we_med, we_sd=we_sd, nu01_median=nu_med, nu01_sd=nu_sd,
                     we_noiseless=pfe["we"], nu01_noiseless=pfe["nu01"],
                     dmc_we_dev_sigma=(p0["we"] - we_med) / we_sd, dmc_nu01_dev_sigma=(p0["nu01"] - nu_med) / nu_sd,
                     frac_synthetic_below_dmc_we=float((syn[:, 0] <= p0["we"]).mean()),
                     frac_synthetic_below_dmc_nu01=float((syn[:, 1] <= p0["nu01"]).mean()),
                     syn_we=syn[:, 0].tolist(), syn_nu01=syn[:, 1].tolist())
    res["B4"]["fluke"] = bool(abs(res["B4"]["dmc_we_dev_sigma"]) <= 2.5 and res["residuals"]["chi2_per_point_vs_exact_after_mean_shift"] <= 4)
    print("B4", {k: v for k, v in res["B4"].items() if not k.startswith("syn")}, flush=True)
    res["time_s"] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, "h2vib_run.json"), "w"), indent=1, default=float)
    print("done", round(res["time_s"]), "s")
