"""h2vib_fix.py - the corrected H2 vibration ("corrected_number" of lab/out/h2vib_prereg.json).
B5 found the cause: the fit model of misura.rs (5 exponentials, E_inf fixed, beta on a 0.04 grid), fed with the
EXACT energies at our 18 R, already gives omega_e -1.6 % and nu01 -0.86 %. Here: candidate models are judged ONLY
on exact data (bias, no noise) and on exact data + our noise (spread); the chosen one is the unbiased model with the
smallest spread; then it is applied to the DMC points. Writes lab/out/h2vib_fix.json. (~4 min, 1 thread)"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from h2vib_run import (Curva, load_points, curve_props, props, levels, EXACT, MU_ATOM, MU_NUC, H2CM, EXP, OUT,
                       local_poly)


class CurvaB(Curva):
    """same model, beta optimised continuously (golden section on the chi2 after the 0.04 scan)"""
    @staticmethod
    def fit(pts, e_inf, nb=5, free_einf=False, bmin=0.2, bmax=1.4):
        def chi_at(beta):
            r, e, de = pts[:, 0], pts[:, 1], pts[:, 2]
            A = np.exp(-np.outer(r, np.arange(1, nb + 1)) * beta)
            y = e - e_inf
            if free_einf:
                A = np.hstack([np.ones((len(r), 1)), A]); y = e
            w = 1 / de
            sol = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)[0]
            cv = Curva(sol[0], beta, sol[1:]) if free_einf else Curva(e_inf, beta, sol)
            return float((((cv.v(r) - e) / de) ** 2).sum()), cv
        bs = np.arange(bmin, bmax, 0.01)
        ch = [chi_at(b)[0] for b in bs]
        i = int(np.argmin(ch))
        a, b = bs[max(i - 1, 0)], bs[min(i + 1, len(bs) - 1)]
        g = (np.sqrt(5) - 1) / 2
        for _ in range(40):
            c1, c2 = b - g * (b - a), a + g * (b - a)
            if chi_at(c1)[0] < chi_at(c2)[0]:
                b = c2
            else:
                a = c1
        chi, cv = chi_at(0.5 * (a + b))
        return cv, chi


MODELS = {
    "rust_nb5_grid": lambda q, ei: Curva.fit(q, ei)[0],
    "nb5_beta_cont": lambda q, ei: CurvaB.fit(q, ei)[0],
    "nb6_grid": lambda q, ei: Curva.fit(q, ei, nb=6)[0],
    "nb6_beta_cont": lambda q, ei: CurvaB.fit(q, ei, nb=6)[0],
    "nb5_Einf_free": lambda q, ei: Curva.fit(q, ei, free_einf=True)[0],
    "nb5_Einf_free_cont": lambda q, ei: CurvaB.fit(q, ei, free_einf=True)[0],
}


def wn(cv):
    mn = cv.minimo()
    lv = levels(lambda r: cv.v(r), MU_ATOM, npt=3000)
    return np.sqrt(max(mn[2], 0) / MU_ATOM) * H2CM, (lv[1] - lv[0]) * H2CM, mn[0]


if __name__ == "__main__":
    t0 = time.time()
    pts, runs, kv = load_points()
    e_inf = 2 * kv["h_E_guided"]
    ex = np.array([(r, float(EXACT(r)), de) for r, _, de in pts])
    pe = props(EXACT, MU_ATOM)
    res = dict(exact_curve=dict(we=pe["we"], nu01=pe["nu01"]), models={})
    rng = np.random.default_rng(777)
    NREP = 150
    noise = [2.3 * ex[:, 2] * rng.standard_normal(len(ex)) for _ in range(NREP)]
    for name, f in MODELS.items():
        we0, nu0, _ = wn(f(ex, -1.0))
        syn = []
        for z in noise:
            q = ex.copy(); q[:, 1] += z
            syn.append(wn(f(q, -1.0))[:2])
        syn = np.array(syn)
        sd = 0.5 * (np.percentile(syn, 84, 0) - np.percentile(syn, 16, 0))
        med = np.median(syn, 0)
        cvd = f(pts, e_inf)
        p = curve_props(cvd)
        res["models"][name] = dict(bias_we=we0 / pe["we"] - 1, bias_nu01=nu0 / pe["nu01"] - 1,
                                   noisy_median_bias_we=med[0] / pe["we"] - 1, noisy_median_bias_nu01=med[1] / pe["nu01"] - 1,
                                   sd_we=float(sd[0]), sd_nu01=float(sd[1]),
                                   dmc=dict(we=p["we"], nu01=p["nu01"], re_A=p["re_A"], beta=cvd.beta, zpe=p["zpe"]))
        print(name, {k: (round(v, 5) if isinstance(v, float) else v) for k, v in res["models"][name].items() if k != "dmc"},
              "DMC:", round(p["we"], 1), round(p["nu01"], 1), round(p["re_A"], 5), f"{time.time()-t0:.0f}s", flush=True)
    # choice (rule fixed in the docstring, applied mechanically): |noiseless bias| <= 0.3 % AND |noisy-median bias|
    # <= 0.3 % for both omega_e and nu01; among those the smallest sd_we
    ok = {k: v for k, v in res["models"].items()
          if max(abs(v["bias_we"]), abs(v["bias_nu01"]), abs(v["noisy_median_bias_we"]), abs(v["noisy_median_bias_nu01"])) <= 0.003}
    choice = min(ok, key=lambda k: ok[k]["sd_we"]) if ok else None
    res["choice"] = choice
    if choice:
        c = res["models"][choice]
        # literature BO -> experiment correction (our own exact-BO numbers vs experiment, B3): nu01 BO(mu_atom) - exp
        corr_nu = EXP["nu01"] - pe["nu01"]
        corr_we = EXP["we"] - pe["we"]
        we_c = c["dmc"]["we"] + corr_we
        nu_c = c["dmc"]["nu01"] + corr_nu
        res["corrected"] = dict(we_BO=c["dmc"]["we"], we_BO_err=c["sd_we"], nu01_BO=c["dmc"]["nu01"], nu01_BO_err=c["sd_nu01"],
                                bo_to_exp_we=corr_we, bo_to_exp_nu01=corr_nu,
                                we=we_c, nu01=nu_c, dev_we=we_c / EXP["we"] - 1, dev_nu01=nu_c / EXP["nu01"] - 1,
                                dev_we_sigma=(we_c - EXP["we"]) / c["sd_we"], dev_nu01_sigma=(nu_c - EXP["nu01"]) / c["sd_nu01"],
                                dev_vs_exactBO_we=c["dmc"]["we"] / pe["we"] - 1, dev_vs_exactBO_nu01=c["dmc"]["nu01"] / pe["nu01"] - 1)
        res["resolved"] = bool(abs(res["corrected"]["dev_we_sigma"]) <= 2 and abs(res["corrected"]["dev_nu01_sigma"]) <= 2)
        print("choice", choice, res["corrected"], "resolved", res["resolved"])
    # local quartic cross-check (omega_e only), window 0.35 / 0.45 bohr
    re0 = Curva.fit(pts, e_inf)[0].minimo()[0]
    res["local_quartic"] = {f"w{w}": dict(dmc=local_poly(pts, re0, w, 4)["we"], exact_points=local_poly(ex, re0, w, 4)["we"])
                            for w in (0.35, 0.45)}
    print("local quartic", res["local_quartic"])
    res["time_s"] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, "h2vib_fix.json"), "w"), indent=1, default=float)
