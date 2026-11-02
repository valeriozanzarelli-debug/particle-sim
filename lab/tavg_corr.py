"""tavg_corr.py - FASE 2 dello studio tavg_: correzione di correlazione IMPARATA dagli esatti della regola.

Correlazione = E(esatta della regola) - E(Hartree-Fock della regola), per i pochi elettroni che sappiamo fare esatti
(ECG: H-, He, He 2^3S, Li+, Li, Li 2P, Be+, ioni He-like Z = 6..10). Si impara una legge E_c[rho_su, rho_giu] da queste
coppie (densita' HF -> correlazione esatta), poi la si applica ALLA CIECA all'ossigeno (e all'azoto).
Famiglie di funzioni (debito dichiarato: la FORMA viene dalla letteratura, i NUMERI dai nostri esatti):
  W: primo termine di Colle-Salvetti/LYP (tipo Wigner, solo coppie a spin opposto): 2 manopole (a, d)
  L: Colle-Salvetti nella forma LYP di Miehlich et al. 1989 (gradienti della densita'): 4 manopole (a, b, c, d)
Scelta pre-registrata: validazione incrociata lascia-uno-fuori (tavg_prereg2.json).
Uso: python tavg_corr.py train | predict | devtest"""
import os, sys, json, math, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_radial as tr
import tavg_blocks as tb

OUT = os.path.join(HERE, "out")
FN = os.path.join(OUT, "tavg_corr.json")
EH_EV = 27.211386245988
CF = 0.3 * (3 * math.pi ** 2) ** (2 / 3)
MA = 1.0003832623372946          # m_kin of the multi_e ECG runs (theta_A); E(theta->0) = E * m0 / m_A
M0 = 1.0000000045603115


def _ecg(name):
    d = json.load(open(os.path.join(OUT, "multi_e_ecg.json")))["prod"]["nr"][name]
    return d["E"] * M0 / d["m_kin"], max(d["conv_err"], 1e-5)


def _tav(name):
    d = json.load(open(os.path.join(OUT, "tav_runs.json")))[name]
    return d["E_theta0"], max(d["barra"], 1e-5)


def training_set():
    """(name, Z, config, E_exact(theta->0), sigma, source) - only exact results of the rule's own machines"""
    T = []
    E, s = _ecg("he1S_Z1"); T.append(("H-", 1, "1s2", E, s, "multi_e_ecg.json he1S_Z1"))
    he = json.load(open(os.path.join(OUT, "multi_e_v5.json")))["he"]["M1v5"]
    T.append(("He", 2, "1s2", he["E_theta0"], max(he["error_bar"], 1e-5), "multi_e_v5.json M1v5"))
    E, s = _ecg("he3S_Z2"); T.append(("He 2^3S", 2, "1s1 2s1", E, s, "multi_e_ecg.json he3S_Z2"))
    E, s = _tav("liP1S"); T.append(("Li+", 3, "1s2", E, s, "tav_runs.json liP1S"))
    E, s = _tav("li2S"); T.append(("Li", 3, "1s2 2s1", E, s, "tav_runs.json li2S"))
    E, s = _tav("li2P"); T.append(("Li 2P", 3, "1s2 2p1", E, s, "tav_runs.json li2P"))
    E, s = _tav("beP2S"); T.append(("Be+", 4, "1s2 2s1", E, s, "tav_runs.json beP2S"))
    for z in (6, 7, 8, 9, 10):
        E, s = _ecg(f"he1S_Z{z}"); T.append((f"Z{z} He-like", z, "1s2", E, s, f"multi_e_ecg.json he1S_Z{z}"))
    return T


_SD = {}


def spin_densities(a):
    if id(a) in _SD:
        return _SD[id(a)]
    _SD[id(a)] = _spin_densities(a)
    return _SD[id(a)]


def _spin_densities(a):
    """rho_up, rho_dn (spherical) and their r-derivatives from the Hund determinant of a radial ROHF atom"""
    r = a.g.r
    det = tr.hund_determinant(a.config)
    nu = {}; nd = {}
    for (s, m, sp) in det:
        (nu if sp > 0 else nd)[s] = (nu if sp > 0 else nd).get(s, 0) + 1
    ru = sum(nu.get(i, 0) * a.P[i] ** 2 for i in range(len(a.config))) / (4 * math.pi * r * r)
    rd = sum(nd.get(i, 0) * a.P[i] ** 2 for i in range(len(a.config))) / (4 * math.pi * r * r)
    if np.isscalar(rd):
        rd = np.zeros_like(r)
    x = a.g.x
    gu = np.gradient(ru, x) / r; gd = np.gradient(rd, x) / r
    return r, a.g.h, ru, rd, gu, gd


def features(a, c, d):
    """X1 = int 4 rho_a rho_b / (rho (1 + d rho^-1/3));  X2 = int omega {...} (Miehlich LYP), so that
    E_c = -a X1 - a b X2"""
    r, hx, ra, rb, ga, gb = spin_densities(a)
    rho = ra + rb
    ok = rho > 1e-14
    rho = np.where(ok, rho, 1e-14); ra = np.where(ok, ra, 0); rb = np.where(ok, rb, 0)
    g = ga + gb
    r13 = rho ** (-1 / 3)
    den = 1 + d * r13
    w = 4 * math.pi * r ** 3 * hx                    # int f d^3r = sum f 4 pi r^2 (r dx)
    X1 = float(np.sum(4 / den * ra * rb / rho * w))
    om = np.exp(-c * r13) / den * rho ** (-11 / 3)
    de = c * r13 + d * r13 / den
    ga2, gb2, g2 = ga * ga, gb * gb, g * g
    br = (ra * rb * (2 ** (11 / 3) * CF * (ra ** (8 / 3) + rb ** (8 / 3)) + (47 / 18 - 7 * de / 18) * g2
                     - (2.5 - de / 18) * (ga2 + gb2) - (de - 11) / 9 * (ra / rho * ga2 + rb / rho * gb2))
          - 2 / 3 * rho ** 2 * g2 + (2 / 3 * rho ** 2 - ra ** 2) * gb2 + (2 / 3 * rho ** 2 - rb ** 2) * ga2)
    X2 = float(np.sum(np.where(ok, om * br, 0.0) * w))
    return X1, X2


def hf_atom(Z, cfg):
    return tb.get_ae(Z, cfg)


def devtest():
    """code check only (declared): LYP with the published constants on our HF He, against the known ~ -0.0437 Eh"""
    a = hf_atom(2, "1s2")
    X1, X2 = features(a, 0.2533, 0.349)
    print("He HF", a.E, "LYP(published constants) =", -0.04918 * X1 - 0.04918 * 0.132 * X2)
    a = hf_atom(10, "1s2 2s2 2p6")
    X1, X2 = features(a, 0.2533, 0.349)
    print("Ne HF", a.E, "LYP(published constants) =", -0.04918 * X1 - 0.04918 * 0.132 * X2)


def fit(rows, model, cgrid, dgrid):
    """rows: list of dict(atom, Ec, sigma). Weighted least squares: linear in (a, ab) [L] or a [W] for each (c, d)"""
    best = None
    for c in (cgrid if model == "L" else [0.0]):
        for d in dgrid:
            F = np.array([features(r["atom"], c, d) for r in rows])
            y = np.array([r["Ec"] for r in rows]); wts = 1 / np.array([r["sigma_c"] for r in rows])
            if model == "L":
                A = -F
            else:
                A = -F[:, :1]
            coef, *_ = np.linalg.lstsq(A * wts[:, None], y * wts, rcond=None)
            res = A @ coef - y
            chi = float(np.sum((res * wts) ** 2))
            if best is None or chi < best["chi2"]:
                best = dict(c=c, d=d, coef=coef.tolist(), chi2=chi, rms_mEh=float(1e3 * np.sqrt(np.mean(res ** 2))))
    return best


def predict(model, p, atom):
    X1, X2 = features(atom, p["c"], p["d"])
    if model == "L":
        return -p["coef"][0] * X1 - p["coef"][1] * X2
    return -p["coef"][0] * X1


def grids(prereg):
    g = prereg["griglie"]
    return np.linspace(*g["c"]), np.linspace(*g["d"])


def train():
    pre = json.load(open(os.path.join(OUT, "tavg_prereg2.json")))
    cg, dg = grids(pre)
    rows = []
    for (name, Z, cfg, E, s, src) in training_set():
        a = hf_atom(Z, cfg)
        # sigma on E_c: exact-run bar plus 1e-6 for HF (radial, ~1e-8 measured)
        rows.append(dict(name=name, Z=Z, config=cfg, E_exact=E, E_HF=a.E, Ec=E - a.E, sigma_c=s + 1e-6, src=src,
                         atom=a))
        print(f"{name:12s} E_exact {E:.7f}  E_HF {a.E:.7f}  Ec {1e3 * (E - a.E):.3f} mEh (sigma {1e3 * s:.3f})",
              flush=True)
    out = dict(training=[{k: v for k, v in r.items() if k != "atom"} for r in rows], models={})
    for model in ("W", "L"):
        t0 = time.time()
        full = fit(rows, model, cg, dg)
        loo = []
        for i in range(len(rows)):
            sub = rows[:i] + rows[i + 1:]
            p = fit(sub, model, cg, dg)
            pred = predict(model, p, rows[i]["atom"])
            loo.append(dict(name=rows[i]["name"], Ec=rows[i]["Ec"], pred=pred, err_mEh=1e3 * (pred - rows[i]["Ec"])))
        cv = float(np.sqrt(np.mean([x["err_mEh"] ** 2 for x in loo])))
        out["models"][model] = dict(fit=full, loo=loo, cv_rms_mEh=cv, manopole=2 if model == "W" else 4,
                                    secondi=time.time() - t0)
        print(model, "fit", full, "CV rms", round(cv, 3), "mEh", flush=True)
    W, L = out["models"]["W"]["cv_rms_mEh"], out["models"]["L"]["cv_rms_mEh"]
    out["scelto"] = "L" if L < 0.9 * W else "W"
    json.dump(out, open(FN, "w"), indent=1)
    print("scelto:", out["scelto"], flush=True)


def blind():
    """blind predictions with the chosen model (no change of anything after seeing these)"""
    res = json.load(open(FN))
    model = res["scelto"]; p = res["models"][model]["fit"]
    targets = [("O", 8, "1s2 2s2 2p4"), ("O+", 8, "1s2 2s2 2p3"), ("N", 7, "1s2 2s2 2p3"), ("N+", 7, "1s2 2s2 2p2")]
    pr = {}
    for name, Z, cfg in targets:
        a = hf_atom(Z, cfg)
        ec = predict(model, p, a)
        pr[name] = dict(Z=Z, config=cfg, E_HF_fine=a.E, Ec_pred=ec, E_pred=a.E + ec)
        print(name, pr[name], flush=True)
        # also the other model, for the record
        om = "L" if model == "W" else "W"
        pr[name]["Ec_altro_modello"] = predict(om, res["models"][om]["fit"], a)
    res["cieco"] = pr
    json.dump(res, open(FN, "w"), indent=1)


if __name__ == "__main__":
    dict(train=train, predict=blind, devtest=devtest)[sys.argv[1]]()
