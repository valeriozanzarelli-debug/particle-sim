"""emiss6_judge.py - verdicts of pre-registration v6 from lab/out/emiss6_res_*.json -> lab/out/emiss6_results.json"""
import json, os, glob
import numpy as np
from zz_core import ALPHA, T_UNIT, MC2_EV
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
HC = 1239.841984
MU = 1 / (1 + 1 / 1836.15267343)
DIRAC = dict(A60=1.0886478591603033e-05, A1=8.074234457993113e-13, f60=0.2536534849842945, w1=1.996959026384637e-05,
             f1=0.13872879997366674 + 0.2774582983611244)
G = np.sqrt(1 - (60 * ALPHA)**2)


def ld(j):
    f = os.path.join(OUT, f"emiss6_res_{j}.json")
    return json.load(open(f)) if os.path.exists(f) else None


def V(ok):
    return "PASS" if ok else "FAIL"


R = {}
ex = {}
g0 = os.path.join(OUT, "emiss6_g0.json")
if os.path.exists(g0):
    d = json.load(open(g0)); ex["G0"] = dict(misurato=d["max_rel_diff"], soglia="< 1e-3", esito=V(d["max_rel_diff"] < 1e-3))
kp = os.path.join(OUT, "emiss_ck", "emiss6_kernelP_128_1.npz")
if os.path.exists(kp):
    K = np.load(kp)["K"]
    r = float(max(np.abs(K[p, -1, :, c]).max() / np.abs(K[p, 1, :, c]).max() for p in (0, 1) for c in range(3)))
    Kb = np.load(os.path.join(OUT, "emiss_ck", "emiss6_kernel_128_1.npz"))["K"]
    rb = float(max(np.abs(Kb[p, -1, :, c]).max() / np.abs(Kb[p, 1, :, c]).max() for p in (0, 1) for c in range(3)))
    ex["GK"] = dict(misurato=r, produzione_cieca=rb, soglia="< 2e-3", esito=V(r < 2e-3))
q = {a: ld(f"q60_{a}") for a in ("0.1", "0.05", "0.025")}
if q["0.025"]:
    d = q["0.025"]
    ex["X1"] = dict(misurato=d["max_dev_exp"], norma_fine=d["norm_end"], soglia="<= 0.005", esito=V(d["max_dev_exp"] <= 0.005))
    ex["X3"] = dict(misurato=d["A_over_W2"], dirac=DIRAC["A60"], scarto=d["A_over_W2"] / DIRAC["A60"] - 1,
                    s_inv=d["A_over_W2"] / T_UNIT, soglia="1.5 %", esito=V(abs(d["A_over_W2"] / DIRAC["A60"] - 1) <= 0.015))
    ex["X4"] = dict(misurato=d["max_p2p_minus_norm"], soglia="<= 2e-3", esito=V(d["max_p2p_minus_norm"] <= 2e-3))
    ex["X5"] = dict(misurato=d["omega_emit_over_rule"] - 1, soglia="2e-4", esito=V(abs(d["omega_emit_over_rule"] - 1) <= 2e-4))
    if "fwhm_over_A" in d:
        ex["X6"] = dict(fwhm_su_A=d["fwhm_over_A"], residuo=d["lor_resid"], centro_su_regola=d["line_centre_over_rule"],
                        soglia="0.05 e 0.05", esito=V(abs(d["fwhm_over_A"] - 1) <= 0.05 and d["lor_resid"] <= 0.05))
if all(q.values()):
    v = {a: q[a]["A_fit_over_rule"] for a in q}
    ex["X2"] = dict(misurato=v, finestre_0025=q["0.025"]["windows_over_rule"], soglia="1 %", esito=V(all(abs(x - 1) <= 0.01 for x in v.values())))
    aa = np.array([0.1, 0.05, 0.025]); yy = np.array([q[a]["A_over_W2"] for a in ("0.1", "0.05", "0.025")])
    M = np.stack([np.ones(3), aa**(2 * G), aa**2], 1)
    c = np.linalg.solve(M, yy)
    ex["X3x"] = dict(estrapolato=c[0], scarto=c[0] / DIRAC["A60"] - 1, esito="info")
rb = {a: ld(f"rabi60_{a}") for a in ("0.1", "0.05", "0.025")}
if all(rb.values()):
    v = {a: rb[a]["Omega_over_rule"] for a in rb}
    ex["B1"] = dict(misurato=v, soglia="2e-3", esito=V(all(abs(x - 1) <= 2e-3 for x in v.values())))
    ex["B1f"] = dict(f32_rabi={a: rb[a]["f32"] for a in rb}, esito="info")
sc = ld("scan60")
if sc:
    ex["B4"] = dict(misurato=sc["scan"], max_diff=sc["max_abs_diff"], soglia="0.02", esito=V(sc["max_abs_diff"] <= 0.02))
pk = {(Z, a): ld(f"pk{Z}_{a}") for Z, aa_ in ((60, ("0.1", "0.05", "0.025")), (1, ("0.5", "0.25", "0.125"))) for a in aa_}
have = {k: v for k, v in pk.items() if v}
if have:
    ex["B7"] = dict(misurato={f"Z{k[0]}_a{k[1]}": v["R_over_rule"] for k, v in have.items()}, soglia="1e-3",
                    esito=V(all(abs(v["R_over_rule"] - 1) <= 1e-3 for v in have.values())) + ("" if len(have) == 6 else " (parziale)"))
if pk[(60, "0.025")]:
    f = pk[(60, "0.025")]["f32"]
    ex["B3"] = dict(misurato=f, dirac=DIRAC["f60"], scarto=f / DIRAC["f60"] - 1, soglia="1.5 %", esito=V(abs(f / DIRAC["f60"] - 1) <= 0.015))
h = {a: ld(f"q1_{a}") for a in ("0.5", "0.25", "0.125")}
p1 = {a: pk[(1, a)] for a in ("0.5", "0.25", "0.125")}
aa = np.array([0.5, 0.25, 0.125]); Mh = np.stack([np.ones(3), aa**2, aa**4], 1)
if all(h.values()):
    v = {a: h[a]["A_fit_over_rule"] for a in h}
    ex["H1"] = dict(misurato=v, finestre={a: h[a]["windows_over_rule"] for a in h}, soglia="1 %", esito=V(all(abs(x - 1) <= 0.01 for x in v.values())))
    A = np.linalg.solve(Mh, np.array([h[a]["A_over_W2"] for a in ("0.5", "0.25", "0.125")]))[0]
    ex["H2"] = dict(per_a={a: h[a]["A_over_W2"] / T_UNIT for a in h}, estrapolato=A, s_inv=A / T_UNIT, vita_ns=1e9 * T_UNIT / A,
                    dirac_s_inv=DIRAC["A1"] / T_UNIT, scarto=A / DIRAC["A1"] - 1, soglia="0.5 %", esito=V(abs(A / DIRAC["A1"] - 1) <= 0.005))
    w = np.linalg.solve(Mh, np.array([h[a]["omega_emit"] for a in ("0.5", "0.25", "0.125")]))[0]
    lam = HC / (w * MC2_EV)
    em = max(abs(h[a]["omega_emit_over_rule"] - 1) for a in h)
    ex["H3"] = dict(omega_emit_su_regola_max=em, lambda_nm=lam, dirac_nm=HC / (DIRAC["w1"] * MC2_EV), scarto=lam / (HC / (DIRAC["w1"] * MC2_EV)) - 1,
                    soglia="1e-4 e 0.02 %", esito=V(em <= 1e-4 and abs(lam / (HC / (DIRAC["w1"] * MC2_EV)) - 1) <= 2e-4))
    ex["H4"] = dict(A_mu_s_inv=A / T_UNIT * MU, nist=6.2649e8, scarto_A=A / T_UNIT * MU / 6.2649e8 - 1, lambda_mu_nm=lam / MU,
                    misurata_nm=121.5668, scarto_lambda=lam / MU / 121.5668 - 1, esito="info")
if all(p1.values()):
    ff = np.linalg.solve(Mh, np.array([p1[a]["f_tot_nr"] for a in ("0.5", "0.25", "0.125")]))[0]
    ex["B2"] = dict(per_a={a: p1[a]["f_tot_nr"] for a in p1}, estrapolato=ff, atteso=0.4162, dirac=DIRAC["f1"], scarto=ff / 0.4162 - 1,
                    soglia="0.5 %", esito=V(abs(ff / 0.4162 - 1) <= 0.005))
    if "H2" in ex:
        lam_c = 2 * np.pi / w
        sig = 4 * np.pi * ALPHA * ff / A
        ex["B6"] = dict(sigma_picco_lC2=sig, atteso=3 * lam_c**2 / (2 * np.pi), rapporto=sig / (3 * lam_c**2 / (2 * np.pi)),
                        sigma_picco_cm2=sig * (3.8615926796e-11)**2, soglia="1 %", esito=V(abs(sig / (3 * lam_c**2 / (2 * np.pi)) - 1) <= 0.01))
vv = os.path.join(OUT, "emiss6_rabi_vivo_res.json")
if os.path.exists(vv):
    ex.update(json.load(open(vv)))
json.dump(ex, open(os.path.join(OUT, "emiss6_results.json"), "w"), indent=1, default=float, ensure_ascii=False)
for k, v in ex.items():
    print(k, v.get("esito"), {kk: vv_ for kk, vv_ in v.items() if kk in ("misurato", "scarto", "estrapolato", "max_diff", "rapporto")})
