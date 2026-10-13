"""vuoto_f2_analyze.py - grades the phase-2 exams against lab/out/vuoto_f2_prereg.json (+ amendments, if any).
Reads lab/out/vuoto_f2_se.json, vuoto_f2_bethe.json, vuoto_f2_breath.json, vuoto_f2_hproj.json.
Writes lab/out/vuoto_f2_exams.json."""
import os, json
import numpy as np
from zz_core import ALPHA
from vuoto_core import pitR_uehling, LOGCOEF, PI

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def load(name):
    p = os.path.join(OUT, name)
    return json.load(open(p)) if os.path.exists(p) else {}


R = {}
A2P = ALPHA / (2 * PI)
# ------------------------------------------------------------------ S: dm log coefficients
SE = load("vuoto_f2_se.json")
for lat, tag in ((False, "S1"), (True, "S2")):
    ths, C, T = [], [], []
    for th in (0.05, 0.025, 0.0125, 0.00625, 0.003125):
        k = "dm|%s|%g|0|std" % ("lat" if lat else "cont", th)
        if k in SE:
            ths.append(th); C.append(SE[k]["dmC_over_m"]); T.append(SE[k]["dmT_over_m"])
    if len(ths) == 5:
        ths = np.array(ths)
        A = np.stack([np.ones(5), np.log(1 / ths), ths], 1)
        bC = np.linalg.lstsq(A, np.array(C), rcond=None)[0][1]
        bT = np.linalg.lstsq(A, np.array(T), rcond=None)[0][1]
        rC, rT, rt = bC / (ALPHA / PI), bT / (ALPHA / (2 * PI)), (bC + bT) / (3 * ALPHA / (2 * PI))
        d = {"theta": ths.tolist(), "dmC_over_m": C, "dmT_over_m": T, "B_C/(alpha/pi)": rC, "B_T/(alpha/2pi)": rT,
             "B_tot/(3alpha/2pi)": rt}
        if tag == "S1":
            d["PASS"] = bool(abs(rC - 1) <= 0.005 and abs(rT - 1) <= 0.005 and abs(rt - 1) <= 0.005)
        else:
            d["S2a_PASS"] = bool(abs(rC - 1) <= 0.02); d["S2b_PASS"] = bool(abs(rT - 1) <= 0.03)
            d["S2c_PASS"] = bool(abs(rt - 1) <= 0.02)
        R[tag] = d


# ------------------------------------------------------------------ G: g - 2
def g2_table(lat):
    out = {}
    for th in (0.05, 0.025, 0.0125):
        vals = {}
        for lf in (1e-3, 3e-4):
            ka, kb = ["g2|%s|%g|%g|std|%s" % ("lat" if lat else "cont", th, lf, g) for g in "AB"]
            if ka in SE and kb in SE:
                pr = dict(SE[ka]["pairs"]); pr.update(SE[kb]["pairs"])
                vals[lf] = dict(SE[ka]); vals[lf]["pairs"] = pr
                vals[lf]["F2_mean_over_alpha2pi"] = float(np.mean([v["F2_over_alpha2pi"] for v in pr.values()]))
        if len(vals) == 2:
            f1, f2 = vals[1e-3]["F2_mean_over_alpha2pi"], vals[3e-4]["F2_mean_over_alpha2pi"]
            f0 = f2 + (f2 - f1) * 3e-4 / (1e-3 - 3e-4)                 # linear in lam -> 0
            # per pair, lam -> 0
            pairs = {}
            for pk in vals[1e-3]["pairs"]:
                a, b = vals[1e-3]["pairs"][pk]["F2_over_alpha2pi"], vals[3e-4]["pairs"][pk]["F2_over_alpha2pi"]
                pairs[pk] = b + (b - a) * 3e-4 / (1e-3 - 3e-4)
            out[th] = {"F2_lam1e-3": f1, "F2_lam3e-4": f2, "F2_lam0": f0, "pairs_lam0": pairs,
                       "spread_lam0": max(pairs.values()) - min(pairs.values()),
                       "dv_over_alpha": vals[3e-4]["dv_over_alpha"], "dm_c_over_m": vals[3e-4]["dm_c_over_m"]}
    return out


for lat, tag in ((False, "G0"), (True, "G1")):
    tb = g2_table(lat)
    if len(tb) == 3:
        ths = np.array([0.05, 0.025, 0.0125])
        f = np.array([tb[t]["F2_lam0"] for t in ths])
        c = np.linalg.solve(np.stack([np.ones(3), ths, ths**2], 1), f)
        d = {"table": {str(k): v for k, v in tb.items()}, "F2_extrap_over_alpha2pi": float(c[0])}
        if tag == "G0":
            d["PASS"] = bool(abs(c[0] - 1) <= 0.005)
        else:
            d["G1_PASS"] = bool(abs(c[0] - 1) <= 0.03)
            d["G2_value"] = tb[0.0125]["F2_lam0"]; d["G2_PASS"] = bool(abs(tb[0.0125]["F2_lam0"] - 1) <= 0.06)
            d["G3_spreads"] = [tb[t]["spread_lam0"] for t in ths]
            d["G3_PASS"] = bool(tb[0.0125]["spread_lam0"] <= 0.5 * tb[0.05]["spread_lam0"])
        R[tag] = d
def mean6(prefix):
    ks = [prefix + "|" + g for g in "AB"]
    if not all(k in SE for k in ks):
        return None
    return float(np.mean([v["F2_over_alpha2pi"] for k in ks for v in SE[k]["pairs"].values()]))


if mean6("g2|lat|0.025|0.001|fine") is not None and mean6("g2|lat|0.025|0.001|std") is not None:
    R["G_grid_check"] = {"std": mean6("g2|lat|0.025|0.001|std"), "fine": mean6("g2|lat|0.025|0.001|fine")}

# ------------------------------------------------------------------ L: Bethe log (prereg: vuoto_f2_prereg_L.json)
BL = load("vuoto_f2_bethe.json")
NR = load("vuoto_f2_nrref.json")
REF = {1: 2.984128556, 2: 2.811769893}
KRY = [30, 100, 300, 1000, 3000]
if NR:
    def nrval(n, kr, field="lnk0_KS"):
        a_, b_ = NR["%d|0.02" % n]["K"][str(kr)][field], NR["%d|0.01" % n]["K"][str(kr)][field]
        return b_ + (b_ - a_) / 3                                   # Richardson h -> 0 (O(h^2))
    Ld = {}
    for th in (0.1, 0.05):
        for n in (1, 2):
            ks = {kp: BL.get("%d|%d|%g" % (n, kp, th)) for kp in (1, -2)}
            if not all(ks.values()):
                continue
            row = {}
            for i, kr in enumerate(KRY):
                L = ks[1]["lnk0_K"][i] / 3 + 2 * ks[-2]["lnk0_K"][i] / 3
                DK = (ks[1]["D_K"][i] / 3 + 2 * ks[-2]["D_K"][i] / 3) / ks[1]["D_sumrule_NR"]
                row[str(kr)] = {"L_walk": L, "L_NR": nrval(n, kr), "diff": L - nrval(n, kr),
                                "DK_walk": DK, "DK_NR": nrval(n, kr, "D_K_over_sumrule"),
                                "hybrid_lnk0": L + (REF[n] - nrval(n, kr))}     # rule below K, NR theory above K
            Ld["%d|%g" % (n, th)] = row
    R["L_table"] = Ld
    if "1|0.05" in Ld and "2|0.05" in Ld:
        R["L1_PASS"] = bool(all(abs(Ld[k][kr]["diff"]) <= 0.02 for k in ("1|0.05", "2|0.05") for kr in ("100", "300")))
        R["L3_PASS"] = bool(all(abs(Ld[k]["300"]["DK_walk"] / Ld[k]["300"]["DK_NR"] - 1) <= 0.01 for k in ("1|0.05", "2|0.05")))
        R["L_hybrid_lnk0_300Ry"] = {k: Ld[k]["300"]["hybrid_lnk0"] for k in ("1|0.05", "2|0.05")}
    if all(k in Ld for k in ("1|0.05", "1|0.1", "2|0.05", "2|0.1")):
        R["L2_PASS"] = bool(all(abs(Ld["%d|0.05" % n]["300"]["L_walk"] - Ld["%d|0.1" % n]["300"]["L_walk"]) <= 0.01 for n in (1, 2)))
    MC2_MHZ = 510998.95 * 2.417989242e14 / 1e6
    R["L_MHz_per_unit_lnk0"] = {n: 4 * ALPHA**5 / (3 * PI * n**3) * MC2_MHZ for n in (1, 2)}

# ------------------------------------------------------------------ B: breathing coin
BR = load("vuoto_f2_breath.json")


def bpit(kind, th, mu, dn, qm):
    return BR.get("%s|%g|%.6g|%s|%g" % (kind, th, mu, dn, qm))


Bd = {}
for th in (0.025, 0.0125):
    mu = float(np.sqrt(th))
    ok = True
    res = {"mu": mu}
    # physical species with / without breathing
    for q in (1, 3):
        rat = []
        for dn in "xyz":
            vals = [bpit("phys", th, m_, dn, x) for m_ in (0.0, mu) for x in (0.05, 0.025, q)]
            if any(v is None for v in vals):
                ok = False; break
            p00, p01, p0q, p10, p11, p1q = vals
            P0a = (4 * p01 - p00) / 3; P0b = (4 * p11 - p10) / 3
            rat.append((P0b - p1q) / (P0a - p0q))
        if ok:
            res["B7_ratio_q%d" % q] = float(np.mean(rat)); res["B7_ratio_q%d_xyz" % q] = rat
    # whole sea
    full = {}
    for dn in "xyz":
        vals = [bpit("full", th, mu, dn, x) for x in (0.05, 0.025, 1, 3)]
        if any(v is None for v in vals):
            ok = False; break
        full[dn] = vals
    if ok:
        P0 = np.mean([(4 * full[d][1] - full[d][0]) / 3 for d in "xyz"])
        pred0 = LOGCOEF * np.log(mu / th)
        res["B5_Pi0"] = float(P0); res["B5_pred"] = float(pred0); res["B5_ratio"] = float(P0 / pred0)
        for i, q in ((2, 1), (3, 3)):
            PRq = np.mean([(4 * full[d][1] - full[d][0]) / 3 - full[d][i] for d in "xyz"])
            pv = pitR_uehling(q) - pitR_uehling(q, m=mu / th)
            res["B6_q%d_ratio" % q] = float(PRq / pv)
        Bd[str(th)] = res
if Bd:
    R["B_table"] = Bd
    if len(Bd) == 2 and all("B5_ratio" in v for v in Bd.values()):
        R["B5_PASS"] = bool(all(abs(v["B5_ratio"] - 1) <= 0.03 for v in Bd.values()))
        R["B6_PASS"] = bool(all(abs(v["B6_q1_ratio"] - 1) <= 0.02 and abs(v["B6_q3_ratio"] - 1) <= 0.02 for v in Bd.values()))
        R["B7_PASS"] = bool(all(abs(v["B7_ratio_q1"] - 1) <= 0.01 and abs(v["B7_ratio_q3"] - 1) <= 0.01 for v in Bd.values()))
if "bands" in BR:
    errs = []
    for th, b in BR["bands"].items():
        mu = np.sqrt(float(th))
        for name, v in b.items():
            errs.append(abs(v - (float(th) if name.startswith("valle") else mu)))
    R["B8_maxerr"] = float(max(errs)); R["B8_PASS"] = bool(max(errs) <= 1e-9)

# ------------------------------------------------------------------ H: hydrogen with the projector on the electron
HP = load("vuoto_f2_hproj.json")
if all("%s|%g" % (n, t) in HP for n in ("1s", "2s") for t in (0.05, 0.025, 0.0125)):
    ref = HP["ref_MHz"]
    ths = np.array([0.05, 0.025, 0.0125])
    H = {}
    for n in ("1s", "2s"):
        v = np.array([HP["%s|%g" % (n, t)]["MHz"] for t in ths])
        c = np.linalg.solve(np.stack([np.ones(3), ths, ths**2], 1), v)
        H[n] = {"theta_atom": dict(zip(map(str, ths), v.tolist())), "extrap": float(c[0]), "ref": ref[n],
                "finest_ratio": float(v[-1] / ref[n]), "extrap_ratio": float(c[0] / ref[n])}
    R["H_table"] = H
    R["H1_PASS"] = bool(abs(H["1s"]["finest_ratio"] - 1) <= 0.01)
    R["H2_PASS"] = bool(abs(H["2s"]["finest_ratio"] - 1) <= 0.01)
    R["H3_PASS"] = bool(abs(H["1s"]["extrap_ratio"] - 1) <= 0.005 and abs(H["2s"]["extrap_ratio"] - 1) <= 0.005)
    f12 = H["2s"]["extrap"] - H["1s"]["extrap"]
    f12ref = ref["2s"] - ref["1s"]
    R["H4"] = {"value_MHz": f12, "ref_MHz": f12ref, "ratio": f12 / f12ref}
    R["H4_PASS"] = bool(abs(f12 / f12ref - 1) <= 0.01)

json.dump(R, open(os.path.join(OUT, "vuoto_f2_exams.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in R.items() if k.endswith("PASS") or k in ("S1", "S2")}, indent=1, default=str))
