"""fw1vuoto_analyze.py - exams of lab/out/fw1vuoto_prereg.json from the raw runs (lab/out/fw1vuoto_raw_w0/1.json).
Output lab/out/fw1vuoto_exams.json. Debiti's numbers (validation E1) are read from origin/master 4f50e2c with
`git show` (no checkout). Run: python lab/fw1vuoto_analyze.py"""
import os
import sys
import json
import subprocess
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fw1vuoto_core import pitR_uehling, LOGCOEF

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
D = {}
for w in (0, 1):
    fn = os.path.join(OUT, f"fw1vuoto_raw_w{w}.json")
    if os.path.exists(fn):
        D.update(json.load(open(fn)))
TH = [0.1, 0.05, 0.025, 0.0125]
L = LOGCOEF
T = 2


def bp(th):
    return 5.0 if th > 0.02 else 6.0


def bc(th):
    return 5.0 if th > 0.02 else 6.5


def raw(meth, rule, mode, th, d, qm, p1="V", p2="V", n=64, beta=None):
    if beta is None:
        beta = bp(th) if mode == "ball" else (6.0 if mode == "cellS" else bc(th))
    k = "%s|%s|%s|%g|%s|%.8g|%s|%s|%d|%g" % (meth, rule, mode, th, d, qm, p1, p2, n, beta)
    v = D.get(k)
    return None if v is None else complex(v[0], v[1])


def pit(meth, rule, mode, th, d, qm, **kw):
    v = raw(meth, rule, mode, th, d, qm, **kw)
    if v is None:
        return np.nan
    q = qm * th
    if meth == "E":
        return -v.real / T / q**2
    return v.real / q**2


def pit0(meth, rule, mode, th, d, **kw):
    return (4 * pit(meth, rule, mode, th, d, 0.025, **kw) - pit(meth, rule, mode, th, d, 0.05, **kw)) / 3


def piR(meth, rule, mode, th, qm, dirs="xyz", **kw):
    return float(np.mean([pit0(meth, rule, mode, th, d, **kw) - pit(meth, rule, mode, th, d, qm, **kw) for d in dirs]))


def logfit(th, y):
    th = np.array(th)
    A = np.stack([-np.log(th), np.ones_like(th), th**2], 1)
    co = np.linalg.lstsq(A, np.array(y), rcond=None)[0]
    return co


def aa(x):
    return np.inf if x is None or (isinstance(x, float) and np.isnan(x)) else abs(x)


def r(x, n=6):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else float(np.round(x, n))


EX = {}

# ---------------------------------------------------------------- E1 validation vs Debiti
try:
    deb = {}
    for f in ("vuoto_3d.json", "vuoto_3d_b.json"):
        deb.update(json.loads(subprocess.run(["git", "show", "4f50e2c:lab/out/" + f], capture_output=True, text=True,
                                             cwd=HERE).stdout))
    diffs = []
    for th in TH:
        for d in "xyz":
            for qm in (0.025, 0.05, 0.3, 1, 3, 10):
                k = "proj|%g|%s|%.8g|64|%g|0.6|1.05" % (th, d, qm, bp(th))
                mine = pit("S", "F0", "ball", th, d, qm)
                if k in deb and not np.isnan(mine):
                    diffs.append(abs(mine / deb[k][0] - 1))
    EX["E1"] = {"n_compared": len(diffs), "max_rel_diff": float(max(diffs)) if diffs else None,
                "PASS": bool(diffs and max(diffs) <= 1e-9)}
except Exception as e:                                   # pragma: no cover
    EX["E1"] = {"error": str(e)}

# ---------------------------------------------------------------- E2 frozen whole sea
E2 = {}
qs = [0.05, 0.1, 0.2, 0.4, 1, 3, 10]
ok = True
for th in (0.1, 0.05, 0.025):
    v = np.cos(th)
    rows = []
    for a, b in zip(qs[:-1], qs[1:]):
        pa, pb = pit("S", "F0", "cellS", th, "x", a), pit("S", "F0", "cellS", th, "x", b)
        meas = pb - pa
        cone = 4 * np.log(b / a) / (6 * np.pi**2 * v)
        pred = cone - 4 * (pitR_uehling(b) - pitR_uehling(a))
        rows.append({"q": [a, b], "meas": r(meas, 8), "pred": r(pred, 8), "err_of_cone_part": r((meas - pred) / cone, 5)})
        if th in (0.05, 0.025) and b <= 3 and not aa((meas - pred) / cone) <= 0.03:
            ok = False
    slope_ir = (pit("S", "F0", "cellS", th, "x", 0.1) - pit("S", "F0", "cellS", th, "x", 0.05)) / np.log(2)
    slope_uv = (pit("S", "F0", "cellS", th, "x", 10) - pit("S", "F0", "cellS", th, "x", 3)) / np.log(10 / 3)
    E2[str(th)] = {"rows": rows, "IR_slope_per_cell_over_L": r(slope_ir / L, 5),
                   "IR_pred_4/cos": r(4 / v, 5),
                   "UV_window_slope_3m_10m_over_L": r(slope_uv / L, 5),
                   "UV_pred_over_L": r((4 * np.log(10 / 3) / v - 4 * (pitR_uehling(10) - pitR_uehling(3)) / L) / np.log(10 / 3), 5),
                   "pit_cell": {str(q): r(pit("S", "F0", "cellS", th, "x", q), 8) for q in qs}}
conv = [pit("S", "F0", "cellS", 0.025, "x", q, n=96) / pit("S", "F0", "cellS", 0.025, "x", q) - 1 for q in (0.05, 3)]
E2["quadrature_n96_vs_n64_rel"] = [r(x, 9) for x in conv]
E2["PASS"] = ok
EX["E2"] = E2

# ---------------------------------------------------------------- E3 FW1 whole sea (stroboscopic)
E3 = {}
slopes = {str(th): r(np.mean([(pit("S", "FW1", "cellG", th, d, 0.05) - pit("S", "FW1", "cellG", th, d, 0.025)) / np.log(2)
                              for d in "xyz"]) / L, 6) for th in TH}
E3["E3a_IR_slope_over_L"] = slopes
E3["E3a_PASS"] = bool(all(aa(s) <= 0.02 for s in slopes.values()))
p0 = [np.mean([pit0("S", "FW1", "cellG", th, d) for d in "xyz"]) for th in TH]
co = logfit(TH, p0)
fd = (p0[3] - p0[2]) / np.log(2)
E3["pit0_avg"] = {str(th): r(p, 8) for th, p in zip(TH, p0)}
E3["E3b_L_fit_over_L"] = r(co[0] / L, 6)
E3["E3b_last_diff_over_L"] = r(fd / L, 6)
E3["E3b_const_c"] = r(co[1], 6)
E3["E3b_fd_steps_over_L"] = [r((p0[i + 1] - p0[i]) / np.log(2) / L, 6) for i in range(3)]
E3["E3b_PASS"] = bool(abs(co[0] / L - 1) <= 0.02 and abs(fd / L - 1) <= 0.02)
ue = {}
for th in TH:
    ue[str(th)] = {str(qm): r(piR("S", "FW1", "cellG", th, qm) / pitR_uehling(qm) - 1, 6) for qm in (0.2, 0.3, 1, 3, 10)}
E3["uehling_err"] = ue
e = ue["0.0125"]
E3["E3c_PASS"] = bool(all(aa(e[str(q)]) <= 0.005 for q in (0.3, 1, 3)) and aa(e["10"]) <= 0.015)
cv = {}
for qm in (1, 3):
    a, b = ue["0.025"][str(qm)] or np.nan, ue["0.0125"][str(qm)] or np.nan
    cv[str(qm)] = {"ratio": r(a / b if b else np.inf, 3), "ok": bool(aa(a / b) >= 3 or aa(b) < 1e-3)}
E3["E3d"] = cv
E3["E3d_PASS"] = bool(all(x["ok"] for x in cv.values()))
E3["E3e_q02_ratio"] = r(1 + (e["0.2"] or np.nan), 6)
E3["E3e_PASS"] = bool(aa(e["0.2"]) <= 0.01)
pd_ = pit0("S", "FW1", "cellG", 0.0125, "d") - pit("S", "FW1", "cellG", 0.0125, "d", 3)
E3["E3f_diag_vs_avg"] = r(pd_ / piR("S", "FW1", "cellG", 0.0125, 3) - 1, 6)
E3["E3f_PASS"] = bool(aa(E3["E3f_diag_vs_avg"]) <= 0.01)
an = {}
for th in TH:
    x, y, z = [pit0("S", "FW1", "cellG", th, d) for d in "xyz"]
    xf, yf, zf = [pit0("S", "F0", "ball", th, d) for d in "xyz"]
    an[str(th)] = {"FW1_(z-x)/L": r((z - x) / L, 6), "FW1_(y-(x+z)/2)/L": r((y - (x + z) / 2) / L, 6),
                   "F0proj_(z-x)/L": r((zf - xf) / L, 6), "F0proj_(y-(x+z)/2)/L": r((yf - (xf + zf) / 2) / L, 6)}
E3["E3g_anisotropy"] = an
E3["E3g_PASS"] = bool(aa(an["0.0125"]["FW1_(z-x)/L"]) > 1e-6)
cq = {}
for th in (0.025, 0.0125):
    cq[str(th)] = [r(pit("S", "FW1", "cellG", th, "x", q, n=80) / pit("S", "FW1", "cellG", th, "x", q) - 1, 9) for q in (0.05, 3)]
E3["quadrature_n80_vs_n64"] = cq
cq2 = [r(pit("S", "FW1mid", "cellG", 0.025, "x", q, beta=5.0) / pit("S", "FW1", "cellG", 0.025, "x", q) - 1, 7) for q in (0.05, 1, 3)]
E3["potential_mid_vs_end_rel_pit"] = cq2
pRmid = (pit0("S", "FW1mid", "cellG", 0.025, "x", beta=5.0) if False else None)
EX["E3"] = E3

# ---------------------------------------------------------------- E4 projector control
E4 = {}
for rule in ("F0", "FW1"):
    p0r = [np.mean([pit0("S", rule, "ball", th, d) for d in "xyz"]) for th in TH]
    cor = logfit(TH, p0r)
    E4[rule + "_proj_L_over_L"] = r(cor[0] / L, 6)
    E4[rule + "_proj_uehling_err_0.0125"] = {str(q): r(piR("S", rule, "ball", 0.0125, q) / pitR_uehling(q) - 1, 6)
                                            for q in (0.2, 0.3, 1, 3, 10)}
    E4[rule + "_proj_pit0"] = {str(th): r(p, 8) for th, p in zip(TH, p0r)}
d1 = {str(q): r(piR("S", "FW1", "ball", 0.0125, q) / piR("S", "F0", "ball", 0.0125, q) - 1, 6) for q in (0.3, 1, 3)}
d2 = {str(q): r(piR("S", "FW1", "cellG", 0.0125, q) / piR("S", "FW1", "ball", 0.0125, q) - 1, 6) for q in (0.3, 1, 3)}
E4["FW1proj_vs_F0proj"] = d1
E4["FW1whole_vs_FW1proj"] = d2
E4["PASS"] = bool(all(aa(v) <= 0.003 for v in d1.values()) and aa((E4["FW1_proj_L_over_L"] or np.nan) - 1) <= 0.02
                  and all(aa(v) <= 0.005 for v in d2.values()))
EX["E4"] = E4

# ---------------------------------------------------------------- E5 sea energy: chi_E, K_T, light speed
E5 = {}
AX = "xyz"


def Kraw(th, i, j, qm, sg=""):
    v = raw("E", "FW1", "cellG", th, sg + AX[i], qm, "A" + AX[j], "A" + AX[j])
    return None if v is None else v.real / T


wk = raw("E", "FW1", "cellG", 0.05, "x", 0.05, "Ax", "Ax", beta=5.0)
kt = Kraw(0.05, 0, 1, 0.05)
E5["E5a_ward_KL_over_KT"] = r(abs(wk.real / T) / abs(kt), 12) if (wk is not None and kt) else None
E5["E5a_PASS"] = bool(E5["E5a_ward_KL_over_KT"] is not None and E5["E5a_ward_KL_over_KT"] <= 1e-6)
mass = {}
pb = {}
for th in TH:
    rows = {}
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            K1, K2 = Kraw(th, i, j, 0.025), Kraw(th, i, j, 0.05)
            if K1 is None or K2 is None:
                continue
            q1, q2 = 0.025 * th, 0.05 * th
            k2 = (K2 - K1) / (q2**2 - q1**2)
            k0 = K1 - k2 * q1**2
            rows[AX[i] + AX[j]] = {"K0_over_K2m2": r(k0 / (k2 * th**2), 9),
                                   "pit_b0": (4 * K1 / q1**2 - K2 / q2**2) / 3}
    mass[str(th)] = {k: v["K0_over_K2m2"] for k, v in rows.items()}
    pb[th] = {k: v["pit_b0"] for k, v in rows.items()}
E5["E5b_photon_mass_K0_over_K2m2"] = mass
E5["E5b_PASS"] = bool(all(aa(v) <= 1e-4 for m in mass.values() for v in m.values() if v is not None))
pbavg = [np.mean(list(pb[th].values())) for th in TH]
cob = logfit(TH, pbavg)
E5["E5c_Lb_fit_over_L"] = r(cob[0] / L, 6)
E5["E5c_fd_steps_over_L"] = [r((pbavg[i + 1] - pbavg[i]) / np.log(2) / L, 6) for i in range(3)]
E5["E5c_PASS"] = bool(abs(cob[0] / L - 1) <= 0.02)
pe = {th: {d: pit0("E", "FW1", "cellG", th, d) for d in "xyz"} for th in TH}
pea = [np.mean(list(pe[th].values())) for th in TH]
coe = logfit(TH, pea)
E5["E5d_Le_fit_over_L"] = r(coe[0] / L, 6)
E5["E5d_uehling_err_0.0125"] = {str(q): r(piR("E", "FW1", "cellG", 0.0125, q) / pitR_uehling(q) - 1, 6) for q in (0.2, 0.3, 1, 3)}
E5["E5d_PASS"] = bool(abs(coe[0] / L - 1) <= 0.02 and all(aa(E5["E5d_uehling_err_0.0125"][str(q)]) <= 0.005 for q in (0.3, 1, 3)))
E5["pit_e0_slice_avg"] = {str(th): {d: r(v, 8) for d, v in pe[th].items()} for th in TH}
E5["pit_b0"] = {str(th): {k: r(v, 8) for k, v in pb[th].items()} for th in TH}
# nu_kk consistency: (q i, A j) and (q j, A i) both give B along k
cons = {}
for th in TH:
    for k in range(3):
        i, j = [a for a in range(3) if a != k]
        a, b = pb[th].get(AX[i] + AX[j]), pb[th].get(AX[j] + AX[i])
        if a is not None and b is not None:
            cons[f"{th}|B{AX[k]}"] = r((a - b) / L, 8)
E5["nu_consistency_(qiAj - qjAi)/L"] = cons
# light speed per (n, e): delta c/(alpha c) = 2 pi [pit_b(B axis) - pit_e(e)]
dc = {}
for th in TH:
    nub = {}
    for k in range(3):
        i, j = [a for a in range(3) if a != k]
        nub[AX[k]] = 0.5 * (pb[th].get(AX[i] + AX[j], np.nan) + pb[th].get(AX[j] + AX[i], np.nan))
    row = {}
    for n in range(3):
        for e_ in range(3):
            if e_ == n:
                continue
            b = 3 - n - e_
            row[f"n{AX[n]}_e{AX[e_]}"] = 2 * np.pi * (nub[AX[b]] - pe[th][AX[e_]])
    row["iso"] = float(np.mean(list(row.values())))
    dc[th] = row
E5["delta_c_over_alpha"] = {str(th): {k: r(v, 6) for k, v in dc[th].items()} for th in TH}
ext = {}
ths = np.array(TH[1:])
for k in dc[TH[0]]:
    y = np.array([dc[th][k] for th in TH[1:]])
    co1 = np.polyfit(ths, y, 1)
    ext[k] = r(co1[1], 6)
E5["delta_c_a_to_0_line_in_theta_3pts"] = ext
EX["E5"] = E5

# ---------------------------------------------------------------- E6 odd responses
E6 = {}
mixed = {}
for th in (0.05, 0.025, 0.0125):
    for i in range(3):
        for m in range(3):
            if m == i:
                continue
            row = {}
            for qm in (0.05, 0.1):
                vp = raw("E", "FW1", "cellG", th, AX[i], qm, "V", "A" + AX[m])
                vm = raw("E", "FW1", "cellG", th, "-" + AX[i], qm, "V", "A" + AX[m])
                if vp is None or vm is None:
                    continue
                q = qm * th
                odd = (vp - vm) / 2 / T
                even = (vp + vm) / 2 / T
                row[str(qm)] = {"odd_over_qL": [r(odd.real / q / L, 10), r(odd.imag / q / L, 10)],
                                "even_over_q2L": [r(even.real / q**2 / L, 8), r(even.imag / q**2 / L, 8)]}
            mixed[f"{th}|q{AX[i]}|A{AX[m]}"] = row
E6["mixed_V_A"] = mixed
anti = {}
for th in (0.05, 0.0125):
    for i in range(3):
        j, k = [a for a in range(3) if a != i]
        for qm in (0.05, 0.1):
            vp = raw("E", "FW1", "cellG", th, AX[i], qm, "A" + AX[j], "A" + AX[k])
            vm = raw("E", "FW1", "cellG", th, "-" + AX[i], qm, "A" + AX[j], "A" + AX[k])
            if vp is None or vm is None:
                continue
            q = qm * th
            odd = (vp - vm) / 2 / T
            even = (vp + vm) / 2 / T
            anti[f"{th}|q{AX[i]}|A{AX[j]}A{AX[k]}|{qm}"] = {"odd_over_qL": [r(odd.real / q / L, 10), r(odd.imag / q / L, 10)],
                                                           "even_over_q2L": [r(even.real / q**2 / L, 8), r(even.imag / q**2 / L, 8)]}
E6["AA_offdiag"] = anti
st = {}
for rule, mode in (("F0", "ball"), ("FW1", "cellG")):
    for i in range(3):
        for m in range(3):
            if m == i:
                continue
            for qm in (0.05, 0.1):
                v = raw("S", rule, mode, 0.05, AX[i], qm, "A" + AX[m], "A" + AX[m], beta=5.0)
                if v is None:
                    continue
                q = qm * 0.05
                st[f"{rule}|q{AX[i]}|A{AX[m]}|{qm}"] = [r(v.real / q**2 / L, 10), r(v.imag / q**2 / L, 10)]
E6["strobo_density_to_A_over_q2L"] = st
imchi = []
for th in TH:
    for d in "xyz":
        for qm in (0.3, 1, 3):
            v = raw("E", "FW1", "cellG", th, d, qm)
            if v is not None:
                imchi.append(aa(v.imag / v.real))
E6["E6d_max_Im_over_Re_chiE"] = float(max(imchi)) if imchi else None
EX["E6"] = E6

json.dump(EX, open(os.path.join(OUT, "fw1vuoto_exams.json"), "w"), indent=1)
print(json.dumps(EX, indent=1)[:20000])
