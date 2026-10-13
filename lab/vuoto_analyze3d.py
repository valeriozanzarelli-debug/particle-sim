"""vuoto_analyze3d.py - grades the pre-registered 3-D exams V2a-V2g from the checkpoint tables. -> lab/out/vuoto_3d_exams.json"""
import os, json
import numpy as np
from vuoto_core import pitR_uehling, LOGCOEF
from vuoto_levels import load_3d

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
D = load_3d()
TH = [0.1, 0.05, 0.025, 0.0125]


def g(th, dn, qm, K=(0.6, 1.05), n=64):
    beta = 6.0 if th < 0.02 else 5.0
    if K != (0.6, 1.05):
        beta = 5.0
    return D["proj|%g|%s|%.8g|%d|%g|%g|%g" % (th, dn, qm, n, beta, K[0], K[1])][0]


def p0(th, dn, K=(0.6, 1.05)):
    return (4 * g(th, dn, 0.025, K) - g(th, dn, 0.05, K)) / 3


def piR(th, dn, qm, K=(0.6, 1.05)):
    return p0(th, dn, K) - g(th, dn, qm, K)


R = {}
# V2a whole sea
fq = sorted((float(k.split("|")[3]), v[0]) for k, v in D.items() if k.startswith("full|"))
fd = dict(fq)
slope = (fd[0.1] - fd[0.05]) / np.log(2)
pred = LOGCOEF / np.cos(0.1)
R["V2a"] = {"pit_whole": fd, "IR_slope": slope, "pred": pred, "rel": slope / pred - 1,
            "PASS": bool(abs(slope / pred - 1) <= 0.02 and all(v < 0 for v in fd.values()))}
# V2b log coefficient
P0 = {th: np.mean([p0(th, d) for d in "xyz"]) for th in TH}
A = np.array([[-np.log(t), 1, t * t] for t in TH])
L = np.linalg.lstsq(A, np.array([P0[t] for t in TH]), rcond=None)[0][0]
fd2 = (P0[0.0125] - P0[0.025]) / np.log(2)
R["V2b"] = {"pit0_xyz": {str(t): P0[t] for t in TH}, "fit_L": L, "fit_rel": L / LOGCOEF - 1,
            "fd_finest": fd2, "fd_rel": fd2 / LOGCOEF - 1,
            "fd_all": {f"{TH[i]}->{TH[i+1]}": (P0[TH[i + 1]] - P0[TH[i]]) / np.log(2) / LOGCOEF - 1 for i in range(3)},
            "PASS": bool(abs(L / LOGCOEF - 1) <= 0.02 and abs(fd2 / LOGCOEF - 1) <= 0.02)}
# V2c shape
tab = {}
for th in TH:
    tab[str(th)] = {str(q): np.mean([piR(th, d, q) for d in "xyz"]) / pitR_uehling(q) - 1 for q in (0.3, 1, 3, 10)}
    tab[str(th)]["per_dir_q1"] = {d: piR(th, d, 1) / pitR_uehling(1) - 1 for d in "xyzd"}
e = tab["0.0125"]
R["V2c"] = {"rel_err": tab, "PASS": bool(all(abs(e[str(q)]) <= 0.005 for q in (0.3, 1, 3)) and abs(e["10"]) <= 0.015)}
# V2d q^2 coefficient at q = 0.2 m (from the table: interpolate in log)
from vuoto_levels import load_pitR, make_interp
qs, tb, _ = load_pitR(0.0125)
f = make_interp(qs, tb)
r2 = float(f(0.2) / pitR_uehling(0.2))
R["V2d"] = {"ratio_q0.2": r2, "PASS": bool(abs(r2 - 1) <= 0.01)}
# V2e isotropy (1,1,1)
iso = piR(0.0125, "d", 3) / np.mean([piR(0.0125, d, 3) for d in "xyz"]) - 1
R["V2e"] = {"rel": iso, "PASS": bool(abs(iso) <= 0.01)}
# V2f anisotropy
an = {str(th): {"(z-x)/L": (p0(th, "z") - p0(th, "x")) / LOGCOEF,
                "slice-avg (xz-y)/L": (0.5 * (p0(th, "x") + p0(th, "z")) - p0(th, "y")) / LOGCOEF} for th in TH}
ch = (p0(0.0125, "z") - p0(0.0125, "x")) / (p0(0.05, "z") - p0(0.05, "x")) - 1
R["V2f"] = {"anisotropy": an, "change_0.05_to_0.0125": ch, "PASS": bool(abs(ch) < 0.30)}
# V2g regulator
K2 = (0.5, 0.9)
reg = {str(q): np.mean([piR(0.025, d, q, K2) for d in "xyz"]) / np.mean([piR(0.025, d, q) for d in "xyz"]) - 1 for q in (1, 3)}
La = (np.mean([p0(0.025, d, K2) for d in "xyz"]) - np.mean([p0(0.05, d, K2) for d in "xyz"])) / np.log(2)
Lb = (P0[0.025] - P0[0.05]) / np.log(2)
R["V2g"] = {"PiR_rel": reg, "L_rel": La / Lb - 1, "pit0_shift_K": np.mean([p0(0.025, d, K2) for d in "xyz"]) - P0[0.025],
            "PASS": bool(all(abs(v) <= 0.005 for v in reg.values()) and abs(La / Lb - 1) <= 0.02)}
# quadrature convergence
conv = {str(q): D["proj|0.0125|x|%.8g|80|6|0.6|1.05" % q][0] / g(0.0125, "x", q) - 1 for q in (0.05, 3)}
R["quadrature_n80_vs_n64"] = conv
json.dump(R, open(os.path.join(OUT, "vuoto_3d_exams.json"), "w"), indent=1, default=float)
print(json.dumps(R, indent=1, default=float))
