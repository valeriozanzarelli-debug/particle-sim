"""emiss_posthoc.py - analyses made AFTER seeing the production results (declared, not pre-registered; they do not
change any verdict of lab/out/emiss_results.json).  -> lab/out/emiss_posthoc.json

1. R1_48 failed (0.674): is it the coupling or the states?  The coupling dipole D(t) = <1s(t)|r|psi(t)> of the run
   SHRINKS in time (the 48^3 levels listened with S_T = 1500 are not single exact levels: valley multiplets, see
   docs/BATTIMENTO_1S.md on master 83485bc).  Here: local rate in windows vs the INSTANTANEOUS rule rate
   (4/3) alpha w^3 |D(t)|^2 / norm(t) W2.
2. S1 failed on the semiclassical side: 1 - p_2p(T) = 0.055.  Is it emission?  Energy given to the light by the
   semiclassical pure 2p vs the quantum one.
3. Lattice dipoles 1s -> 2p (48^3, listened levels) vs Dirac (continuum): angular ratios and size.
4. The Laboratory decay after the cure of the 1s (B2 of BATTIMENTO_1S: p_exc 0.097 at t = 450 -> 0.064 at t = 1600;
   semiclassical, K = 10, Gaussian average sigma = 2 cells, coupling as v2) vs the rates of this study.
"""
import json, os
import numpy as np
from zz_core import ALPHA
from emiss_core import Dirac3D, Light3D, CK, OUT, norm2
os.environ["EMISS_CFG"] = "none"
import emiss_configs as EC
from emiss_analysis import load

L48 = json.load(open(os.path.join(OUT, "emiss_levels_48.json")))
EX = json.load(open(os.path.join(OUT, "emiss_exact.json")))
res = {}
w = L48["2p3/2+3"]["omega"] - L48["1s+"]["omega"]
W2 = Light3D(128, 4.0, sponge=0, sigma=1.0).form_factor2(w)

# 1
win = []
for name in ("P_R1_48", "P_K10"):
    if not os.path.exists(os.path.join(CK, f"run_{name}.npz")):
        continue
    d = load(name)
    t, n, D = d["ser_t"], d["ser_norm"], d["ser_Dc"]
    K2 = 100.0 if name == "P_K10" else 1.0
    Ar = 4 / 3 * ALPHA * w**3 * np.sum(np.abs(D)**2, 1) / n * K2
    rows = []
    for a, b in ((100, 600), (600, 1200), (1200, 2000), (2000, 3000)):
        m = (t >= a) & (t < b)
        loc = -np.polyfit(t[m], np.log(n[m]), 1)[0]
        rows.append(dict(window=[a, b], A_local=loc, A_rule_inst_W2=float(np.mean(Ar[m]) * W2), ratio=float(loc / (np.mean(Ar[m]) * W2))))
        print(f"{name} t {a}-{b}: local rate {loc:.4e}, instantaneous rule rate x W2 {np.mean(Ar[m]) * W2:.4e}, ratio {loc / (np.mean(Ar[m]) * W2):.4f}")
    res[name + "_windows"] = rows
    res[name + "_D_abs_start_end"] = [np.abs(D[0]).tolist(), np.abs(D[-1]).tolist()]

# 2
dq, ds = load("P_R1_48"), load("P_S1sc")
res["S1_energy_to_light"] = dict(qe=float(dq["ser_eph"][-1] + dq["ser_abs"][-1]), sc=float(ds["ser_eph"][-1] + ds["ser_abs"][-1]))
print("energy given to the light by the pure 2p: quantum", res["S1_energy_to_light"]["qe"], " semiclassical", res["S1_energy_to_light"]["sc"])

# 3
D48 = Dirac3D(48, 0.5)
s = EC.level(48, "1s+")
R = EX["Z60"]["2p3/2"]["R_length"]
cont = {"2p1/2+": R**2 / 9, "2p3/2+": 2 * R**2 / 9, "2p3/2+3": R**2 / 3, "2p1/2-": 2 * R**2 / 9, "2p3/2-": R**2 / 9}
dips = {}
for st, c in cont.items():
    dd = D48.dip_matrix(s, EC.level(48, st))
    dips[st] = dict(abs=np.abs(dd).tolist(), sum2=float(np.sum(np.abs(dd)**2)), continuum_sum2=c, ratio=float(np.sum(np.abs(dd)**2) / c),
                    omega_level=L48[st.replace("-", "+")]["omega"] if st.endswith("-") and st.replace("-", "+") in L48 else L48.get(st, {}).get("omega"))
    print(f"{st}: |<1s|r|{st}>| = {np.round(np.abs(dd), 4)}, sum |d|^2 = {dips[st]['sum2']:.4f} vs continuum Dirac {c:.4f} (ratio {dips[st]['ratio']:.3f})")
res["dipoles_48"] = dips

# 4
pl0, pl1, T = 0.097, 0.064, 1600 - 450
A_lab = -(np.log(pl1 / (1 - pl1)) - np.log(pl0 / (1 - pl0))) / T
W2lab = Light3D(128, 4.0, sponge=0, sigma=2.0).form_factor2(w)
pre = 4 / 3 * ALPHA * w**3 * 100 * W2lab
res["lab_decay"] = dict(A_logistic=A_lab, W2_sigma2=W2lab,
                        pred_smooth_dev_states=pre * 1.70, pred_lattice_2p3_2_z=pre * dips["2p3/2+"]["sum2"],
                        pred_dirac_continuum=pre * cont["2p3/2+"])
print(f"Laboratory after the cure: logistic rate {A_lab:.3e}; predictions K^2 A_rule W2(sigma 2) with |D|^2 of: "
      f"lattice 1s->2p3/2 (z) {pre * dips['2p3/2+']['sum2']:.3e}, dev states (S_T 300) {pre * 1.70:.3e}, Dirac continuum {pre * cont['2p3/2+']:.3e}")
json.dump(res, open(os.path.join(OUT, "emiss_posthoc.json"), "w"), indent=1, default=float)

# 5. SC10 against the logistic with the RULE rate at t = 0 (K^2 A_rule W2, no fit) instead of the K10 fit
d = load("P_SC10")
t, p = d["ser_t"], np.array(d["ser_pexc"])[:, 0]
A0 = 100 * 4 / 3 * ALPHA * w**3 * dips["2p3/2+3"]["sum2"] * W2
pl = p[0] * np.exp(-A0 * t) / (1 - p[0] + p[0] * np.exp(-A0 * t))
res["SC10_rule"] = dict(A0=A0, max_dev=float(np.max(np.abs(p - pl))), p_end=float(p[-1]), logistic_end=float(pl[-1]))
print(f"SC10 vs logistic with K^2 A_rule(0) W2 = {A0:.3e} (no fit): max deviation {res['SC10_rule']['max_dev']:.4f}; p(T) {p[-1]:.4f} vs {pl[-1]:.4f}")
# 6. K10: bright/dark estimate. If the 2p is (bright fraction f, full coupling) + (dark, no coupling) the norm stalls at 1 - f
dk = load("P_K10")
nk = dk["ser_norm"]
res["K10_norm_end"] = float(nk[-1])
res["K10_D2_ratio_end"] = float(np.sum(np.abs(dk["ser_Dc"][-1])**2) / np.sum(np.abs(dk["ser_Dc"][0])**2))
print(f"K10: norm(T) {nk[-1]:.4f}; |D|^2(T)/|D|^2(0) {res['K10_D2_ratio_end']:.3f}")
json.dump(res, open(os.path.join(OUT, "emiss_posthoc.json"), "w"), indent=1, default=float)
