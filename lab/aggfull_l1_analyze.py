"""aggfull_l1_analyze.py - verdicts of PHASE 3 (lab/out/aggfull_l1_prereg.json) from lab/out/aggfull_l1.json.
Refits (same seeds, deterministic) the T and M rules of the helium-scale jobs for the direct <G_M - G_T> on 1s(1.6875).
Writes lab/out/aggfull_l1_verdict.json"""
import json, os
import numpy as np
import aggfull_l1 as L
import agg_zz_core as azc

OUT = L.OUT
R = json.load(open(os.path.join(OUT, "aggfull_l1.json")))
V = []


def check(n, ok, msg):
    V.append((n, bool(ok), msg)); print(f"[{'PASS' if ok else 'FAIL'}] {n}: {msg}", flush=True)


# M1
mm = [abs(v["m_over_mfine_minus_1"]) for k, v in R.items() if k.endswith(("_T_s1", "_T_s2", "_T_s3"))]
th = np.array([L.THETA[c] for c in ("A2", "A", "B")])
m = np.array([R[f"{c}_fine"]["m_fine"] * (1 + R[f"{c}_Z1_w2_T_s1"]["m_over_mfine_minus_1"]) for c in ("A2", "A", "B")])
mU = np.array([R[f"{c}_fine"]["m_fine"] * (1 + R[f"{c}_Z1_w2_U_s1"]["m_over_mfine_minus_1"]) for c in ("A2", "A", "B")])
Vd = np.stack([np.ones(3), th ** 2, th ** 4], 1)
m0, m0U = np.linalg.solve(Vd, m)[0], np.linalg.solve(Vd, mU)[0]
check("M1 mass pinned to the fine walk", max(mm) <= 1e-9 and abs(m0 - 1) <= 3e-6,
      f"max |m/m_fine - 1| over {len(mm)} T fits {max(mm):.1e}; m0 (theta -> 0) = {m0:.9f} (old learner, same fit: {m0U:.7f})")
# M2
z = [(s, R[f"A_Z10_w4_T_s{s}"]["check_Z10"]) for s in (1, 2, 3)]
check("M2 Z = 10 Pauli level within the 3e-5 budget (3 seeds)", all(c["ok_budget_3e-5"] for _, c in z),
      "; ".join(f"s{s}: energy {c['energy_rel']:.2e} rest {c['rest_rel']:.1e} packet {c['packet']:.1e} bz {c['bz_rel90']:.2f}" for s, c in z)
      + " | old learner U: " + ", ".join(f"{R[f'A_Z10_w4_U_s{s}']['check_Z10']['energy_rel']:.2e}" for s in (1, 2, 3)))
# K1
fr = [R[f"A_Z10_w4_T_s{s}"]["k4frac"] for s in (1, 2, 3)]
rhe = [R[f"A_Z10_w4_T_s{s}"]["R_band_He_Zeff1.6875_mEh"] / R[f"A_Z10_w4_T_s{s}"]["R_dirac_He_mEh"] for s in (1, 2, 3)]
check("K1 k^4 carried by the band", all(abs(f[x] - 1) <= 0.01 for f in fr for x in ("axis", "face", "body")) and all(abs(r - 1) <= 0.02 for r in rhe),
      f"axis/face/body {[round(f['axis'], 5) for f in fr]} {[round(f['face'], 5) for f in fr]} {[round(f['body'], 5) for f in fr]}; "
      f"band/Dirac on 1s(1.6875) {[round(r, 4) for r in rhe]} (old level-1 learner: {R['A_Z1_w2_U_s1']['R_band_He_Zeff1.6875_mEh'] / R['A_Z1_w2_U_s1']['R_dirac_He_mEh']:.3f})")
# P1
th_A = L.THETA["A"]; prev = azc.Level("fine3", 3, "dirac4", 4, 1, th_A, period=2)
cf = np.array(R["A_fine"]["c2c4"]); D = L.dirs13(); h_lc = 16 * th_A
dG = []
for s in (1, 2, 3):
    kap, Gs, sw, phi0 = L.samples(prev, 16, 2 / L.A0 * th_A * 16, np.random.default_rng(s), 8000, 0.25, 1e-3)
    gT, oT = L.band_of(L.fit(kap, Gs, sw, phi0, 4, "T", cons=cf, D=D))
    gM, oM = L.band_of(L.fit(kap, Gs, sw, phi0, 4, "M", cons=cf, D=D))
    k = azc.sample_state_k(np.random.default_rng(9), 200000, 3, 1.6875)
    kp = k * h_lc / L.A0
    dG.append(float(np.mean(L.gband(gM, oM, kp) - L.gband(gT, oT, kp)) / h_lc * L.A0 ** 2 * 1e3))
pm = [R[f"A_Z2_w4_M_s{s}"] for s in (1, 2, 3)]
negmax = max(max(p["kernel_neg"].values()) for p in pm)
check("P1 Markov level: positive kernel at every beta", all(p["wrong_sign_taps_tol"] == 0 for p in pm) and negmax <= 1e-12
      and all(p["max_c2_rel_dev"] <= 1e-6 and p["max_c4_rel_dev"] <= 1e-6 for p in pm) and all(abs(x) <= 0.01 for x in dG),
      f"wrong taps {[p['wrong_sign_taps_tol'] for p in pm]}; max negative weight beta = 1e-6..1e-2: {negmax:.1e}; c2/c4 dev "
      f"{max(p['max_c2_rel_dev'] for p in pm):.1e}/{max(p['max_c4_rel_dev'] for p in pm):.1e}; <G_M - G_T> on 1s(1.6875) "
      f"{[round(x, 5) for x in dG]} mEh; (Z = 10 check of M, reported: {[round(p['check_Z10']['energy_rel'], 6) for p in pm]}, "
      f"bz {[round(p['check_Z2']['bz_rel90'], 2) for p in pm]})")
# P2
dev = {w: max(max(R[f"A_Z2_w{w}_M_s{s}"]["max_c2_rel_dev"], R[f"A_Z2_w{w}_M_s{s}"]["max_c4_rel_dev"]) for s in (1, 2, 3)) for w in (2, 3, 4)}
check("P2 hop-length bound w >= sqrt(3)/h = 3.24", dev[2] > 1e-3 and dev[3] > 1e-3 and dev[4] <= 1e-6,
      f"max constraint deviation w=2 {dev[2]:.2e}, w=3 {dev[3]:.2e}, w=4 {dev[4]:.1e}")
# PA1
odds = {k: max(v["parity_odd_rel"].values()) for k, v in R.items() if isinstance(v, dict) and "parity_odd_rel" in v}
oTM = max(o for k, o in odds.items() if "_T_" in k or "_M_" in k)
oU = {k: round(o, 7) for k, o in odds.items() if "_U_" in k}
check("PA1 parity exact in T and M", oTM <= 1e-12, f"max odd/even over all T and M fits {oTM:.1e}; old learner U (raw): {oU}")
bm = {}
for k in ("A_Z1_w2_U_s1", "A_Z10_w4_U_s1", "A_Z10_w4_T_s1", "A_Z2_w4_M_s1"):
    kn = R[k]["kernel_neg"]
    pos = [float(b) for b, v in kn.items() if v <= 1e-12]
    bm[k] = min([b for b in sorted(map(float, kn)) if all(kn[x] <= 1e-12 for x in kn if float(x) >= b)], default=None)
print("beta_min (a.u.) by level:", bm)
json.dump(dict(verdict=V, beta_min=bm, dG_M_T_mEh=dG, m0=m0, m0_old=m0U), open(os.path.join(OUT, "aggfull_l1_verdict.json"), "w"), indent=1, default=float)
print(f"=== {sum(v[1] for v in V)} PASS / {sum(not v[1] for v in V)} FAIL")
