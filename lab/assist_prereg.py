"""assist_prereg.py - pre-registration of the dynamically assisted Schwinger study (writes lab/out/assist_prereg.json).

Written BEFORE the production runs. Contains: configurations, the theory predictions computed now (worldline
instanton / WKB / exact Sauter-pulse formula / LCFA), the PASS/FAIL exams with thresholds, free parameters (0).
Dev runs done before this file (declared, NOT production): dev1 single Sauter pulse (E, tau) = (0.15, 22.2),
(0.05, 2), (0.3, 3); dev2/4 assisted E1 = 0.125, eps = 0.05; dev5/6 continuum-only slope fits at eps = 0.01
(Sauter, cos) and eps = 0.001 (cos); dev7/8 continuum Floquet at xi = 0.05, 0.1, 0.2, 0.3.
What dev taught (and shaped the exams): (i) the finite-eps WKB with the turning point next to the pole of the
Sauter pulse is NOT asymptotic at reachable E1 (local adiabaticity c = eps m^2/(gamma^2 E1) < 1), so the threshold
exam uses the eps -> 0 response Q = [f(eps) + f(-eps) - 2 f(0)] / (2 eps^2), whose exponent is the eps -> 0
worldline instanton action; (ii) Floquet splittings below ~1e-8 are at the numerical floor -> xi = 0.1, 0.2.
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from assist_core import (S_wi, s_wi_sauter, s_linder_cos, s_sauter_single, sauter_pulse_exact, wkb_pairs, Field,
                         qbar, OUT, PI)
from scipy.optimize import brentq

GS = 0.2                                   # slow-field Keldysh parameter m Om / E1 (fixed)
E1_FIT = [0.16, 0.125, 0.1, 0.08, 0.0625]  # E1 / E_S for the exponent fits (5 points)
GAM_S = [0.5, 1.0, 1.3, 1.5, 1.6, 1.7, 1.8, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0]
GAM_C = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0]
THETAS = [0.2, 0.1, 0.05, 0.025]

pre = dict(study="assist_ (Schwinger assistito dinamicamente)", date="2026-10-13",
           units="hbar = c = m = 1, e absorbed: E in units of E_S = m^2/e; lattice spacing a = theta / m",
           free_parameters=0,
           inputs="m (through theta = m a), field shapes and amplitudes (the experiment being simulated)",
           rule="lab/assist_core.py walk_pairs: zig-zag 1+1 D (R1 method), mode by mode, Peierls phase at mid-tick",
           reference="lab/assist_core.py dirac_pairs: continuum 1+1 Dirac, 4th-order Magnus, dt = 0.01 (0.005 check)",
           richardson="f_extrap = (4 f(theta=0.025) - f(theta=0.05)) / 3  (assumes O(theta^2), checked by the 0.1/0.05 ratio)")

# ---------------------------------------------------------------- V1 Schwinger (slow Sauter pulse, gamma_s = 0.05)
v1 = dict(configs=[dict(E=E, gs=0.05) for E in (0.25, 0.16, 0.1)], thetas=[0.05, 0.025, 0.0125],
          p_grid="p in [-16, 16], dp = 0.02",
          exams=dict(V1a="max_p |f_lat(theta=0.0125)/f_exact - 1| <= 0.03 over modes with f_exact > 1e-3 f_exact(0), "
                          "and Richardson(0.025, 0.0125) <= 0.01",
                     V1b="|N_lat(theta=0.0125) / N_LCFA - 1| <= 0.10, N_LCFA = int dt E(t)/(2 pi) exp(-pi/E(t)) "
                          "(local Schwinger-Nikishov rate)"))
for c in v1["configs"]:
    E, gs = c["E"], c["gs"]; Om = gs * E
    t = np.linspace(-12 / Om, 12 / Om, 400001)
    Et = E / np.cosh(Om * t) ** 2
    c["N_LCFA"] = float(np.trapezoid(Et / (2 * PI) * np.exp(-PI / np.maximum(Et, 1e-300)), t))
    p = np.arange(-16, 16.0001, 0.02)
    c["N_exact"] = float(np.trapezoid(sauter_pulse_exact(p, E, 1 / Om), p) / (2 * PI))
    c["f0_exact"] = float(sauter_pulse_exact(0.0, E, 1 / Om))
    c["f0_schwinger"] = float(np.exp(-PI / E))
pre["V1_schwinger"] = v1

# ---------------------------------------------------------------- V2 Keldysh crossover (single Sauter pulse)
v2 = dict(configs=[dict(gamma=g, E=E) for g, E in ((0.1, 0.2), (0.3, 0.2), (1.0, 0.2), (3.0, 0.1), (10.0, 0.05))],
          thetas=[0.1, 0.05, 0.025], p_grid="p in [-(1/gamma + 3), 1/gamma + 3], dp = 0.01",
          exams=dict(V2="for every gamma: |f_lat(0)/f_exact(0) - 1| <= 0.03 and |N_lat/N_exact - 1| <= 0.03 at "
                        "theta = 0.025; Richardson(0.05, 0.025) within 0.01"))
for c in v2["configs"]:
    tau = 1 / (c["gamma"] * c["E"])
    c["f0_exact"] = float(sauter_pulse_exact(0.0, c["E"], tau))
    c["f0_schwinger"] = float(np.exp(-PI / c["E"]))
pre["V2_keldysh"] = v2

# ---------------------------------------------------------------- V3 multiphoton (oscillating field, Floquet, p = 0)
v3 = dict(field="A(t) = xi m sin(omega t), p = 0; resonances n = 1..5 located as minima of the Floquet splitting",
          xis=[0.1, 0.2], walk_ticks_per_period="N even, a ~ 0.025 (also a ~ 0.05)", dirac_steps=3000,
          exams=dict(V3a="selection rule at p = 0: Delta_n(even) / Delta_{n+1}(odd) <= 1e-2 for n = 2, 4 (both xi), walk",
                     V3b="multiphoton law Delta_n ~ xi^n: log2(Delta_n(0.2)/Delta_n(0.1)) within 5% of n, n = 1, 3, 5, walk",
                     V3c="|omega_n(walk, a~0.025) / omega_n(Dirac) - 1| <= 1e-3 for n = 1..5"),
          report="omega_n vs the adiabatic dressed-energy estimate n omega = 2 qbar (not graded)",
          qbar={str(x): qbar(x) for x in (0.1, 0.2)})
pre["V3_multiphoton"] = v3

# ---------------------------------------------------------------- A: assisted, weak Sauter pulse
Ssch = S_wi(0.1, GS)   # below threshold = single slow Sauter pulse exponent (x E1)
gc_pred = brentq(lambda g: S_wi(g, GS) - Ssch * (1 - 1e-9), 1.0, 3.0)
A = dict(field="E(t) = E1 sech^2(gs E1 t) + eps E1 sech^2(gamma E1 t)  (gamma = m omega/E1 combined Keldysh)",
         gs=GS, gammas=GAM_S, E1_fit=E1_FIT, eps_response=0.01, eps_finite=[0.01, 0.1], thetas=THETAS,
         observables=dict(
             f0="pair number in the mode p = 0 (crosses the gap at the pulse peak)",
             Q="eps -> 0 response Q = [f(+eps) + f(-eps) - 2 f(0)] / (2 eps^2), eps = 0.01 (0.005 check)",
             S_fit="fit ln Q = -S/E1 + nu ln(1/E1) + c over E1_fit (3 parameters, 5 points)",
             N="int dp/2pi f(p), p in [-4, 4], dp = 0.01, at E1 = 0.1; with E1 only, eps pulse only, both"),
         predictions=dict(S_wi={str(g): S_wi(g, GS) for g in GAM_S},
                          S_wi_gs0={str(g): float(PI * s_wi_sauter(g)) for g in GAM_S},
                          gamma_crit_gs0=PI / 2, gamma_crit_gs=gc_pred,
                          note="S_wi = eps->0 worldline instanton (Schuetzhold-Gies-Dunne 2008) with the slow field "
                               "sech^2(Om t) kept: instanton cut at tau = pi/(2 omega); gs -> 0 gives eq. 15"),
         exams=dict(
             A1="rule vs continuum: |f_lat(theta=0.025)/f_Dirac - 1| <= 0.03 and Richardson <= 0.01 for f(0) at every "
                "(gamma, E1, eps in {0, +-0.01, 0.1}); and |Q_lat/Q_Dirac - 1| <= 0.05 (Q from Richardson values)",
             A2="exponent: |S_fit(rule)/S_wi(gamma, gs) - 1| <= 0.03 for gamma <= 6 (gamma = 10, omega = m: reported)",
             A3="threshold: one-parameter fit S(gamma) = S_wi(gamma * gc_pred/gc, gs) to S_fit(rule), gamma <= 6 -> "
                "|gc/gc_pred - 1| <= 0.05 (gc_pred = %.4f; = pi/2 for gs -> 0)" % gc_pred,
             A4="report only (prediction of the model): f(0), N, enhancement N_both/(N_E1 + N_eps) at E1 = 0.1, "
                "eps = 0.01, 0.1, all gamma and gamma = 20 (omega = 2m); compared with WKB (finite eps) and with "
                "zeroth+first order in eps (Torgrimsson et al. 2017: should suffice for Sauter)"))
pre["A_sauter"] = A

# ---------------------------------------------------------------- O: assisted, weak oscillating field
from scipy.special import iv
O = dict(field="A(t) = (E1/Om) tanh(Om t) + (eps E1/omega) sin(omega t) sech^2(Om t), Om = gs E1",
         gs=GS, gammas=GAM_C, E1_fit=E1_FIT, eps=[0.001, 0.01], thetas=THETAS,
         observables=dict(S_fit="fit ln f(0) = -S/E1 + nu ln(1/E1) + c over E1_fit; S0_fit same with eps = 0"),
         predictions=dict(), exams=dict(
             O1="rule vs continuum: |f_lat(0.025)/f_Dirac - 1| <= 0.03, Richardson <= 0.01, all configs",
             O2="where the first-order correction 2 eps I1(gamma)/gamma <= 0.3: |S_fit/S_WKB - 1| <= 0.03 "
                "(S_WKB = exact-field WKB exponent at p = 0, dominant turning point)",
             O3="where 2 eps I1(gamma)/gamma <= 0.05 (first order trustworthy): |(S_fit/S0_fit) / s_Linder - 1| <= 0.03 "
                "(Linder et al. 2015 eq. 37)",
             O4="report only: points beyond (gamma ~ |ln eps|, nonperturbative / omega ~ m)"))
for eps in O["eps"]:
    for g in GAM_C:
        F = Field(0.1, GS * 0.1, eps * 0.1, g * 0.1, "cos")
        F0 = Field(0.1, GS * 0.1)
        f, al = wkb_pairs(F, 0.0, return_all=True)
        f0, al0 = wkb_pairs(F0, 0.0, return_all=True)
        O["predictions"][f"eps{eps}_g{g}"] = dict(
            first_order_corr=float(2 * eps * iv(1, g) / g), graded_O2=bool(2 * eps * iv(1, g) / g <= 0.3),
            graded_O3=bool(2 * eps * iv(1, g) / g <= 0.05),
            S_WKB=float(4 * al[0][0] * 0.1), S0_WKB=float(4 * al0[0][0] * 0.1),
            s_ratio_WKB=float(al[0][0] / al0[0][0]), s_Linder=float(s_linder_cos(g, eps)), n_turning_points=len(al))
pre["O_cos"] = O

pre["limits_report"] = ("1+1 vs 3+1: via m -> m_perp identity (spin factor 2, integral over p_perp) for one config; "
                        "real lasers: table of E/E_S, omega/m, gamma and exponent for ELI-class + XFEL (not graded)")
pre["amendments"] = "any change after this commit goes to lab/out/assist_amend.json, declared in the document"

json.dump(pre, open(os.path.join(OUT, "assist_prereg.json"), "w"), indent=1)
print(json.dumps(pre["A_sauter"]["predictions"], indent=1))
print("gamma_crit(gs=0.2) =", gc_pred)
print(json.dumps(pre["V1_schwinger"]["configs"], indent=1))
