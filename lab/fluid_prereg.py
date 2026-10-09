"""fluid_prereg.py - PRE-REGISTRATION of F1-F4 ('fluido in negativo'). Writes lab/out/fluid_prereg.json.
Committed BEFORE the production runs (fluid_f1.py, fluid_f2.py, fluid_f3.py, fluid_f4.py).

DEV (declared, looked at, NOT production numbers):
 F1  - speed: N = 2567 balls, ~1 ms/step. Antithetic kicks (+dV/-dV from the same state) on R_w = 6, L = 60, M0 = rho pi R_w^2:
       the disk velocity decays monotonically, NO plateau (2D: viscous entrainment never stops). 30 replicas, t <= 4:
       P/V(t) passes M0 + N_miss (~189) at t ~ 0.8. Weak harmonic trap fits (omega 0.3 and 0.6) gave a strongly damped,
       non-oscillator-like response (M from Omega^2 + Gamma^2 = 155 and 237: frequency-dependent friction, the fit model is
       wrong) -> trap method dropped. Laplace-domain mass M(s) = P/(s Vhat(s)) on 24 replicas, t <= 8: monotone in s
       (512 at s = 0.1 ... 103 at s = 5), no plateau; at s = c/R_eff (0.94) M - M0 = 81 vs N_miss = 104 (ratio 0.78);
       ratio 1 is reached near s ~ 0.75 c/R_eff. Window estimate: viscous correction ~ sqrt(nu/s)/R, compressible ~ s R/c,
       both >~ 25 % at R ~ 6: the 2D MD fluid has no clean added-mass window at these sizes. Criteria below set
       accordingly (bracket + point at the NATURAL acoustic rate s* = c/R_eff, not at the s where dev gave ratio 1).
       Sound speed (thermodynamic route, rho ~ 0.75, T = 1): c_T = 4.8, c_S = 6.3.
 F2  - radial GP solver: slope of E(ln R) = 0.9993 pi; constant A = 2.078 (Fetter 1.464 with xi = hbar/sqrt(2 m g n),
       i.e. 1.464 sqrt 2 = 2.070 in our xi = hbar/sqrt(m g n)); N_miss(R) also grows like pi ln R (seen in dev; that is
       why E2c below exists - it was suggested by the dev, declared here).
 F3  - h = 0.25: w gamma = 0.998 ... 0.983 (v0 0 ... 0.8), E/gamma 7.993 ... 8.090; v0 = 0.9 starts to break (0.955).
       h = 1: kinks radiate and stall at v ~ 0.4. Kink-antikink: h = 0.25 pass through at all v tried; h = 1 captured
       into a long-lived breather at v <= 0.4 (centre energy 11.6 -> 11.4 over 70 time units).
"""
import json, os, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
S0 = 20261010

PRE = {
    "F1": {
        "model": "2D MD, WCA balls (m = sigma = eps = 1, cut 2^(1/6)), periodic L x L, triangular start at rho0 = 0.75, "
                 "T = 1 (velocity rescale every 10 steps during equilibration 2000 steps and during waits; NVE during "
                 "kicks), dt = 0.005, Verlet list skin 0.4. Hole = free disk of radius R_w whose wall is a WCA in "
                 "s = r - R_w (smooth, radial force only = slip wall), bare mass M0.",
        "sizes": {"R_w": [3, 6, 10], "L": [40, 60, 90], "replicas": [40, 24, 10],
                  "seeds": [S0 + 1, S0 + 2, S0 + 3]},
        "kick": "from an equilibrium state: two copies, disk V += +dV x and -dV x (dV = 1), response "
                "V(t) = (V+ - V-)/2, t_max = L / c_S; replicas separated by 400 thermostatted steps",
        "M0": "M0 = rho0 pi R_w^2 (nominal, formula fixed in advance); bare-mass check: R_w = 6 with M0 = 0.25 rho0 pi R_w^2 "
              "and dV = 2, 24 replicas, seed %d" % (S0 + 4),
        "c_S": "thermodynamic route on pure fluid L = 32: c_T^2 = dP/drho|T (rho0 * 0.97 / 1.03), "
               "c_S^2 = c_T^2 + T (dP/dT|rho)^2 / (rho^2 c_v) (T = 0.95 / 1.05), 1500 + 4000 steps each, virial pressure",
        "N_miss": "radial profile around the disk during waits; rho_b = mean density in r in [R_w + 8, L/2 - 1]; "
                  "N_miss = rho_b pi (R_w + 6)^2 - <count in r < R_w + 6>; R_eff = sqrt(N_miss / (pi rho_b))",
        "estimator": "Laplace-domain fluid mass M_fl(s) = P / (s Vhat(s)) - M0, Vhat(s) = int_0^tmax V e^{-st} dt "
                     "(trapezoid), s* = c_S / R_eff",
        "exams": {
            "E1_bracket": "m N_miss in [M_fl(2 c/R_eff), M_fl(c/(2 R_eff))] for all 3 sizes",
            "E2_point": "M_fl(s*) / (m N_miss) in [0.75, 1.25] for all 3 sizes",
            "E3_scaling": "log-log slope of M_fl(s*) vs N_miss over the 3 sizes in [0.85, 1.15]",
            "E4_bare_mass": "R_w = 6: |M_fl(s*; M0) - M_fl(s*; M0/4)| / M_fl(s*; M0) < 0.10",
        },
        "prediction": "ideal 2D potential flow: added mass of a cylinder = rho pi R^2 = mass of the missing balls (ratio 1). "
                      "Expected deviation in MD: viscous (raises M at low s) and compressible (lowers at high s).",
    },
    "F2": {
        "model": "radial Gross-Pitaevskii vortex, psi = f(r) e^{i theta}, hbar = m = g n = 1 (xi = 1, c = 1, n = 1); "
                 "Newton on r in (0, 400], dr = 0.02, f(400) = 1",
        "R": [10, 20, 40, 80, 160],
        "exams": {
            "E2a_log": "slope of E(R) vs ln R = pi (pi n hbar^2/m) within 1 %",
            "E2b_core": "A in E = pi ln(A R / xi): within 2 % of 1.464 sqrt(2) = 2.070 (Fetter 2009 convention converted)",
            "E2c_missing_balls": "slope of [E_swirl(R) - m c^2 N_miss(R)] vs ln R: |.| < 0.02 pi (swirl energy = missing "
                                 "balls x c^2 up to an R-independent constant: Bernoulli). Constant reported.",
        },
        "vortex_mass": "NOT measured dynamically (debt). Reported: candidate masses E/c^2 (Popov / Duan, log-divergent), "
                       "core missing mass, literature disagreement stated.",
    },
    "F3": {
        "model": "Frenkel-Kontorova in sine-Gordon units: phi_n'' = (phi_{n+1} - 2 phi_n + phi_{n-1})/h^2 - sin phi_n; "
                 "free ends; c = 1 = sound speed of the chain without substrate; velocity Verlet dt = 0.2 h",
        "h_main": 0.25, "h_secondary": [0.5, 1.0],
        "single_kink": "start = continuum moving kink 4 atan exp(gamma (x - x0)), Lx = 400, T = 60; fit "
                       "4 atan exp((x-X)/w) every time unit for t > 10; v = slope of X(t); E = energy within |x - X| < 15",
        "v0": [0, 0.2, 0.4, 0.6, 0.8, 0.9],
        "exams": {
            "E3a_contraction": "h = 0.25: |w gamma / w(0) - 1| < 0.05 for every v0 <= 0.8",
            "E3b_energy": "h = 0.25: |E / (gamma E(0)) - 1| < 0.05 for every v0 <= 0.8",
            "E3c_limit": "h = 0.25 driven-damped phi'' += F - alpha phi', alpha = 0.05, F in {0.05, 0.1, 0.2, 0.4, 0.8}, "
                         "kink from rest, T = 150, v from last 50: all v < 1, v(F=0.8) > 0.95, and |v / v_MS - 1| < 0.02 "
                         "with McLaughlin-Scott v_MS = 1/sqrt(1 + (4 alpha / (pi F))^2)",
            "E3d_annihilation": "kink at -10 + antikink at +10 moving towards each other; h in {0.25, 1.0}; "
                                "v in {0.05, 0.1, 0.2, 0.4, 0.8}; run d/v + 1000. Final: lumps = greedy windows of "
                                "+-10 with energy > 2; frac_rad = 1 - E_lumps/E. Class: 'annihilated into phonons' if "
                                "frac_rad > 0.5; 'bound (breather)' if one lump at |x| < 10; 'pass-through' if two lumps. "
                                "PREDICTION: h = 0.25 pass-through at all v (near-integrable sine-Gordon); h = 1 bound "
                                "at v <= 0.4; NO case annihilates within the run. Breather decay rate reported.",
        },
    },
    "F4": {
        "what": "draw only the defects: F1 = isoline rho_cg = rho_b / 2 of the Gaussian-smoothed density (width 1.5, grid "
                "0.5) of one snapshot of R_w = 10; F2 = isoline |psi|^2 = 1/2 (one circle); F3 = one marker per kink "
                "(phi crossing pi + 2 pi k)",
        "exams": {
            "E4a_compression": "F1: drawn primitives (isoline segments) / balls < 0.01",
            "E4b_fidelity": "F1: centroid of the deficit field within 1 sigma of the true disk centre, and deficit count "
                            "from the field within 20 % of N_miss (profile)",
        },
    },
    "parameters": {
        "F1": "physical: rho0, T (2) + sizes R_w (the variable); numerical: dt, skin, dV, waits, r-window for N_miss (+6)",
        "F2": "0 free (units hbar = m = g n = 1)",
        "F3": "h (1 physical: discreteness); alpha, F for E3c only",
        "F4": "smoothing width, grid, iso level (3 rendering choices)",
    },
}


def meta():
    try:
        git = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    except Exception:
        git = "unknown"
    return {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": git}


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "fluid_prereg.json"), "w") as f:
        json.dump({"meta": meta(), "prereg": PRE}, f, indent=1)
    print("written lab/out/fluid_prereg.json")
