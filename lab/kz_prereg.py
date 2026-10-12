"""kz_prereg.py - PRE-REGISTRATION of K1-K3 ('sincronia': metronomes that synchronise, Kibble-Zurek defects).
Writes lab/out/kz_prereg.json.  Committed BEFORE the production runs (kz_k1.py, kz_k2.py, kz_k3.py).

DEV (declared, looked at, NOT production numbers):
 - timing (torch CPU, 2 threads): Metropolis sweep 7.6 ms (L=32), 15 ms (L=48), 44 ms (L=64); acceptance 0.83 at Tc.
   Swendsen-Wang + over-relaxation: 14 ms per sweep at L=32.
 - one seed, L=32: counting the RAW plaquette windings at t = t_hat = tauQ^0.573 after Tc gives 0.60 / 0.56 / 0.52 / 0.49
   per site for tauQ = 32 / 128 / 512 / 2048: almost flat, because just below Tc the field is full of small THERMAL
   vortex loops (lattice-size, equilibrium).  Gradient-flow cooling for a FIXED time (2 or 8) removes them but then
   the cooling time does not scale with t_hat.  Fix adopted: cool for a time b * t_hat (for z = 2 dynamics, t_hat ~
   xi_hat^2, so a cooling time proportional to t_hat smooths up to a length proportional to xi_hat: the count is a
   pure scaling quantity n = xi_hat^-2 F(a, b) if KZ holds).
 - with that protocol, one seed, L = 32, tauQ = 64 and 1024, a = 1,2,4, b = 0.25,0.5,1: densities 0.088->0.015,
   0.046->0.0085, 0.026->0.0047 (a=1), ... i.e. a two-point slope around 0.6-0.65 (one seed, L=32 already small at
   tauQ = 1024: mean line spacing ~ 15).  This is the only alpha-like number seen before the pre-registration.
"""
import json, os, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

PRE = {
    "model": "classical XY model H = -sum cos(theta_i - theta_j), J = 1, periodic simple cubic lattice (2D square for K3). "
             "Tc(3D) = 2.20184 (literature beta_c = 0.454165), used as the input crossing point of the quench; K2 checks it.",
    "theory": {
        "KZ": "n ~ tauQ^-alpha, alpha = (D - d) nu / (1 + z nu); D = 3, line defects d = 1.",
        "mean_field": "nu = 1/2, z = 2: alpha = 0.5",
        "3D_XY_modelA": "nu = 0.6717, z = 2: alpha = 2*0.6717/(1+1.3434) = 0.5733 (the honest target: our dynamics is model A)",
        "3D_XY_modelF_helium": "z = 3/2: alpha = 0.668 (NOT what our dynamics can test)",
        "freeze_out": "t_hat = tauQ^(z nu/(1+z nu)) = tauQ^0.573 (time unit 1 sweep, prefactor 1 by convention)",
    },
    "K1": {
        "dynamics": "checkerboard Metropolis (proposal +-1 rad), linear quench T(t) = Tc (1 - t/tauQ), start at 1.2 Tc after "
                    "200 sweeps of equilibration there; t = 0 at Tc crossing.",
        "count": "vortex line length per site = sum over all plaquettes of |winding| / L^3, on a COPY of the field cooled by "
                 "gradient flow (T = 0) for b * t_hat; counted at t = a * t_hat after Tc.",
        "main_protocol": {"L": 48, "seeds": 8, "tauQ": [64, 128, 256, 512, 1024, 2048, 4096], "a": 2, "b": 0.5,
                          "t_hat_exponent": 0.573},
        "fit": "least squares of ln n vs ln tauQ over all 7 tauQ (seed means); 95% CI by bootstrap over seeds "
               "(resampling the 8 seeds independently at each tauQ, 4000 resamples, percentile interval).",
        "PASS_prediction_beyond_mean_field": "95% CI of alpha (main protocol) EXCLUDES 0.5 AND INCLUDES 0.5733.",
        "sensitivity_reported": [
            "a in {1, 2, 4} x b in {0.25, 0.5, 1}: 9 alphas (same runs)",
            "L = 32 (8 seeds, same tauQ) and L = 64 (4 seeds, tauQ 256-4096) finite-size check",
            "fit window tauQ >= 256 only",
            "mean-field protocol: t_hat = tauQ^0.5 (a = 2, b = 0.5), L = 48, 6 seeds",
        ],
        "ROBUST_flag": "in addition to PASS: all 9 (a,b) alphas within 0.05 of the main alpha, AND the L = 64 alpha "
                       "within 0.05 of the L = 48 alpha on the same tauQ window, AND the mean-field-protocol alpha "
                       "also excludes 0.5. If PASS but not ROBUST, the claim is downgraded to 'indicative'.",
    },
    "K2": {
        "method": "equilibrium Swendsen-Wang (embedded Ising, Wolff reflection) + 1 over-relaxation per sweep, L = 8, 12, "
                  "16, 24, 32 at beta = 1/2.20184; 1000 thermalisation sweeps, then measurements (L=8-16: 20000, "
                  "L=24: 12000, L=32: 8000).  Binder ratio U4 = <m^4>/<m^2>^2; derivatives d ln<|m|>/d beta and "
                  "d ln<m^2>/d beta by fluctuation formulas ~ L^(1/nu); U4(T) by single-histogram reweighting.",
        "nu_estimate": "mean of the two log-derivative power-law fits over L = 8..32; error = jackknife (20 blocks) "
                       "combined with half the difference between the two estimators and between fit windows "
                       "8..32 vs 12..32.",
        "PASS": "|nu - 0.6717| <= 0.05 AND the U4 crossing temperature of (16, 32) within 0.5% of 2.20184. "
                "Report also how far from mean field (0.5) and the crossing value of U4 (literature ~1.243).",
    },
    "K3": {
        "a_pair": "2D, L = 128, vortex at x0 - d/2, antivortex at x0 + d/2, d = 24, arctan ansatz relaxed by gradient flow "
                  "for t = 3; then inertial dynamics theta'' = -gamma theta' + sum sin(theta_j - theta_i), gamma = 0, "
                  "Verlet dt = 0.05. Lattice wave speed c = 1 (omega^2 = sum 2(1 - cos k)).",
        "a_PASS": "(i) the pair approaches and annihilates (vortex count -> 0) before t = 300; (ii) total energy conserved "
                  "to 1e-3 relative; (iii) maximum vortex speed < c = 1; (iv) after annihilation the energy (>= 90% of "
                  "the initial pair energy) leaves as a wave whose front moves at 1.0 +- 0.1.",
        "b_kick": "same pair given a Lorentz-boosted initial velocity v = 0.5 perpendicular to the separation "
                  "(theta(x, y - v t) with gamma factor, theta' = -v gamma d_y theta): report the measured initial "
                  "centroid speed (expected ~0.5) and whether it stays < 1. Report only (no pass/fail).",
        "c_KZ_2D": "2D quench L = 128, tauQ = 256, stop at T = 0.2 (t = 0.909 tauQ), then inertial dynamics gamma = 0.02 "
                   "from the frozen configuration (momenta 0): vortex number vs time. PASS: vortex number drops by >= 50% "
                   "in t = 200 (defects move and annihilate) while kinetic energy (waves) appears.",
        "d_KZ_3D": "a K1 final 3D configuration (L = 48, tauQ = 1024, cooled) under inertial dynamics gamma = 0.02: line "
                   "length vs time (report only).",
    },
    "K4": "docs/SINCRONIA.md (simple Italian) + lab/out/kz_figure.png",
}

if __name__ == "__main__":
    try:
        sha = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=HERE).decode().strip()
    except Exception:
        sha = "?"
    PRE["written"] = datetime.datetime.now().isoformat(timespec="seconds")
    PRE["parent_commit"] = sha
    with open(os.path.join(OUT, "kz_prereg.json"), "w", encoding="utf-8") as f:
        json.dump(PRE, f, indent=1)
    print("written", os.path.join(OUT, "kz_prereg.json"))
