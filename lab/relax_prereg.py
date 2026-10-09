"""relax_prereg.py - PRE-REGISTRATION of the 'Born rule as equilibrium' tests (written BEFORE any run).

Model: de Broglie-Bohm surfers (hbar = m = 1), v = Im(grad psi / psi).
Writes lab/out/relax_prereg.json with git hash and time. Committed before relax_box.py / relax_slit.py run.
"""
import json, os, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

PREREG = {
    "test1_box": {
        "setup": "2-D box [0,pi]^2, psi = sum of the M lowest (n,m) in a square block n,m = 1..sqrt(M), equal "
                 "amplitudes, random phases (seed fixed), energies (n^2+m^2)/2; M = 1, 4, 16, 25; N = 20000 surfers; "
                 "RK4 with per-surfer sub-cycling near nodes; t in [0, 8 pi]",
        "starts": "rho0 = |psi0|^4 normalised, and rho0 = uniform; control: rho0 = |psi0|^2 (M = 16)",
        "H": "Hbar = sum_cells rhobar ln(rhobar / |psi|^2bar), cells 16x16 (main), also 4x4, 8x8, 32x32; "
             "finite-sample floor F = (K-1)/(2N)",
        "tau": "tau = first time Hbar - F <= (Hbar(0) - F)/e (linear interpolation); rate = 1/tau",
        "vortices": "count of phase-winding plaquettes of psi on a 128^2 grid, time average over 20 snapshots",
        "criteria": {
            "C1_control_M1": "M = 1: Hbar(t_end) / Hbar(0) > 0.9 for both starts (no motion, no relaxation)",
            "C2_equivariance": "M = 16 started at |psi|^2: max_t Hbar <= 3 F",
            "C3_relax": "M = 16 and M = 25, both starts: Hbar(t_end) - F < 0.1 (Hbar(0) - F)",
            "C4_order": "rate(4) < rate(16) <= rate(25)*1.2 i.e. more modes -> not slower; and vortex count "
                        "increases with M; rank of rate equals rank of vortex count over M = 4, 16, 25",
            "C5_exponential": "ln(Hbar - F) vs t linear fit over [0, 3 tau] has R^2 >= 0.9 (M = 16, 25) - "
                              "'roughly exponential'",
        },
        "expected_failure_risk": "M = 4 may relax only partially (few vortices, quasi-periodic): reported, not hidden",
    },
    "test2_scale_lag": {
        "2a": "same runs as test1 (M = 16, 25, start |psi0|^4): tau at coarse-graining 4x4 vs 16x16 cells. "
              "Valentini-like prediction: tau(4x4) > 1.2 tau(16x16) (large scales lag)",
        "2b": "M = 16, start rho0 = |psi0|^2 (1 + 0.8 cos(k x)), long k = 2 vs short k = 12; tau on 32x32 cells. "
              "Prediction: tau(long) > 1.2 tau(short)",
        "honesty": "this is only an analogue: in Valentini's CMB argument the lag comes from modes that leave the "
                   "horizon (no dynamics at all), not from wavelength inside a box",
    },
    "test3_double_slit": {
        "3a_1D": "paraxial double slit in 1 transverse dimension (z = time), incident Gaussian (no nodes before "
                 "the slits). Wrong start = uniform across the open slits instead of |psi|^2. Prediction: NO "
                 "relaxation (1-D Bohm trajectories cannot cross, no vortices): fringe correlation of the wrong "
                 "start stays below the ideal-sampler band (mean - 3 sigma) at every screen distance",
        "3b_2D": "2 transverse dims (x,y), z = time. Upstream: square multimode guide [0,pi]^2 with M modes "
                 "(M = 16 vortex-rich, M = 1 none) for a length Z0, surfers started uniform (wrong) at z = 0; "
                 "at z = Z0 a double-slit mask in x, then free propagation (FFT) to the screen. Score: correlation "
                 "of the x-histogram of clicks with the exact x-marginal of |psi|^2 (bins of 8 fine bins). "
                 "Prediction: M = 16, Z0 = 8 pi -> correlation within 3 sigma of the ideal sampler; Z0 = 0 -> "
                 "below the band; M = 1 -> below the band for every Z0",
    },
    "test4_CHSH": "optional, skipped if over time budget",
}

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    try:
        git = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    except Exception:
        git = "unknown"
    PREREG["_meta"] = {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head": git}
    with open(os.path.join(OUT, "relax_prereg.json"), "w", encoding="utf-8") as f:
        json.dump(PREREG, f, indent=1, ensure_ascii=False)
    print("pre-registration written", PREREG["_meta"])
