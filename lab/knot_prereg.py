"""knot_prereg.py - PRE-REGISTRATION of K1-K5 ('nodo 3D': the Skyrme model on a lattice of balls).
Writes lab/out/knot_prereg.json. Committed BEFORE the production runs (knot_k1.py ... knot_k5.py).

DEV (declared, looked at, NOT production numbers):
 - radial ODE (shooting, massless): E/12pi^2 = 1.23148, E2 = 72.82, E4 = 73.03, B = 1.000000, Lambda = 103.90 Skyrme units
   (= (2pi/3) * 49.6 in ANW units; ANW quote 50.9), r(f = pi/2) = 0.875, f ~ 1.99 / r^2 tail.
 - ball of radius R with f(R) = 0 (1D minimisation): E/12pi^2 = 1.2525 (R=4), 1.2382 (R=6), 1.2345 (R=8): finite-size
   excess ~ 1/R^3, not ln R.
 - 4th-order CENTRAL differences: the hedgehog unwound completely (E -> 0, B -> 0) under arrested Newton flow at h = 0.3
   (zero-energy checkerboard mode, as already met in 2D in LUCE_INTRAPPOLATA). Switched to a forward/backward symmetric
   scheme (sigma term: mean of forward and backward |d phi|^2; Skyrme term: |g_i ^ g_j|^2 averaged over the 4 forward /
   backward combinations): no zero modes, 2nd order. With this scheme the hedgehog STILL unwinds at h = 0.3 and h = 0.4:
   lattice Derrick scan E(mu) shows no barrier against shrinking at h >= 0.3 (the Skyrme term is under-counted when the
   core spans < 3 balls). At h = 0.2 a barrier exists; relaxed E/12pi^2 = 1.2058 at h = 0.2, box half-size 4.8
   (-2.1 % from the continuum, i.e. O(h^2) error) and B = 0.9956.
 - dynamics (RK4 + projection, h = 0.2, dt = 0.1 h): energy conserved to 1e-5 over t = 1; boosted ODE-profile at v = 0.6
   gives E / (gamma E0) = 0.996 but the centroid moved at ~0.556 instead of 0.6 in the first time unit: the continuum
   ODE profile is not the lattice minimum (forces at t = 0). Production therefore boosts the RELAXED LATTICE profile.
 Therefore the production spacings are h = 0.2 and h = 0.15 (r_half/h = 4.4 and 5.8 balls); h >= 0.3 is reported as
 the 'too few balls' failure, not used for exams.
"""
import json, os, subprocess, datetime, math

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

PRE = {
    "model": "Skyrme model, massless pions (m = 0), Skyrme units: energy F_pi/(4e), length 2/(e F_pi) (ANW convention "
             "F_pi = 2 f_pi). Static density |d phi|^2 + 1/2[(tr D)^2 - tr D^2]; kinetic |v|^2(1+|d phi|^2) - sum (v.d_i phi)^2. "
             "Lattice: forward/backward symmetric differences, Dirichlet phi = (1,0,0,0) on the walls. Bogomolny bound 12 pi^2 |B|.",
    "parameter_count": {
        "Lagrangian": "2 (F_pi, e); pion mass not used (would be a 3rd)",
        "shape_and_ratio_exams_in_Skyrme_units": 0,
        "numerical_choices": "h, box, dt, taper, sponge (convergence parameters, not physics)",
        "K4_physical_numbers": "only dimensionless e matters for dimensionless ratios; it is fixed from ONE measured ratio "
                               "M_Delta/M_N (a 1-input fit). F_pi then from M_N (2nd input). Everything else declared below as "
                               "prediction, fit, or identity.",
    },
    "K1": {
        "runs": "B = 1 hedgehog, arrested Newton flow, dt = 0.15 h, stop when relative energy change over 100 iterations < 1e-9 "
                "or 2500 iterations; L^3 = 32^3, 48^3, 64^3 at h = 0.2 (half-box 3.2, 4.8, 6.4) and h = 0.15 (2.4, 3.6, 4.8); "
                "start = ODE profile; extra: crude start f = pi exp(-r/0.6) at h = 0.2, L = 48; informative: h = 0.25, L = 38.",
        "literature": "E/12pi^2 = 1.2315 (massless; ANW M = 36.5 F_pi/e = 1.232 x 3 pi^2) - verified by our ODE (1.23148)",
        "exams": {
            "K1a_continuum": "h^2-Richardson from (h=0.2, L=48) and (h=0.15, L=64) (same half-box 4.8): |E/12pi^2 - 1.2315| < 1 %",
            "K1b_box": "h = 0.2: |E(L=64) - E(L=48)| / E(L=64) < 0.5 %",
            "K1c_power_not_log": "h = 0.2: [E(32) - E(48)] / [E(48) - E(64)] >= 2.5 (1/L^3 gives ~4, ln L gives 1.4)",
            "K1d_topology": "|B - 1| < 0.01 for all production lattices and the crude start reaches the same E within 0.2 %",
        },
    },
    "K2": {
        "runs": "relaxed lattice profile (h = 0.2) boosted along x: phi(x) = phi_s(gamma (x - x0)), v_field = -v d_x phi; "
                "box 112 x 48 x 48, dt = 0.1 h, v = 0, 0.3, 0.6, 0.8, run time such that the knot travels >= 3 radii (t = 6, "
                "v = 0.3: t = 8); v = 0.8 repeated at h = 0.15 (box 150 x 64 x 64) if time allows; reference E0 = lattice "
                "energy of the unboosted same profile in the same box.",
        "exams (each v = 0.3, 0.6, 0.8 at the finest h run)": {
            "K2a_energy": "|E / (gamma E0) - 1| < 3 %",
            "K2b_contraction": "baryon density: (sigma_x / sigma_y) averaged over the second half of the run = 1/gamma within 3 % "
                               "(rest run v = 0 sets the reference ratio, which must be 1 within 1 %)",
            "K2c_shape_kept": "energy density sigma_x at end / at t = 0 within 3 % (contracted shape is a travelling solution)",
            "K2d_speed": "centroid velocity (baryon density, linear fit over second half) within 3 % of v",
            "K2e_min_balls": "report the smallest contracted r_half / h for which K2a-d pass",
            "energy_drift": "< 0.5 % (else the run is declared numerically invalid)",
        },
        "informative": "energy-density sigma_x/sigma_y vs 1/gamma (T00 is not a Lorentz scalar: not expected to be exactly 1/gamma)",
    },
    "K3": {
        "runs": "product ansatz U = U1 U2, h = 0.2, box 64^3 (half 6.4), sponge in the outer 1.4 (sigma <= 2), centres +-2 on x, "
                "dt = 0.1 h. B = 1 + B = -1 (antiSkyrmion = hedgehog with pi -> -pi): v = 0 (each at rest), 0.3, 0.6 head-on; "
                "t_max = 12. B = 1 + B = 1 in the attractive channel (2nd isorotated by pi about z), v = 0.5 each, t_max as fits.",
        "exams": {
            "K3a_unknotting": "S + antiS: integral |b| d^3x falls below 0.2 (from ~2) before t_max, for every speed",
            "K3b_radiation": "S + antiS: [E_in(0) - E_in(t_max)] / E_in(0) >= 70 % with E_in = energy inside the sphere r < 4.5",
            "K3c_B2": "S + S: B stays 2 within 2 %, and a 90-degree scattering (lumps leave along an axis perpendicular to x) "
                      "or a bound toroidal B = 2 state (reported, literature: torus, 90-degree scattering)",
            "control": "a single Skyrmion in the same lattice and time keeps B within 1 % (else annihilation could be a lattice artefact)",
        },
    },
    "K4": {
        "inputs": "Lambda from the relaxed lattice hedgehog (h = 0.2 and 0.15) by rigid isorotation, kinetic energy = Lambda w^2 / 2",
        "exams": {
            "K4a_inertia": "lattice Lambda (Richardson in h^2) within 2 % of the radial ODE value 103.90",
            "K4b_mu_ratio": "|mu_p/mu_n| (e fixed from M_Delta/M_N = 1232/938.92: a 1-input fit) within 3 % of experiment 1.460; "
                            "ANW obtained 1.43",
            "K4c_expected_failures": "g_A within 30 % of 1.2754 and F_pi within 30 % of 184 MeV: literature (ANW 0.61, 129 MeV) says both FAIL",
            "spin_half": "NOT computed: Finkelstein-Rubinstein (pi_1 of configuration space = Z2) is cited, not derived",
        },
        "classification": "0-input: E/12pi^2, Lambda (Skyrme units), r_E/r_B, rotor ratio (M_{5/2}-M_{3/2})/(M_{3/2}-M_{1/2}) = 5/3, "
                          "mu_NDelta = (mu_p-mu_n)/sqrt2 (identity of the collective quantization). 1-input (M_Delta/M_N): "
                          "|mu_p/mu_n|, g_A, F_pi/M_N, M_N r. 2-input: F_pi, radii in fm, absolute moments.",
    },
    "K5": {
        "render": "isosurface (marching tetrahedra) of the baryon density at 10 % of its maximum, relaxed B = 1, h = 0.2, L = 64",
        "exams": {
            "K5a_cost": "drawn triangles / number of balls < 1 %",
            "K5b_fidelity": "isosurface centroid within h of the true centre and principal-radius anisotropy < 5 %",
        },
    },
}

if __name__ == "__main__":
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=HERE).stdout.strip()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "knot_prereg.json"), "w", encoding="utf-8") as f:
        json.dump({"meta": {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": head},
                   "prereg": PRE}, f, indent=1, ensure_ascii=False)
    print("written lab/out/knot_prereg.json")
