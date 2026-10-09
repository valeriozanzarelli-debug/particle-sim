"""kk_shape_prereg.py - PRE-REGISTRATION of the 'forma nascosta' tests (written and run BEFORE any test).

Writes lab/out/kk_shape_prereg.json with a timestamp. Every later script reads its pass/fail criteria from
here; nothing below is changed after the tests are run (if a criterion turns out badly posed, the doc says
so instead of editing it).
"""
import json, datetime, os

PREREG = {
    "written": None,
    "shapes": {
        "control": "disk",
        "ellipse_aspects": [1.5, 2.0, 3.0],
        "triangle": "equilateral",
        "stadium": "Bunimovich, straight half-length a = cap radius r (a/r = 1)",
        "equal_area": "all 2-D shapes have the same area (pi R^2 sites); ratios do not depend on the scale",
        "boundary": "mass wall: flip angle theta = pi/2 outside the shape (infinitely heavy outside), theta_in = 0 inside",
    },
    "T1_spectra": {
        "claim": "the rest masses of the ball are the eigenmodes of the walk on the hidden shape",
        "check_disk": "walk levels (one doubler copy, M=+1 block) agree with the continuum infinite-mass Dirac "
                      "disk J_l(kR) = +-J_{l+1}(kR) (Berry-Mondragon 1987), using R_eff = sqrt(sites/pi) (no fit); "
                      "PASS if relative error of the lowest 6 levels < 5% at R = 20 and DECREASES from R = 10 to 20 to 40",
        "check_laplacian": "5-point Dirichlet Laplacian: disk vs j_{l,m}^2/R^2, equilateral triangle vs "
                           "(16 pi^2/9 s^2)(m^2+mn+n^2); PASS if lowest 6 within 3% at R = 30",
        "check_product": "circle x disk, charge n = 1, 2: E_n^2 = (k0_n)^2_lattice + E_shape^2 within 2% where "
                         "k0 enters through the lattice dispersion of the circle axis alone (sin-type) -- reported both ways",
        "check_convergence_other_shapes": "ellipse/triangle/stadium: lowest 6 levels x R_eff change < 3% from R=20 to 40",
    },
    "T1_charge": {
        "loop_only": "hidden space = ONE closed loop of N = 96 sites; charge n = momentum around it. Shape enters "
                     "only through the curvature kappa(s) of the embedding: (a) flip theta(s) = 0.1 kappa(s)/<kappa> "
                     "('mass = curvature', Valerio), (b) da Costa potential V = -kappa^2/(8 m) (imposed, debt). "
                     "Prediction: circle drift < 1e-12; ellipse 2:1, trefoil eps=0.2, stadium drift > 0.05 over 400 ticks.",
        "plane_only": "hidden space = 2-D shape only; candidate charge = J_z = L_z + sigma_x/2 about the centre. "
                      "Prediction: disk drift small (lattice square symmetry only, < 0.1), ellipse/triangle/stadium "
                      "drift > 0.3 over 400 ticks, starting from <J_z> ~ 1.5.",
        "product": "hidden space = circle (8 sites) x shape: charge n on the circle drift < 1e-12 for every shape.",
    },
    "T2_statistics": {
        "observable": "mean spacing ratio <r~> (Atas et al. 2013): Poisson 0.386, GOE 0.531, GUE 0.600",
        "laplacian": "desymmetrized quarter shapes, Dirichlet, levels 50..(k h < 0.6): quarter disk and golden "
                     "rectangle -> <r~> < 0.45 (Poisson-like); quarter stadium -> <r~> in [0.49, 0.57] (GOE)",
        "walk": "quarter stadium mass-wall walk (one doubler copy, M=+1 block), levels E in (0.05, 0.8): "
                "PASS if <r~> > 0.5 (Wigner-Dyson; GOE or GUE not predicted in advance)",
        "scar": "full stadium Dirichlet Laplacian, lowest 400 modes: exists a mode whose weight in a tube of "
                "half-width pi/(2k) around the UNSTABLE horizontal axis orbit, measured in the caps |x| > a, "
                "is >= 3x the tube area fraction",
    },
    "T3_leptons": {
        "targets": "m_mu/m_e = 206.768, m_tau/m_e = 3477.2, Koide Q = 0.66666",
        "pass_ratio": "within 1% on both ratios",
        "E0_shape_modes": "0 parameters per fixed shape (ellipse aspect scan = 1 parameter): lowest three distinct "
                          "rest masses of the walk on the shape. Prediction (Ashbaugh-Benguria: lambda2/lambda1 <= 2.539 "
                          "in 2-D for the Laplacian): m2/m1 < 2 for every shape -> FAIL",
        "E1_curvature_mass": "1 parameter mu: hidden loop (in the product circle x loop, charge sector fixed) with "
                             "flip theta(s) = mu kappa(s)/<kappa>. Fixed shapes chosen now: ellipse 2:1, trefoil eps=0.2, "
                             "stadium a=r, egg (0.15, 0.10). Fit mu to m2/m1 = 206.768 if reachable, then predict m3/m1 "
                             "and Q. Expected: m2/m1 not reachable -> FAIL",
        "E1pm_signed": "same with theta(s) = mu (kappa/<kappa> - 1) (signed mass, domain walls). Expected: index "
                       "theorem gives exactly massless modes (integral of theta = 0) -> report",
        "E2_localization": "exponential hierarchy from domain walls on the walk: three negative-mass segments of "
                           "lengths d, 2d, 3d in a +M loop (1 parameter M d, shape ratio 1:2:3 fixed now). Fit M d to "
                           "m2/m1, predict m3/m1 (~ (m2/m1)^2) and Q. Free segment lengths (2 parameters) = fit, no prediction.",
        "parameter_rule": "a ratio is a PREDICTION only if it was not used to fix any parameter",
    },
}

if __name__ == "__main__":
    PREREG["written"] = datetime.datetime.now().isoformat(timespec="seconds")
    os.makedirs(os.path.join(os.path.dirname(__file__), "out"), exist_ok=True)
    p = os.path.join(os.path.dirname(__file__), "out", "kk_shape_prereg.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(PREREG, f, indent=1, ensure_ascii=False)
    print("pre-registered:", p, PREREG["written"])
