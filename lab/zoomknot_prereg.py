"""zoomknot_prereg.py - pre-registration of the ZOOM-KNOT study (Z0-Z4), written BEFORE the production runs.
Writes lab/out/zoomknot_prereg.json. Engine: lab/zoomknot_core.py (committed together with this file)."""
import json, os, subprocess, datetime

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
head = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()


def lv(h, N, R):
    return dict(h=h, N=N, R=R)


HC = [lv(0.2, 24, None)]                       # coarse h = 0.2, 49^3, half box 4.8 (P4 used 48^3, half 4.7)
CONF = {
  "U0": HC,
  "Z1_h0.05_R0.6": [lv(0.05, 14, 0.6)] + HC,
  "Z1_h0.05_R1.0": [lv(0.05, 22, 1.0)] + HC,
  "Z1_h0.05_R1.4": [lv(0.05, 30, 1.4)] + HC,
  "Z1_h0.0667_R1.0": [lv(1 / 15, 17, 1.0)] + HC,          # optional (time permitting), h-dependence
  "Z1_h0.0667_R1.4": [lv(1 / 15, 23, 1.4)] + HC,          # optional
  "Z2": [lv(0.05, 14, 0.6), lv(0.1, 14, 1.2), lv(0.2, 14, 2.4), lv(0.4, 14, 4.8), lv(0.8, 14, 9.6),
         lv(1.6, 14, 19.2), lv(3.2, 14, 38.4), lv(6.4, 14, None)],
}
sites = {k: sum((2 * l["N"] + 1) ** 3 for l in v) for k, v in CONF.items()}

P = {
 "meta": {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": head,
          "hardware": "i5 4 cores shared with two other agents; 1 torch thread per run"},
 "question": "Valerio: 'where you zoom in you see the truth (fine balls), elsewhere precise agglomerates that introduce a "
             "small, controlled error'. Apply to the B=1 Skyrme knot whose core pressure on a uniform h=0.2 lattice is off "
             "by 57% (P4 of PRESSIONE.md). Goal: p(0) within 5% and E/12pi^2 within 0.3% of the exact ODE at a cost <= a "
             "64^3 uniform lattice (262144 sites).",
 "exact_reference": {"E_over_12pi2": 1.231446, "p0": 12.236, "Lambda": 106.8, "r0": 1.076,
                     "source": "knot_core.radial_bvp / press_ode.json (P1 of PRESSIONE.md)"},
 "engine": {
   "file": "lab/zoomknot_core.py",
   "grids": "nested cubic vertex grids, finest first, origin on a grid point; level l owns |x|_inf <= R_l; ghost layer of "
            ">= 2 points prolonged from the next coarser level by tensor cubic Lagrange interpolation + renormalisation; "
            "coarse points inside a finer box are injected from it (restriction). One field, no double counting.",
   "energy": "E = sum_l h_l^3 sum W_l [k2_l e2 + k4_l e4], e2/e4 = knot_core.static_density split (fb differences), "
             "W_l = trapezoid(box R_l) - trapezoid(box R_{l-1}) (conservative partition of volume).",
   "relaxation": "all levels jointly, L-BFGS (as knot_core.relax_lbfgs, 20 inner, history 30, strong Wolfe), variables "
                 "scaled by sqrt(h_l) per level; start = exact ODE hedgehog sampled on every level.",
   "stop_rule": "relative E change over 5 outer steps < 2e-8, OR 45 outer steps, OR 30 min total wall time per config "
                "(checkpointed every <= 9 min). Declared in the results.",
   "cost": "sites = all grid points of all levels including ghosts and covered points (upper bound of the work)."},
 "configs": CONF, "sites": sites, "budget_sites_64cubed": 64 ** 3,
 "learned_couplings_Z2": {
   "what": "for every level l >= 1 of Z2: two numbers k2_l, k4_l multiplying the bare e2, e4 lattice energy (16 numbers "
           "in total for 8 levels; finest level k = 1).",
   "how": "agglomerator principle (ZOOM.md: block, learn, check). Training configurations = exact hedgehog scaled "
          "f(r/lam), lam in {0.7, 0.85, 1.2, 1.4} (lam = 1 held out). Fine truth = bare lattice energy at spacing h_l/4 on "
          "the shell owned by level l (trapezoid weights). Blocked configuration = 2 steps of vertex full-weighting "
          "(1/4,1/2,1/4 per axis, the vertex-grid version of a 2x2x2 block), i.e. 4x4x4 fine balls -> 1 agglomerate, "
          "renormalised to the sphere. Least squares per term: k2 = sum E2f E2c / sum E2c^2 (same for k4).",
   "check": "held-out lam = 1: per-level energy of the exact profile with bare and learned couplings vs the continuum "
            "integral of the exact density over the same shell (trapezoid at h_l/8 of the analytic T00)."},
 "exams": {
   "Z0": "code check: U0 relaxed E/12pi^2 within 0.2% of 1.20578 (P4, half box 4.7 vs 4.8) and core max deviation of "
         "p (4th-order central, shells of 0.1, r < r0/2) > 30% of p(0). PASS expected.",
   "Z1a": "nested grids: at least one Z1 config with sites <= 262144 has |p0 - 12.236|/12.236 < 5% AND core max "
          "deviation (central) < 5%.",
   "Z1b": "error (core max deviation) decreases monotonically with the fine box R_f in {0.6, 1.0, 1.4} at h_f = 0.05.",
   "Z1c": "E/12pi^2 within 0.3% of 1.231446 for the best Z1. Expected FAIL (Dirichlet box 4.8 and bare h = 0.2 outside).",
   "Z2a": "Z2 (learned couplings, 195112 sites): |p0 - 12.236|/12.236 < 5% AND core max deviation < 5%.",
   "Z2b": "Z2: |E - E_exact|/E_exact < 2% (target) ; stretch goal < 0.3% (judged separately).",
   "Z2c": "Z2: |Lambda - 106.8|/106.8 < 3% (kinetic term bare, Dirichlet wall at 89.6).",
   "Z2d": "learning helps: (i) on held-out lam = 1, |learned per-level error| < |bare per-level error| for >= 6 of the 7 "
          "levels l >= 1; (ii) Z2-learned E error < Z2-bare E error (same grid, couplings = 1, relaxed separately).",
   "Z2e": "agglomerate error is controlled: summed static error of all levels >= 1 (learned, held-out lam = 1, signed "
          "sum) < 0.3% of E_exact.",
   "Z3": "error budget table (fine-core sites, agglomerate sites, p0 error, E error, Lambda error) for U0, Z1, Z2-bare, "
         "Z2-learned. 'Zoom principle' = U0 core error > 30%, best fine-core error <= 5%, agglomerates add < 2% in E and "
         "< 5% in p0 (pre-registered target X = 2% in E, 5% in p0). Judged from Z0/Z1a/Z2a/Z2b.",
   "Z4": "boosted v = 0.3 on the nested grid, E = gamma E0 within 3%: OPTIONAL, run only if budget remains; otherwise "
         "declared SKIPPED.",
   "reported_not_judged": "virial (E4-E2)/E, von Laue int p / int |p| (fb and central), B, r0, profile p(r)."},
 "pilot_disclosure": "Before writing this file ONE pilot run was made (Z1_h0.05_R1.0, 12 outer steps, 9 min, not converged "
                     "to the stop rule): E/12pi^2 = 1.21721, p0(central) = 11.35 (-7.3%), core max dev 6.8%, virial -1.3%. "
                     "The predictions below use it.",
 "predictions": {
   "Z0": "PASS", "Z1a": "FAIL (pilot -7% at R_f = 1.0; R_f = 1.4 exceeds the budget with h_f = 0.05)",
   "Z1b": "PASS", "Z1c": "FAIL (about -1.2%)",
   "Z2a": "uncertain, lean PASS (learned couplings should fix the global scale set by the coarse levels)",
   "Z2b": "2%: PASS; 0.3%: uncertain (50/50)", "Z2c": "PASS (missing tail ~ (16 pi/3) C^2/R = 0.8%)",
   "Z2d": "PASS", "Z2e": "PASS", "Z4": "SKIPPED (budget)"},
}
json.dump(P, open(os.path.join(OUT, "zoomknot_prereg.json"), "w"), indent=1)
print(json.dumps(sites, indent=1))
