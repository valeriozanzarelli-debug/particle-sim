"""knot2_prereg.py - pre-registration of the knot2 study (S1-S4), written BEFORE any production run.
Dev seen before writing (declared): K1-K5 results of knot_*.py (knot_results.json, NODO_3D.md); the analytic vacuum
dispersion of the 2nd-order and 4th-order staggered stencils (math, no simulation); a timing test of one RK4 step
(old 1.1 s, c4 1.8 s on 112x48x48); a check that M^-1 (M v) = v to 1e-14 for old/c2/c4; the radial profiles of the
OLD relaxed lattice knots (h = 0.2 L48/L64, h = 0.15 L64: f(r) is 10-15 % below the ODE at r ~ 2, i.e. the lattice
knot is SMALLER than the continuum one) and the old Lambda(R) numbers (known since K4: 63.3/72.3/77.1, 58.0/69.5/76.5).
"""
import json, os, subprocess, datetime
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=HERE).stdout.strip()
P = {
 "schemes": {
  "old": "K1-K5 engine: static energy with forward/backward 2nd-order differences; kinetic matrix uses the CENTRAL gbar=(g+ + g-)/2 in the -(v.gbar)^2 Skyrme piece (mixed stencils)",
  "c2": "same static energy, kinetic matrix built like the static one: -(1/2) sum_s (v.g^s)^2 (consistent, 2nd order)",
  "c4": "4th-order STAGGERED forward/backward differences g+ = (-phi[x+2]+27phi[x+1]-27phi[x]+phi[x-1])/24h, g-(x)=g+(x-1), used identically in static and kinetic terms (consistent, 4th order, no checkerboard zero mode)"},
 "S1_dispersion": {
  "S1a": "small-amplitude (eps = 1e-3, float64) standing plane pion waves on a PERIODIC lattice, h = 0.2, N = 64 along the wave, k h = 2 pi m / 64, m = 1..32, along an axis (box 64x2x2) and along the (1,1,0) diagonal (box 64x64x2); omega measured by fitting A(t) = A0 cos(omega t) to the projected amplitude over >= 1.5 periods of RK4 (dt = 0.05 h); phase velocity omega/|k|, group velocity d omega/d|k| by finite differences of the MEASURED omega. Code check (PASS/FAIL): measured omega equals the analytic lattice dispersion within 0.5 % for every m. Reported: max |v_phase - 1| and max |v_group - 1| for k h <= pi/4, <= pi/2 and the full band, with the k where they occur, for 2nd order (old = c2 in vacuum) and c4.",
  "S1b": "mismatch on a knot background: for perturbations with x-gradient only on the relaxed OLD rest knot (h = 0.2, L = 48), local light speed along x c_eff^2 = min over tangent delta of delta.K.delta / delta.M.delta with K the static x-stiffness built from the transverse (y,z) gradients and M the kinetic matrix without x-gradient terms; c_eff = 1 identically for c2 and c4 (same stencil). Reported: min c_eff for old and where.",
  "S1c": "Cherenkov thresholds (analytic): smallest k h where omega(k)/k < v along an axis; old: v=0.3 none, v=0.6 none, v=0.8 kh=2.26; c4: none, none, kh=2.89. Reported: fraction of the contracted knot's x-spectrum (|FT pi_a|^2 along k_x) above that k h, h = 0.2 and 0.15.",
  "S1_prediction": "if the K2 slowdown is driven by the stencil MISMATCH, c2 at h = 0.2 fixes v = 0.6 as well as c4; if driven by the vacuum DISPERSION of short waves, c2 does not and c4 does (ablation runs below decide; stated before the runs)."},
 "S2_boost": {
  "runs": "scheme c4. Relaxed c4 knots by L-BFGS (start = old relaxed field) at h = 0.2 L = 48 and h = 0.15 L = 64. Boost exactly as K2 (knot_k2.py): phi(gamma(x-x0)), v_field = -v d_x phi, boxes 112x48x48 (h = 0.2) and 150x64x64 (h = 0.15), dt = 0.25 h, float32; h = 0.2: v = 0 (T = 6), 0.3 (T = 8), 0.6 (T = 6), 0.8 (T = 6); h = 0.15: v = 0.6, 0.8 (T = 6), and v = 0.3 (T = 8) if time allows. Rest-ratio reference at h = 0.15 = the h = 0.2 c4 rest run (no h = 0.15 rest run, cost). Ablation (informative, decides S1_prediction): c2 at h = 0.2, v = 0.6 and 0.8.",
  "exams (each v, each h)": {
   "E": "|E/(gamma E0) - 1| < 3 %",
   "contraction": "baryon sigma_x/sigma_y averaged over the 2nd half, divided by the rest ratio, = 1/gamma within 3 %",
   "speed": "baryon centroid velocity (linear fit, 2nd half) within 3 % of v",
   "shape": "energy-density sigma_x(end)/sigma_x(0) within 5 %",
   "valid": "energy drift < 0.5 %"},
  "min_balls": "report the smallest contracted r_half/h (r_half of the relaxed profile / gamma / h) for which all four pass"},
 "S3_Lambda": {
  "definition": "Lambda = 2 T for a rigid isorotation (w = 1) of the relaxed knot, T = v.M.v with the scheme's kinetic matrix",
  "tail_correction": "linear tail theory: outside the core each pion component is harmonic, pi_a = C (D_a - H_a), D_a = x_a/r^3, H_a harmonic in the cube with H_a = D_a on the walls (Dirichlet vacuum). Missing inertia dLambda = (4/3) C^2 Q / R_eff, Q = R [ int_out |D|^2 + int_cube (2 D.H - |H|^2) ] (pure number for a cube, computed numerically once; a sphere gives 4 pi 9/5), R_eff = (L+1) h / 2 (first vacuum ghost layer). C fitted on the lattice field itself: least squares of pi_a = C (D_a - H_a) over 2.5 <= r <= 3.5. Lambda_corr = Lambda_lat + dLambda. Remainder expected O(1/R^3).",
  "primary": "c4, h = 0.15, L = 64: Lambda_corr within 2 % of the ODE value 106.77 (PASS/FAIL). Consistency (reported): c4 h = 0.2 L = 48 and L = 64 corrected values.",
  "secondary": "(i) the same tail correction applied to the OLD fields (h^2 Richardson of h = 0.2 L48 and h = 0.15 L64, as K4); (ii) fit Lambda(R) = Lam_inf - a/R_eff - b/R_eff^3 to c4 h = 0.2 L = 32, 48, 64 (3 points, exact) and compare a with (4/3) C^2 Q.",
  "mu_ratio": "|mu_p/mu_n| with ANW formulas (one input M_Delta/M_N fixes e), from the primary lattice Lambda_corr and the lattice E and <r^2>_B of the same field; compare with experiment 1.460 (PASS if within 3 %) and with the ODE 1.417."},
 "S4_annihilation": {
  "runs": "scheme c4, h = 0.2, box 80^3 (half 7.9), sponge in the outer 1.4 (strength 2), knot at x = -2 and anti-knot (pi -> -pi) at +2, product ansatz, taper (3.5, 4.3), relaxed c4 profile (h = 0.2), v = 0.6 and 0.3 each head-on, T = 24 (30 if time allows), dt = 0.25 h.",
  "exams": {
   "radiation": "E_in(r < 4.5) drops by >= 70 % of E_in(0) by t = 24 (PASS/FAIL); report the time it crosses 70 %",
   "unknotting": "integral |b| < 0.2 (as K3a)",
   "massless": "pion fields pi_a and d pi_a/dt in the shell 4.5 <= r <= 6.3 (smooth cos^2 mask), 3D FFT, per |k| bin omega_est(k) = sqrt(sum|FT dpi/dt|^2 / sum|FT pi|^2), averaged over t = 16, 20, 24. Fit omega_est^2 = c^2 k^2 + m^2 over the bins holding 90 % of the spectral energy: PASS if |c - 1| < 10 % and m^2 < (0.25)^2 (i.e. |m| < 0.25 Skyrme units, about 1/4 of the dominant k) and the energy-weighted mean of omega_est / omega_lattice(k) is within 10 % of 1. Report the dominant wavelength 2 pi / <k>."}},
 "S5": "docs/NODO_3D.md: new section appended (Italian, numbers, failures, debts) and figure lab/out/knot2_figure.png",
 "budget": "each script <= 10 min wall (resumable chunks of 330 s), <= 2 threads, ~3 h total; runs dropped for time are declared, not hidden"}
res = dict(meta=dict(written=datetime.datetime.now().isoformat(timespec="seconds"), git_head_before=head), prereg=P)
json.dump(res, open(os.path.join(OUT, "knot2_prereg.json"), "w"), indent=1)
print(json.dumps(res, indent=1)[:400])
