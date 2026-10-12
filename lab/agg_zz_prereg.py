"""agg_zz_prereg.py - pre-registration of AGGLOMERATO ZIG-ZAG (A0-A3), written BEFORE the production runs.
Writes lab/out/agg_zz_prereg.json. Engine: lab/agg_zz_core.py (committed together with this file)."""
import json, os, subprocess, datetime

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)
head = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
A0 = 137.035999          # Bohr radius in lambda_C (= 1/alpha)

P = {
 "meta": {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": head,
          "hardware": "i5 4 cores shared with other agents; <= 2 threads"},
 "question": "Valerio: 'real hydrogen in the zig-zag needs ~1e10 cells; where did the agglomerator go?' Answer by doing it: "
             "learn the coarse zig-zag rule from the fine rule (block, learn, step, check, repeat), run hydrogen at the Bohr "
             "scale with few cells, keep the fine zig-zag only where we zoom.",
 "units": "hbar = c = m_e = 1; length lambda_C, time hbar/(m c^2). theta = m a (zz_core, frozen).",
 "engine": {
   "file": "lab/agg_zz_core.py",
   "block": "b cells per axis -> 1 (mean x b^{d/2}); b steps of the level below -> 1 step (tau = B fine ticks: c kept)",
   "learn_1D": "kernel psi'(X) = sum_{|j|<=w} K_j psi(X+j) (2x2), least squares (normal equations by FFT, ridge 1e-11) on "
               "(blocked state, blocked state one coarse step later) of random fields evolved by the level below (the fine "
               "level: the zig-zag walk run tick by tick); data spectrum Gaussian of width ks = clip(0.3 m B, 0.05, 0.6) per "
               "coarse cell inside the coarse Brillouin cube; the kernel's unitary polar factor is used ('the balls are "
               "conserved').",
   "learn_3D": "plane-wave least squares (same information as evolved packets, by linearity and translation invariance, "
               "without a box): the level below applied to e^{ik.x} (fine level: the walk's two-tick Bloch matrix, checked "
               "equal to the running walk to 6e-16), b steps, seen through the new descriptor. Pauli descriptor: the learned "
               "object is the GENERATOR G(k) = sum_j G_j e^{ik.j} (Hermitian stencil, |j|_inf <= w), one step = e^{i phi0} "
               "e^{-iG}; relative least squares (lstsq) on k drawn like the observable state (75%) + uniform over the "
               "Brillouin cube (25%, weight 0.03); motion phase unwrapped into [-0.5, 2pi - 0.5). 4-component descriptor: "
               "kernel fit (learn_k) with polar factor.",
   "descriptors": "1-D: (R, L) block means, both bands. 3-D: 'pauli2' (large component phi = (L+R)/sqrt2 of positive-energy "
                  "states; tried first) or 'dirac4' (both bands).",
   "check": "in k-space (exact for linear translation-invariant rules = running the packets): (i) packets after "
            "ceil(1000/B) coarse steps (= 1000 fine ticks), |dpsi|^2/|psi|^2 <= 1e-3; 1-D: sigma = 4 cells, |mean k| = 0.5 per "
            "cell; 3-D: physical packets of the observable (mean k = k_s/2, sigma = 2/(k_s B) cells). (ii) energy: 1-D: "
            "max |d eps| over k <= 0.3 m relative to eps(0.3 m) <= 1e-4; 3-D: <|d eps|>/<eps> with k drawn from the "
            "momentum density of a 1s-like state of size 1/k_s (k_s = 1/a0 up to the cell a0/16, then 1/(2 a0) for the n = 2 "
            "level a0/8) <= eps_E = 3e-4 per level; rest-phase error in the same units <= eps_E; (iii) 3-D tameness: the "
            "generator has no negative motion energy anywhere in the Brillouin cube (<= 1% of <eps>) and the 90th "
            "percentile of the relative motion-phase error over the cube <= 0.2.",
   "auto_rule": "AUTOMATIC level selection: candidates b in {2,3,4} (1-D) or {16,8,4,3,2} (3-D) (b x the lower level's step "
                "must be a multiple of its period: 2 ticks for the fine 3-D walk), w in {1..6} (1-D) / {1..4} (3-D), the "
                "descriptor (3-D); choose: fewest components, then largest b, then smallest w, among candidates passing every "
                "check. Stop at the observable's cell (zoom <-> level: cell <= a0/16 for 1s, <= a0/8 for n = 2). No admissible "
                "candidate -> the hierarchy stops there (flag).",
   "repeat": "level l+1 is learned from the LEARNED level l, never from the fine rule.",
 },
 "chains": {
   "A": {"theta_fine": A0 / 16 / 256, "note": "fine cell a = a0/4096 = 0.03346 lambda_C (kinetic-mass lattice error "
         "theta^2/3 = 3.7e-4), so that a block of 256 = a0/16"},
   "A1_ladder": {"theta_fine": A0 / 16 / 256, "b": [2, 2, 2, 2, 2], "w": "smallest passing, same check",
                 "note": "the fixed b = 2 ladder 2,4,8,16,32 requested in the task, for the A1 table"},
   "B_zoom": {"theta_fine": A0 / 16 / 32, "note": "coarser fine rule (a = 0.2676 lambda_C, resolves the Zitterbewegung "
              "period pi/m in 11.7 ticks) for the hybrid zoom run A3 only; its own lattice error (theta^2/3 = 2.4% on the kinetic "
              "mass) is shared by hybrid and uniform-coarse runs and does not enter the A3 exam"},
 },
 "exams": {
   "A1a": "1-D, every level of the b=2 ladder and of the automatic hierarchy: coarse rule vs the TRUE fine walk on Gaussian "
          "packets (sigma = 4 coarse cells, |mean k| < 0.5 per coarse cell, random spinor), |dpsi|^2/|psi|^2 < 1% after 1000 "
          "fine ticks (rounded down to whole coarse steps). Cost reduction (cells x ticks) >= b^2 = B^2 reported, also with "
          "the stencil size (taps x C^2 per cell per step) included.",
   "A1b": "kinetic mass from the learned positive-branch dispersion: fit eps(k) - eps(0) = k^2/(2 m_kin) - c4 k^4 + c6 k^6 "
          "over |k| <= min(0.3 m, 0.8/B): |m_kin/m - 1| < 1% at every level with B >= 32 (where the window holds >= 0.27 "
          "per cell). The rest mass is reported as the rest phase m tau mod 2 pi vs theta B mod 2 pi (it becomes a pure "
          "phase at the coarse level).",
   "A1c": "c4 vs 1/(8 m^3): within 10% at every level with B >= 32 (the fine rule's own c4 is reported as well).",
   "A1d": "automatic hierarchy (1-D and 3-D): built without hand intervention; levels, b, w, descriptor, check error and "
          "cost exported to assets/atomo/agglomerato_zz.json. Informative (no threshold) besides A1a-c on its levels.",
   "A2a": "3-D, chain A, deepest auto level with cell <= a0/16: coarse rule W = e^{-iV tau/2} e^{-iG} e^{-iV tau/2} "
          "(V = cell average of -alpha/r, nucleus at a cell corner, rest phase divided out), periodic box 192^3 (+-6 a0); "
          "eigenvector by preconditioned inverse iteration (LOBPCG) on H = G/tau + V, energy read from the rule: "
          "E = -arg<x|W|x>/tau (residual reported). |E_1s/(-13.605693 eV) - 1| < 0.1% (infinite nuclear mass).",
   "A2b": "2s and 2p (deepest auto level with cell <= a0/8, box 192^3 = +-12 a0): each within 1% of -3.401423 eV and "
          "|E_2s - E_2p| < 1% of (E_2p - E_1s).",
   "A2c": "BONUS fine structure 2p1/2-2p3/2 vs 10.94 GHz within 20% -- only if the learned kernel carries spin. Declared in "
          "advance: a pauli2 kernel learned from free data cannot carry spin-orbit (it is a property of the coupling to the "
          "potential through the small component, not of free propagation); if pauli2 is chosen the bonus is reported as "
          "NOT AVAILABLE together with the measured spin-dependent part of the kernel.",
   "A3a": "hybrid zoom (chain B): learned coarse rule everywhere, the fine zig-zag (4 components, point Coulomb phase per tick) "
          "in the coarse cell holding the proton (nucleus at its centre; cube of side ~8.6 lambda_C, inscribed sphere 4.3 "
          "lambda_C), halo of one coarse cell prolonged from the coarse field each coarse step (spectral interpolation of "
          "the block means, deconvolved, small component from the positive-energy relation), restriction = block mean of "
          "phi. 1s energy from the phase slope of <psi_0|psi(t)> on the coarse field vs the uniform-coarse run on the same "
          "box: |E_hyb/E_coarse - 1| < 0.1%.",
   "A3b": "Zitterbewegung only in the fine region: the time series of the fine-region observable (dipole <z> of the fine "
          "patch interior) has a spectral peak within 10% of 2 E_1s ~ 2 m c^2/hbar; the coarse descriptor has no negative "
          "band (by construction) and its time step tau = 8.6 hbar/mc^2 > period pi hbar/mc^2: the trembling cannot exist there.",
   "A3c": "cost of hybrid vs uniform fine (cells x ticks) on the same box and time: reported.",
 },
 "development_runs_declared": [
   "1-D: b = 2 is exact (two ticks of the walk are a trigonometric polynomial in 2k: the blocked map closes); b = 4 with "
   "w = 2..4 reproduces 1000 fine ticks at 1e-5..1e-9 with band-limited training.",
   "3-D first try with sharp band-limited (sphere) training: out-of-band kernel |K| up to 150 (unstable) -> training spectrum "
   "changed to the full coarse cube with Gaussian weight (engine above). pauli2 at b = 2 from the fine walk: training error "
   "2e-3 (w = 1), 1.3e-4 (w = 2): the positive-band projection of an ultra-relativistic cell (B << lambda_C) is non-local.",
 ],
}
json.dump(P, open(os.path.join(OUT, "agg_zz_prereg.json"), "w"), indent=1)
print("written", os.path.join(OUT, "agg_zz_prereg.json"))
