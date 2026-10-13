"""aggfull_prereg.py - pre-registration of AGGLOMERATORE COMPLETO, written BEFORE the production runs.
Writes lab/out/aggfull_prereg.json. Engine: lab/aggfull_core.py, aggfull_learn.py, aggfull_proton.py, aggfull_atom.py,
aggfull_run.py, aggfull_fineref.py (committed together with this file)."""
import json, os, subprocess, datetime
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
head = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
P = {
 "meta": {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": head,
          "hardware": "cloud container, 4 cores shared with other agents; 1 thread, <= 2.5 GB"},
 "question": "Make the coarse rule learned by the agglomerator free of hand-written pieces: (1) the coarse potential, "
             "(2) spin-orbit/Darwin at the coarse level, (3) the 2s at 1 % without a huge lattice, (4) the proton recoil "
             "from the proton moving as a heavy ball of the same rule (debts #8, #9, #12 of CONTESTO_CLOUD.md section 4).",
 "engine": {
   "coupling": "fine 3-D zig-zag walk (theta = a0/4096, chain A) in a weak potential plane wave; first-order response over "
               "B ticks seen through the positive band in the Foldy-Wouthuysen-aligned basis (unitary polar factor of the "
               "old large-component descriptor); Duhamel inversion -> coupling of the coarse GENERATOR, "
               "<k'|V_c|k> = Vt(k'-k) Gamma(k',k). Gamma fitted by relative least squares with a Hermitian sandwich stencil "
               "sum_jl C_jl e^{-ik'j} e^{ikl}; candidates point (1 site), star (7), cube (27), chosen automatically: smallest "
               "passing the check. Pairs (k,k'): 75 % drawn like the observable state (size a0 up to cell a0/16, 2 a0 for a0/8) "
               "+ 25 % uniform in the cube (weight 0.03). Repeat: level l+1 learned from the LEARNED level l; a level with no "
               "admissible stencil is flagged and the next level learns from the level below it.",
   "coupling_check": {"eps_abs": 1e-4, "eps_rel_kdep": 2e-2, "tame_max": 1.5},
   "field": "coarse field = band-limited restriction (coarse Fourier modes, no aliasing: 8 half-cell sub-samples per coarse "
            "cell) of the fine input -Z alpha/r (nucleus on a grid point; translation-covariant, no corner/centre choice).",
   "atom": "H = T_e(k) [+ T_p(-k)] + gamma0 V_c (eigenproblem, LOBPCG, parity sector kept), then dGamma = Gamma - gamma0 "
           "(Darwin-like, spin-orbit, rest) in first order; the kinetic remainder <T_e - k^2/2m_kin> (m_kin from the learned "
           "generator at k = 1/64..1/16 per cell) is the learned relativistic kinetic term.",
   "proton": "fine zig-zag walk with theta_p = m_p a_p, a_p = (a0/16)/2^19 (theta_p = 0.029995); same auto_level of "
             "agg_zz_core with b in {16,8,4,2}; observable: proton at v = alpha c (k = alpha m_p, the same 2.2e-4 per fine "
             "cell as chain A), capped at 1/16 per cell; final cells a0/16 and a0/8 (= the electron's). Generators "
             "symmetrised k <-> -k (fine walk parity checked exact to 2e-16).",
   "recoil": "pair state in the zero-total-momentum sector (exact for translation-invariant rules): T_e(k) + T_p(-k); "
             "the reduced mass is NOT an input.",
 },
 "configurations": {
   "a0/16": {"level_B": 256, "N": 256, "box_a0": "+-8"},
   "a0/8": {"level_B": 512, "N": 256, "box_a0": "+-16", "box_checks": "N = 128 (+-8 a0) for 1s, N = 192 (+-12 a0, the old box) for 1s and 2s"},
   "moving_atom": {"level_B": 512, "N": 128, "P": "n x 4 pi/L, n = 0,1,2 (P/2 a grid vector), split 1/2 (no masses used)"},
   "old_restriction_diagnostic": "cell average of -alpha/r (agg_zz_atom), with the learned kinetic, N = 192 and 256",
 },
 "exams": {
   "L1": "learned coupling vs the fine walk DIRECTLY at the same B (held-out state pairs), levels B = 256 and 512: mean "
         "relative error of the k-dependent part <= 2 % and mean absolute error <= 1e-4",
   "L2": "Darwin read off the learned stencil: gk / (-1/8) and gx / (+1/4) within 2 % (B = 256 and 512)",
   "L3": "spin-orbit read off the learned stencil: gso / (1/4) within 2 % (B = 256 and 512)",
   "E1": "1s1/2 at a0/16: E0 + first-order dGamma vs Dirac (infinite mass) within 0.1 %; also reported vs the fine level "
         "(zig-zag radial walk at the same theta, aggfull_fineref.py)",
   "E2": "1s1/2, 2s1/2, 2p1/2, 2p3/2 at a0/8 (E0 + first order) vs Dirac within 0.1 % each",
   "S1": "2s1/2 at a0/8 vs Dirac within 1 % (the old A2b exam, redone)",
   "F1": "fine structure 2p3/2 - 2p1/2 (6x6 first-order matrix of the learned dGamma on the coarse 2p) vs Dirac "
         "(infinite mass, 10.94966 GHz) within 1 %",
   "F2": "pattern: lower pair and upper quartet each degenerate within 1 % of the splitting; Lande interval ratio "
         "(upper - mean)/(mean - lower) = 0.5 within 1 %",
   "F3": "alpha^4 shifts (Dirac - Schroedinger) of 1s, 2s, 2p1/2, 2p3/2 at a0/8 (and 1s at a0/16) = <kinetic remainder> "
         "+ <dGamma> (first order) within 5 % each",
   "R1": "proton learned kinetic mass (parity-even generator at k = 1/64..1/16 per cell) at cells a0/16 and a0/8: "
         "|m_kin/m_p - 1| <= 1e-3",
   "R2": "recoil shift (E_pair/E_static - 1) / (-m_e/(m_e + m_p)) within [0.99, 1.01] for 1s (a0/16), 1s, 2s, 2p (a0/8)",
   "R3": "moving atom: M from E(P) - E(0) = P^2/2M; M/(m_e + m_p) within 1e-3",
 },
 "not_testable_declared": "beyond-reduced-mass recoil (Barker-Glover / Breit, ~ alpha^4 m_e^2/m_p ~ 1e-9 relative) is below "
                          "the resolution and needs the transverse photon exchange that the rule (static Coulomb phase) "
                          "does not have; QED (Lamb) absent as in ZIGZAG.md.",
 "development_runs_before_this_file": [
   "aggfull_dev_pot.py and scratchpad scripts: band-limited potential WITH aliasing (pseudo-spectral) gave the 1s at a0/16 "
   "+0.08 % (too bound) -> replaced by the alias-free Galerkin field; exact k^2/2 kinetic + Galerkin: 1s at a0/8 "
   "-1.24e-4 (box +-8..12 a0), cell average -0.46 %; 2s at a0/8 in the old box +-12 a0: +1.17 % (box effect suspected: "
   "periodic images of the proton); box effect of the 1s at +-6 a0 ~ +5e-4",
   "dev_aggfull_learn: first version with uniform pairs in plain least squares failed (level 1 relativistic cube); "
   "relative least squares adopted; B = 16 has no admissible stencil (tameness 1.83 for the cube stencil) -> flagged; "
   "B = 256 star: gso 0.2503, gk -0.1252, gx 0.2502; B = 512 star: gso 0.2505",
   "dev_aggfull_proton: k_obs = 1/a0 (1e-7 per proton fine cell) unresolvable; own-cell packets too demanding at "
   "relativistic levels; adopted v = alpha c observable capped at 1/16 per cell: hierarchy complete, m_kin even part "
   "1.0002; learned odd part (fake drift ~1e-10 c) found and removed by parity symmetrisation",
 ],
}
json.dump(P, open(os.path.join(OUT, "aggfull_prereg.json"), "w"), indent=1)
print("written", os.path.join(OUT, "aggfull_prereg.json"))
