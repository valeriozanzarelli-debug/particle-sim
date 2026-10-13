"""protetta_prereg.py - pre-registration of the study 'VELOCITA_PROTETTA' (prefix protetta_), written BEFORE the
production runs.  Question: is there a ball rule in which loops CANNOT separate the limiting speed of the electron from
that of light (Collins et al. PRL 93, 191301 (2004) tuning problem)?  Candidate tested: light = composite pair of the
same Weyl balls that make the electron (Bisio-D'Ariano-Perinotti, arXiv:1407.6928).
Run: python3 lab/protetta_prereg.py  -> lab/out/protetta_prereg.json"""
import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = {
 "study": "protetta_ (docs/VELOCITA_PROTETTA.md)", "date": "2026-10-13",
 "free_parameters": 0,
 "inputs": ["c, hbar, m_e (theta = m a), alpha", "frozen rule lab/zz_core.py (not modified)",
            "Debiti engine lab/vuoto_f3_core.py @35afcc5 copied read-only as lab/protetta_vf3.py (+ protetta_vcore, protetta_vf2)"],
 "sign_convention": "gap = c_gamma - c_e in units of alpha (c_e = dressed limiting speed of the electron, dv/alpha)",
 "declared_development_runs": [
   "dv_exact(theta=0.003125, R, K12/3, pd=0.05, axis x) = 0.048221 alpha (reproduces Debiti 35afcc5 exactly)",
   "dv_exact(theta=0.003125, R, K12/3, pd=10, axis x) = 0.047943 alpha (T -0.688, C +0.736: parts move, sum stable)"],
 "theory_predictions": {
   "T1_free_identity": "the rule's light walker (spin 1, S_k = -i eps_k, order xyz|zyx) is EXACTLY the triplet part of "
       "two R-Weyl walkers of the rule each at momentum k/2 (zero relative momentum): exp(-i k S) = exp(-i k/2 sigma) x "
       "exp(-i k/2 sigma) on the symmetric subspace, axis by axis, any order. Hence omega_gamma(k) = 2 omega_W(k/2) and "
       "c_gamma = c_W at every k (not only k->0). Same for the cyclic order C.",
   "T2_noncollinear": "a pair with relative momentum q perpendicular to k has energy k + 2 q_perp^2/k: the composite is a "
       "particle of mass 2|q_perp| and its speed is 1 - 2 (q_perp/k)^2 (no alpha, no a): kinematic cost of the composite.",
   "T3_massive_constituents": "if the constituents flip (theta > 0, i.e. they ARE electrons), the q = 0 pair at rest has "
       "quasi-energy 2 theta per tick: a 'photon' of mass 2 m_e. Excluded by design.",
   "T4_dressed_inheritance": "with charged flipless constituents, c_gamma(composite) = c_W(dressed). The dim-4 LV "
       "coefficient is soft-insensitive to the mass, so dv(p >> m) = dv(p << m) + O(theta^2, (p a)^2, alpha^2): "
       "D_j(p/m) = dv_j(p) - dv_j(0.05 m) -> 0 for a -> 0 at fixed p/m. Then the composite gap at the self-energy level "
       "is -D -> 0.",
   "T5_neutral_constituents": "neutral constituents (neutrino theory of light, de Broglie/Jordan) have no O(alpha) "
       "self-energy: c_gamma = 1 (bare), gap = -dv_e, i.e. the same as against the bare light walker (no gain).",
   "T6_not_covered": "the exchange (vertex) term between the two charged constituents is O(alpha) and is required by gauge "
       "invariance together with the self-energies; it is NOT computed here and no symmetry of the rule sets its LV part "
       "to zero. Prediction (qualitative): inheritance is exact only at the self-energy level."},
 "exams": {
   "L1": {"what": "max |eig(W_light(k)) - eig(triplet of W_R(k/2) x W_R(k/2))| (phases) over 2000 random k in the "
          "Brillouin zone, rules R and C; plus operator identity with a fixed intertwiner", "PASS_if": "< 1e-12"},
   "L2": {"what": "group velocity of light at k vs Weyl at k/2, axes x, y, z, (1,1,0), (1,1,1), |k| a = 1e-3 .. 0.5",
          "PASS_if": "< 1e-9"},
   "L3": {"what": "speed of a q_perp pair: (1 - v)/(2 (q_perp/k)^2) at k a = 0.05, q_perp/k = 0.01, 0.02, 0.05",
          "PASS_if": "within 5 % of 1"},
   "L4": {"what": "rest quasi-energy of the flipping q=0 pair / (2 theta) at theta = 0.05, 0.0125",
          "PASS_if": "within 1e-3 of 1 (prediction: the composite is massive -> not a photon)"},
   "D1": {"what": "D_j(pd) = dv_j(p = pd*theta) - dv_j(0.05 theta), rule R, cuts K12=(0.6,1.05) and K12/3=(0.2,0.35), "
          "axes x and y, pd = 1, 3, 10, theta in {0.0125, 0.00625, 0.003125, 0.0015625, 0.00078125} with p <= 0.04; "
          "extrapolated a -> 0 at fixed pd (fit d0 + d1 theta + d2 theta^2 if >= 3 points, else last point)",
          "PASS_if": "|d0| < 1e-3 alpha for every (cut, axis, pd) with >= 3 theta points",
          "meaning": "PASS = the charged composite inherits the electron's dressed c at the self-energy level"},
   "D2": {"what": "same for the cyclic order C, cut K12/3, axis x, pd = 10, theta = 0.003125, 0.0015625, 0.00078125",
          "PASS_if": "|D| < 1e-3 alpha (descriptive)"},
   "D3": {"what": "massless limit: theta = 3.90625e-4, p = 0.01, 0.02, 0.03 (pd = 25.6, 51.2, 76.8) vs dv(0.05 theta), R, both cuts, x, y",
          "PASS_if": "|D| < 2e-3 alpha"}},
 "requests_not_done_here": [
   "Debiti fase 4: c of the DRESSED elementary light walker of the rule (isotropic and along y) - not duplicated here "
   "(coordinator, 13/10). Until then every gap 'vs light' is against the BARE light walker.",
   "Exchange/vertex term inside the composite photon (two-loop-like): not done."],
 "verdict_rules": {
   "inherits": "D1 PASS -> 'the composite photon inherits c_e at the self-energy level (gap < 1e-3 alpha)', NOT 'solved': "
               "the experimental requirement is |gap| < ~3e-18 alpha and the vertex term is missing.",
   "fails": "D1 FAIL -> the inheritance is broken already at the self-energy level by mass effects."}
}
json.dump(P, open(os.path.join(OUT, "protetta_prereg.json"), "w"), indent=1, ensure_ascii=False)
print("written")
