"""aggmol_frat_prereg.py - PRE-REGISTRATION of study aggmol_frat_ ("test frattale con H2"): up the level loop of
docs/CICLO_LIVELLI.md with the H2-H2 force of aggmol phase 3: molecule -> molecule-ball -> cluster -> blob-ball -> drop.
Written after the development runs listed below and BEFORE the production runs. Writes lab/out/aggmol_frat_prereg.json."""
import json, os, subprocess, datetime
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = {
 "written": datetime.datetime.now().isoformat(timespec="seconds"),
 "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip(),
 "principle (Valerio)": "a rule is a rule only if it is the same at every zoom; we choose only how much each step holds. "
                        "Here the same form at every level: balls + a learned pair force + friction/noise obeying "
                        "fluctuation-dissipation (zero noise where nothing thermal is integrated out)",
 "declared up front": [
   "this is LIQUID H2 STANDING IN FOR WATER until oxygen passes the many-electron rule (docs/MOLECOLE.md §8.2)",
   "nuclei are CLASSICAL (no zero-point motion, no path integrals): they OVERBIND. Known size of the gap for H2: the "
   "zero-point kinetic energy in the liquid (~60-70 K per molecule) is about twice the pair well (~33 K); classical "
   "H2 with an accurate potential is expected to come out ~1.2-1.5x too dense and ~1.5-2.5x too cohesive, with a "
   "shorter first-neighbour distance and a larger surface tension. Thresholds below are 'with classical nuclei'",
   "classically a 1000-molecule cluster of this potential melts near T* = kT/eps ~ 0.5-0.6, i.e. 17-20 K: at 15 K it "
   "is expected SOLID (a crystallite), at 20 K liquid or close to melting (reported, not an exam)",
   "the fine level (rigid molecules, rank-2 anisotropic force from 5 orientations) inherits all debts of aggmol phase 3 "
   "(perturbative force, missing exchange-dispersion, rigid r = 1.4 a0, C6 -9.5 %) and truncates the anisotropy at rank 2",
 ],
 "development (declared)": [
   "fifth orientation S (parallel, 45 deg) computed with the phase-3 machine (aggmol_frat_aniso.py), C9 = 40.14 a.u. "
   "from the same pseudo-states (C6 by quadrature 10.861 = Casimir-Polder sum)",
   "MD engine checks: Mackay 13 / 55 / 147 minima -112.7 / -164.2 / -190.6 K per molecule (LJ with our eps: -114 / -170 / "
   "-199); NVE drift 0.011 K/molecule over 200 steps of 4.4 fs; 7-8 ms per step for 1000 molecules",
   "dev cube of 1000 at 20 K (gamma 1/t_u, 48 ps): stays one drop, asphericity 0.004-0.011, Epot/N ~ -156 K, a few "
   "molecules leave (vapour); at 15 K Epot/N ~ -176 K, nothing leaves",
   "dev rotors at R = 6.5 a0: ball V0 -33.414 K, quantum rotors -33.502 K (J<=4 and J<=6 identical), classical rotors at "
   "20 K free energy -35.60 K, classical minimum (T-shape) -52.13 K; equilateral trimer quantum -100.50 vs ball -100.24 K",
 ],
 "config": {
   "units": "MD in Angstrom, K, amu; t_u = 1.0967 ps; m(H2) = 2 m_p + 2 m_e = 2.01565 amu; B_rot = 1/(m_p r^2) = 87.74 K",
   "ball potential": "V0 = a0 + 2a1/3 + a2/3 + (a3+a4)/9 on the phase-3 grid (cubic spline 5-12 a0), exponential wall "
                     "below 5 a0, -C6/R^6 - C8/R^8 tail (phase-3 fit) above 12 a0, cut and shifted at 24 a0 (12.7 A)",
   "integrator": "BAOAB Langevin, dt = 0.004 t_u (4.4 fs), gamma = 0.2 / t_u (thermostat; dynamics reported, not examined)",
   "step 1": "dimers R = 5.5..10 a0: ball V0, quantum para rotors (even J <= 6), classical rotors (free energy over "
             "orientations at 10, 20, 50 K); equilateral trimers R = 6..7.5 a0: ball, quantum rotors (J <= 4), classical "
             "rotors at 20 K; Axilrod-Teller C9 term; B2(T) ball vs classical anisotropic at 20, 50, 100, 300 K",
   "step 2": "Mackay 13 / 55 / 147, local minimum, then MD at 5, 10, 15, 20 K (110 ps each, last 88 ps averaged): E/N, "
             "Lindemann index. Blob level: k-means blobs of ~13 molecules (k = round(N/13), centres carried from frame to "
             "frame), learned on the 147 cluster at 20 K (220 ps, frames every 0.44 ps): blob pair force by force matching "
             "(linear hats every 0.5 A from 2 to 30 A, ridge 1e-6), blob friction from the decay time of the blob-velocity "
             "autocorrelation, noise from fluctuation-dissipation; blob mass = members x m",
   "step 3": "1000 molecules from a 10x10x10 simple cube at the real liquid spacing (3.62 A), T = 20 K and 15 K, "
             "165 ps (checkpoints), averages over the last 110 ps: largest cluster (link 5.2 A), asphericity, core "
             "density (r < R_e/2), core g(r) first peak, core potential energy per molecule u_core, surface energy, "
             "Axilrod-Teller energy on 5 snapshots (perturbative, reported). Coarse level: the blob rule from step 2 on "
             "77 blobs from the same cube, same T, same time",
 },
 "exams": {
   "S1": "dimer, R = 6.0, 6.5, 7.0 a0: quantum para rotors minus ball, |dV| <= 2 % of |V0| (the ball is the right "
         "molecule for para-H2 at these T)",
   "S2": "equilateral trimer R = 6.5, 7.0 a0: quantum rotors minus ball <= 2 %",
   "S3": "(reported) classical rotors at 20 K minus ball, dimer and trimer: what a classical fine level would add",
   "S4": "(reported) Axilrod-Teller three-body / pair energy at the trimer minimum and in the drop",
   "C1": "(reported) cluster minima and melting by size",
   "D0": "drop forms by itself from the cube at 20 K: largest cluster >= 95 % of molecules and asphericity <= 0.02 "
         "averaged over the last 110 ps",
   "D1": "core density / 0.0708 g/cm3 in [1.10, 1.50] at 20 K (classical overbinding band)",
   "D2": "first g(r) peak of the core at 20 K in [3.15, 3.90] A (real ~3.7; classical expected shorter)",
   "D3": "Delta H_vap (classical) = -u_core - tail correction + T, ratio to 0.90 kJ/mol (108.3 K per molecule) in [1.4, 2.6]",
   "D4": "surface tension (rough) = (U_drop - N u_core) / (4 pi R_e^2), R_e = (3N / 4 pi rho_core)^(1/3); ratio to "
         "2.0 mN/m in [1.2, 4.0]",
   "L1": "LOOP EXAM: the blob rule learned on the 147 cluster predicts the 1000 drop at 20 K: radius of gyration of the "
         "blob centres within 10 % of the fine MD (k-means blobs of the fine drop, same k)",
   "L2": "LOOP EXAM: first peak of the blob-blob g(r) within 10 % of the fine MD",
   "L3": "(reported) the same at 15 K, and the blob rule applied back to the 147 cluster (self-consistency)",
 },
}
json.dump(P, open(os.path.join(OUT, "aggmol_frat_prereg.json"), "w"), indent=1)
print("written")
