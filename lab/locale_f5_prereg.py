"""locale_f5_prereg.py - pre-registration of PHASE 5 of REGOLA LOCALE, written BEFORE lab/locale_f5_run.py and
lab/locale_f5_split.py. Earlier pre-registrations are not modified. Writes lab/out/locale_f5_prereg.json."""
import json, os, subprocess, datetime
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
head = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
P = {
 "meta": {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": head},
 "question": "Breathe the light too, so that electron and light have the same limiting speed (c cos mu), while the "
             "pi-valleys of the electron get a mass. Plus: scaling with the lattice step of the 1s spin-doublet split "
             "from the broken time reversal (phase 4, X1).",
 "light_rule": "the photon as the pair of Riemann-Silberstein fields F+ = E + icB (shift with S.n, as the frozen light walk) "
               "and F- = E - icB (shift with -S.n), i.e. 6 components, like L and R of the electron; between the moves of "
               "each tick a local coin exp(-i s_t mu X6), X6 swaps F+ <-> F-, s = (-1, +1) on the two ticks (average zero: "
               "no photon mass). Note found while designing (stated before the runs): a continuous local mixing of F with "
               "F* cannot keep E and B real (the unitary matrices [[a, b], [b*, a*]] need a b = 0), so this is the light "
               "as a one-photon wave function (complex, positive frequency), not a classical real field. Without coin it "
               "is two independent copies of the frozen light walk.",
 "electron_rule": "lab/locale_f4_core.walk(theta, mu_x = mu) (breathing coin of VUOTO), optionally + twisted Wilson hu = 0.1.",
 "configs": "mu in {0.3, 0.5}; theta in {0.1, 0.02}; with and without hu = 0.1.",
 "exams": {
  "E1_speed": "limiting speed of the electron c_e (from E^2 = E0^2 + c_e^2 p^2, fit p in [0.02, 0.1] m, 4 directions) vs "
              "light c_l (phase velocity at k = 1e-3, 4 directions): |c_e/c_l - 1| <= 1e-4 at theta = 0.02 (the frozen "
              "pair has its own lattice mismatch ~ theta^2/6: reported; at theta = 0.1 the exam is on the CHANGE with "
              "respect to the frozen mismatch, <= 1e-4).",
  "E2_pi_mass": "electron: half-gap of the 4 pi-valleys and the zone minimum from pi >= 0.3 / a (theta = 0.1); with hu = 0.1 "
                "also the N valleys and zone minimum from 0 >= 0.3 / a (then ALL doublers).",
  "E3_electron": "|E(0)/theta - 1| <= 1e-12; kinetic mass consistent with the common speed: |m_kin c_l^2 / E0 - 1| minus "
                 "the same quantity for the frozen pair, <= 1e-4 (theta = 0.1) - i.e. the breathing adds nothing beyond "
                 "the common rescaling of c.",
  "E4_maxwell": "light N_1% (worst of phase, group, diffraction; errors measured against the measured c_l) in 3-D <= 1.1 x "
                "frozen (36.3 -> 40); no photon mass |omega(k = 0)| <= 1e-12.",
  "E5_charge": "unitarity of electron and light period operators <= 1e-12; real-space norm conservation over 50 periods on "
               "24^3 <= 1e-12 (both).",
  "E6_info": "helicity purity of the light bands at k = 2 pi / 20 (weight of F+ in the + band) and its dependence on mu.",
  "verdict": "the breathing pair (electron + light) is ADOPTABLE AS A PROPOSAL if E1-E5 pass for some mu with E2.",
  "T1_split": "spin-doublet split of the 1s (reduced model of phase 3/4: Z = 60, nucleus on the corner, twisted Wilson "
              "h = 0.1, initial spin-up on one copy, Ritz on the 2-tick operator) for theta in {0.7, 0.5, 0.35, 0.25} "
              "(box 48^3 for theta >= 0.35, 64^3 for 0.25; periods ~ 350/theta). Doublet = the two largest-weight levels in "
              "the window whose spin expectations are opposite (cosine <= -0.9). Report Delta(theta) and the power law; "
              "information only (no threshold), but 'vanishes with the step' if the fitted exponent is >= 1."},
}
with open(os.path.join(OUT, "locale_f5_prereg.json"), "w") as f:
    json.dump(P, f, indent=1, ensure_ascii=False)
print("written")
