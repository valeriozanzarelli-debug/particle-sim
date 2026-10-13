"""locale_f5_amend_prereg.py - DECLARED AMENDMENT to phase 5 (locale_f5_prereg.json, 5183ecf, is not modified), after the
correction from SURFISTI (docs/SURFISTI_ESATTI.md 12.8): the frozen rule has an exact time reversal with Theta^2 = -1.
Written BEFORE lab/locale_f5_kramers.py. Writes lab/out/locale_f5_amend_prereg.json."""
import json, os, subprocess, datetime
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
head = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
P = {
 "meta": {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": head},
 "reasoning_before_runs": "T = (i sigma_y on L and on R) K maps every shift to its inverse and exp(-i a X) to its inverse "
     "(X is real), so T W T^-1 = G W^-1 G^-1 with a local G: an exact time reversal (Kramers, T^2 = -1). The twisted "
     "Wilson uses exp(-i h Y) with Y = sigma_y (x) 1 on chirality: -iY is REAL, so T maps it to ITSELF, not to its "
     "inverse: the twisted Wilson breaks time reversal (as twisted mass does in lattice QCD). Time-reversal-even "
     "chirality-mixing generators: X (x) 1 (the Wilson along the mass, F2b of phase 1: kinetic mass +9 %) and "
     "Y (x) sigma_j. Candidate: the 'spin-link' twisted Wilson, link along e_j rotates with Y (x) sigma_j (cubic covariant, "
     "anticommutes with the mass X); the on-site part is undone exactly by the inverse of the link product at k = 0.",
 "exams": {
  "K1_control": "Kramers pairing of the period operator on 4^3 with a random potential (theta = 0.37): F0 max pair "
                "splitting <= 1e-10 (positive control, reproduces SURFISTI); F2 (Y Wilson, h = 0.1) splitting > 1e-6 "
                "(prediction: broken); F2b (X Wilson) <= 1e-10.",
  "K2_spinlink": "spin-link Wilson SL (h = 0.1, theta = 0.1): Kramers <= 1e-10; N-valley and zone-from-0 half-gap >= 0.3 / a; "
                 "|E(0)/m - 1| <= 1e-12; drift <= 1e-9; |m_kin / m_kin(frozen) - 1| <= 1e-4 (4 directions); unitarity <= 1e-12.",
  "K3_1s": "if K2 passes: reduced 1s model (48^3, Z = 60, theta = 0.5, corner, spin-up start, Ritz): the 1s doublet split "
           "<= 1e-5 (Y Wilson: 1.1e-3) and the main level has >= 90 % of its Fourier weight in the Gamma valley.",
  "doc": "par. 13.4 of docs/REGOLA_LOCALE.md is corrected (old sentence struck through with a note)."}}
with open(os.path.join(OUT, "locale_f5_amend_prereg.json"), "w") as f:
    json.dump(P, f, indent=1, ensure_ascii=False)
print("written")
