"""nucleo1s_f5_amend.py -> lab/out/nucleo1s_f5_amend1.json : DECLARED AMENDMENT to nucleo1s_f5_prereg.json (commit 89e99f2),
written after the first bare scans and BEFORE any learning at Z = 20. No threshold changed.
Why: at Z = 20 the pre-registered bare listening (T = 500) shows ONE line only (the 4-valley multiplet is not resolved:
at Z = 20 the binding is 0.0107 mc^2 and the splits are much smaller than at Z = 60), so the nucleus-blind level U is not
defined. The bare scan at Z = 20 is redone with T = 2000 (the listening length of the learning). The T = 500 record is
kept ('bare_Z20_0.5'). Also seen (reported in the doc): bare U at Z = 40 is +9.2 % from Dirac, at Z = 80 +0.6 %; the
"U within 1 %" of Z = 60 is not universal.
Junction: the case k = 0.2 with sigma = 10 has a momentum spread of 25 %, which mixes in nearly stopped components
(the 1-D prototype already showed R ~ 2 % at k = 0.1). Added cases with sigma = 25 at normal incidence (k = 0.2, 0.4),
REPORTED next to the pre-registered ones; J2 is judged on the pre-registered cases as written."""
import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
json.dump(dict(amends="lab/out/nucleo1s_f5_prereg.json (commit 89e99f2)", bare_T={"20": 2000.0},
               junction_extra=[dict(a=0.25, k=[0.2, 0, 0], Ny=1, Nz=1, sig=25.0),
                               dict(a=0.25, k=[0.4, 0, 0], Ny=1, Nz=1, sig=25.0)]),
          open(os.path.join(OUT, "nucleo1s_f5_amend1.json"), "w"), indent=1)
print("written")
