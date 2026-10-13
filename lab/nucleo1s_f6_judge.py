"""nucleo1s_f6_judge.py - exams A1-A3 of phase 6 -> lab/out/nucleo1s_f6_judged.json"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA, MC2_EV  # noqa: E402
from nucleo1s_core import dirac_binding  # noqa: E402
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "nucleo1s_f6_prereg.json")))
R = json.load(open(os.path.join(OUT, "nucleo1s_f6_results.json")))
th = np.array(P["config"]["theta_out"])
J = {"shifts": {}}
for Z in (1, 2):
    Za = Z * ALPHA
    for n in (1, 2, 3):
        y = np.array([R[f"Z{Z}_th{t}_n{n}"]["B"] for t in th])
        ext = np.polyfit(th**2, y, 2)[-1]
        ext2 = np.polyfit(th[1:]**2, y[1:], 1)[-1]
        S = Za**2 / (2 * n * n); D = dirac_binding(n, -1, Za)
        J["shifts"][f"Z{Z}_{n}s"] = dict(dirac_shift_eV=(D - S) * MC2_EV, rule_shift_eV=(ext - S) * MC2_EV,
                                         rel=(ext - S) / (D - S) - 1, rel_2pt=(ext2 - S) / (D - S) - 1,
                                         level_rel_err=ext / D - 1)
sh = J["shifts"]
J["A1"] = dict(PASS=bool(all(abs(sh[f"Z1_{n}s"]["rel"]) <= 0.03 for n in (1, 2))))
J["A2"] = dict(PASS=bool(all(abs(sh[f"Z2_{n}s"]["rel"]) <= 0.03 for n in (1, 2))))
J["A3"] = dict(PASS=bool(all(abs(sh[f"Z{Z}_{n}s"]["rel"] - sh[f"Z{Z}_{n}s"]["rel_2pt"]) <= 0.02
                             for Z in (1, 2) for n in (1, 2))))
J["cells"] = {k: v["cells"] for k, v in R.items() if k.endswith("_n1")}
json.dump(J, open(os.path.join(OUT, "nucleo1s_f6_judged.json"), "w"), indent=1)
for k, v in sh.items():
    print(k, f"Dirac {v['dirac_shift_eV']:.4e} eV  rule {v['rule_shift_eV']:.4e}  rel {100*v['rel']:+.2f} %  2pt {100*v['rel_2pt']:+.2f} %  level {v['level_rel_err']:+.2e}")
print({k: J[k]["PASS"] for k in ("A1", "A2", "A3")})
