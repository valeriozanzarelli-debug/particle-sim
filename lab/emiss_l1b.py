"""emiss_l1b.py - test L1b of lab/out/emiss_prereg2.json: the old beat in the regenerated Laboratory export.
Run: python lab/emiss_l1b.py -> lab/out/emiss_l1b.json"""
import json, os
import numpy as np
from emiss_core import OUT

D = "/mnt/project-files/ripresa/dati/atomo_accoppiato"


def amp(t, x, lo=0.010, hi=0.016):
    m = t > 400
    t, x = t[m], x[m] - x[m].mean()
    h = np.sin(np.pi * (t - t[0]) / (t[-1] - t[0]))**2
    om = np.arange(lo, hi, 1e-5)
    S = 2 * np.abs(np.exp(1j * np.outer(om, t)) @ (h * x)) / np.sum(h)
    k = int(np.argmax(S))
    return float(om[k]), float(S[k])


d = np.genfromtxt(os.path.join(D, "serie.csv"), delimiter=",", names=True)
man = json.load(open(os.path.join(D, "manifest.json")))
w_new, a_new = amp(d["t"], d["p_exc"])
lb = json.load(open(os.path.join(OUT, "emiss_lab1s.json")))
fr = lb["free"]
t_old = np.array(fr["t"])
p_old = np.array(fr["p_run"][:len(t_old)])
w_old, a_old = amp(t_old, p_old)
res = dict(export_omega_1s=man.get("omega_1s"), new=dict(omega=w_new, amplitude=a_new), old_free_run=dict(omega=w_old, amplitude=a_old),
           p_exc_range_new=[float(d["p_exc"][d["t"] > 400].min()), float(d["p_exc"][d["t"] > 400].max())], pass_=bool(a_new < 0.01))
json.dump(res, open(os.path.join(OUT, "emiss_l1b.json"), "w"), indent=1)
print(f"L1b: new export (omega_1s {man.get('omega_1s')}): largest line in [0.010, 0.016] at {w_new:.5f}, amplitude {a_new:.4f} "
      f"(threshold 0.01) -> {'PASS' if a_new < 0.01 else 'FAIL'}; old 1s free run: {w_old:.5f}, amplitude {a_old:.4f}; "
      f"p_exc range t > 400 now {res['p_exc_range_new']}")
