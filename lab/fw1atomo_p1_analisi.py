"""fw1atomo_p1_analisi.py - reads the P1 results in the scratch folder and prints, per configuration, the exact
levels (residual < 1e-6, weight >= 1 %) in [0,84; 0,92]: omega, weight, valleys, <r>, w(r<1), spin, sublattice.
Run: python lab/fw1atomo_p1_analisi.py [prefix]"""
import os, sys, json, glob
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import fw1atomo_core as fc

pre = sys.argv[1] if len(sys.argv) > 1 else ""
for fn in sorted(glob.glob(os.path.join(fc.SCR, "p1", pre + "*.json"))):
    d = json.load(open(fn))
    print("==", os.path.basename(fn)[:-5], {k: d["config"][k] for k in ("rule", "h", "n", "a", "pos")}, f"{d['s_per_period']:.3f} s/per")
    for l in d["levels"]:
        if l["resid"] < 1e-5 and max(l["w_start"]) >= 0.01 and 0.80 <= l["omega"] <= 0.93:
            cl = int(np.argmax(l["classes"]))
            print(f"  {l['omega']:.5f} res {l['resid']:.0e} w {max(l['w_start']):.3f} G {l['val_Gamma']:.3f} N {l['val_N']:.3f} pi {l['val_pi']:.3f}"
                  f" r {l['r_mean']:.2f} in1 {l['w_r_lt_1']:.3f} S {np.round(l['spin'], 2)} cls {cl}:{max(l['classes']):.2f}")
