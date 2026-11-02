"""emiss6_g0.py - exam G0 (phase 4): the light entered through its response at the atom (kernel) gives the SAME run as
the full 128^3 light box of the production code (lab/emiss_run.py), for the production coupling (centered + transverse).
Z = 60, radial rule a = 0.1, pure 2p3/2, K = 1, T = 1500.  Measure max |norm_kernel - norm_box| / (1 - norm_box(T)).
Run in pieces: python lab/emiss6_g0.py  (exit 3 = again)."""
import os, sys, json
import numpy as np
os.environ["EMISS_CFG"] = "g0_box"
import emiss_configs as EC
from emiss_run import Runner
from emiss6_core import qe_run, levels, OUT

Z, A, RMAX, T = 60, 0.1, 40.0, 1500.0
lev = levels(Z, A, RMAX)
lam = 2 * np.pi / lev["w"]
cfg = EC.radial_cfg("g0_box", Z, A, RMAX, T / lam, cells_per_lambda=30)
r = Runner(cfg)
if not r.run(float(sys.argv[1]) if len(sys.argv) > 1 else 500):
    sys.exit(3)
nb = np.array(r.ser["norm"]); tb = np.array(r.ser["t"])
res, meta = qe_run("g0_kernel", Z, A, RMAX, None, T_abs=tb[-1], parity=False)
if res is None:
    sys.exit(3)
nk = res["norm"]
m = min(len(nk), len(nb))
assert np.allclose(res["t"][:m], tb[:m]), (res["t"][:5], tb[:5])
d = float(np.max(np.abs(nk[:m] - nb[:m])) / (1 - nb[m - 1]))
out = dict(exam="G0", T=float(tb[m - 1]), loss_box=float(1 - nb[m - 1]), loss_kernel=float(1 - nk[m - 1]),
           max_rel_diff=d, threshold=1e-3, verdict="PASS" if d < 1e-3 else "FAIL", b=meta["b"], sub=meta["sub"])
print(out)
json.dump(out, open(os.path.join(OUT, "emiss6_g0.json"), "w"), indent=1)
