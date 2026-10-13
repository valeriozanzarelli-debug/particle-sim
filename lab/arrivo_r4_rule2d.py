"""arrivo_r4_rule2d.py KIND SPIN a - production runs of the rule in 2+1 D (parallel-plate guide, axes y, z).
KIND 'free': surfers (Bohm-Dirac trajectories of the rule, NO spin term added), samples NR_C (first 2000);
KIND 'det' : the rule's cell detector at z = L. SPIN in x, y, z. D = 30, omega = 20, L = 10 (DD units).
Resumable (<= 540 s per call). Result: scratch r4/<tag>.npz"""
import os
import sys
import time
import numpy as np
import arrivo_rule2d as G

SCR = os.environ.get("ARRIVO_SCRATCH", "/tmp/claude-0/-home-claude-particle-sim/9de435bc-26c2-5d88-b70d-07efb8033e5b/scratchpad")
os.makedirs(os.path.join(SCR, "r4"), exist_ok=True)
kind, spin, a = sys.argv[1], sys.argv[2], float(sys.argv[3])
SPIN = dict(x=(1, 1), y=(1, 1j), z=(1, 0))
tag = "%s_%s_a%g" % (kind, spin, a)
outp = os.path.join(SCR, "r4", tag + ".npz")
if os.path.exists(outp):
    print("DONE (exists)", tag)
    sys.exit(0)
ck = os.path.join(SCR, "r4", tag + "_ckpt.npz")
D, om, Lz = 30, 20.0, 10.0
if kind == "free":
    smp = np.load(os.path.join(SCR, "arrivo_r1_NR_C.npz"))
    YZ0 = np.array([smp["Y0"][:2000], smp["Z0"][:2000]])
    cfg = dict(D=D, a=a, omega=om, Lz=Lz, ymax=0.8, mode="free", spinor=SPIN[spin], tau_end=8.0, zc=2 * Lz,
               ell=0.6 * Lz, YZ0=YZ0)
else:
    tau_end = 60.0 if a > 0.2 else 10.0
    cfg = dict(D=D, a=a, omega=om, Lz=Lz, ymax=0.8, mode="detector", spinor=SPIN[spin], tau_end=tau_end, YZ0=None)
done, res = G.run_2d(cfg, ck, budget_s=520, verbose=False)
if not done:
    print("PARTIAL", tag)
    sys.exit(0)
save = dict(tau=res["tau"], norm=res["norm"], prep_resid=res["prep"]["resid"], prep_Om=res["prep"]["Om"],
            spin0=res["prep"]["spin"])
if kind == "free":
    save.update(T=res["T"], flux=res["flux"])
else:
    save.update(F=np.cumsum(res["absorbed"]))
np.savez(outp, **save)
if os.path.exists(ck):
    os.remove(ck)
if kind == "free":
    T = res["T"]
    print("DONE", tag, "arrived %.4f mean %.4f max %.4f" % (np.isfinite(T).mean(), np.nanmean(T), np.nanmax(T)))
else:
    print("DONE", tag, "F_end %.5f" % np.sum(res["absorbed"]))
