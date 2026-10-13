"""arrivo_r2_rule1d.py KIND L D a - production runs of the rule in 1+1 D (configs of lab/out/arrivo_prereg.json).
KIND: 'z'   = spin along the guide (alpha = 0): flux at L + hybrid surfers (samples NR_B for L=10)
      'x'   = spin along x (alpha = pi/2): hybrid surfers (+ flux), samples NR_B (L=10) / NR_A first 5000 (L=100)
      'det' = the rule's cell detector at z = L (no trajectories)
      'zramp' = as 'z' with the absorbing ramp instead of the light cone (edge check, D = 20)
Resumable: each call works <= 540 s; prints DONE when finished. Results: scratch r2/<tag>.npz"""
import os
import sys
import time
import numpy as np
import arrivo_rule1d as R1

SCR = os.environ.get("ARRIVO_SCRATCH", "/tmp/claude-0/-home-claude-particle-sim/9de435bc-26c2-5d88-b70d-07efb8033e5b/scratchpad")
os.makedirs(os.path.join(SCR, "r2"), exist_ok=True)

kind, Lz, D, a = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
tag = "%s_L%g_D%d_a%g" % (kind, Lz, D, a)
outp = os.path.join(SCR, "r2", tag + ".npz")
if os.path.exists(outp):
    print("DONE (exists)", tag)
    sys.exit(0)
ck = os.path.join(SCR, "r2", tag + "_ckpt.npz")
t0 = time.time()
if kind == "det":
    tau_end = 100.0 if D < 80 else 60.0
    if os.path.exists(ck):
        c = np.load(ck)
        run = R1.Run1D(D, a, Lz * D, "detector", tau_end * D**2)
        run.R[:], run.L[:] = c["R"], c["L"]
        run.tick = int(c["tick"])
        run.absorbed = [c["ab"]]
    else:
        run = R1.Run1D(D, a, Lz * D, "detector", tau_end * D**2)
    n_end = int(np.ceil(tau_end * D**2 / a))
    while run.tick < n_end and time.time() - t0 < 540:
        run.advance(min(50000, n_end - run.tick))
    ab = np.concatenate(run.absorbed)
    if run.tick < n_end:
        np.savez(ck, R=run.R, L=run.L, tick=run.tick, ab=ab)
        print("PARTIAL", tag, run.tick, n_end, "%.0f s" % (time.time() - t0))
        sys.exit(0)
    F = np.cumsum(ab)
    tau = (np.arange(ab.size) + 1) * a / D**2
    np.savez(outp, tau=tau[::max(1, ab.size // 20000)], F=F[::max(1, ab.size // 20000)], norm_left=run.norm(),
             F_end=F[-1], omega_box=run.omega_box)
    if os.path.exists(ck):
        os.remove(ck)
    print("DONE", tag, "F_end %.5f norm_left %.5f" % (F[-1], run.norm()))
    sys.exit(0)

# trajectories
if Lz == 10.0:
    smp = np.load(os.path.join(SCR, "arrivo_r1_NR_B.npz"))
    nS = 4000
    omega = 1000.0
else:
    smp = np.load(os.path.join(SCR, "arrivo_r1_NR_A.npz"))
    nS = 5000
    omega = 1000.0
XYZ0 = (smp["X0"][:nS], smp["Y0"][:nS], smp["Z0"][:nS])
alpha = np.pi / 2 if kind == "x" else 0.0
if Lz == 10.0:
    tau_end = 6.0 if kind == "x" else (40.0 if D < 80 else 10.0)
else:
    tau_end = 50.0
edge = "cone" if (D == 20 and kind != "zramp") else "ramp"
ramp = (3 * Lz * D, Lz * D, 0.05)
res = R1.run_hybrid(D, a, Lz, omega, alpha, XYZ0, tau_end, edge=edge, ramp=ramp, do_traj=True,
                    ckpt=ck, budget_s=500, verbose=False, chunk=20000)
if res is None:
    print("PARTIAL", tag, "%.0f s" % (time.time() - t0))
    sys.exit(0)
step = max(1, res["tau"].size // 40000)
# keep the flux on a coarse grid (bin averages) for the record
nb = res["tau"].size // step
fl = res["flux"][:nb * step].reshape(nb, step).mean(axis=1)
tt = res["tau"][:nb * step].reshape(nb, step).mean(axis=1)
np.savez(outp, tau=tt, flux=fl, T=res["T"], norm=res["norm"], omega_box=res["omega_box"], edge=edge,
         flux_int=np.sum(res["flux"]) * a / D**2)
print("DONE", tag, "arrived %.4f mean %.4f max %.4f (%.0f s)" % (np.isfinite(res["T"]).mean(), np.nanmean(res["T"]),
                                                               np.nanmax(res["T"]), time.time() - t0))
