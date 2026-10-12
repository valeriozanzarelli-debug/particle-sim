"""press_lattice_posthoc.py - POST-HOC (not pre-registered) diagnosis of the P4 FAIL.
usage: python press_lattice_posthoc.py h L
For a relaxed lattice knot: (1) core slope a = -df/dr at the centre vs the ODE (2.0075);
(2) p(r) measured with the SAME forward/backward discretisation the lattice energy uses (kc.static_density):
    p_fb = (e4_fb - e2_fb)/3, vs the pre-registered 4th-order-central measurement;
(3) core deviation of both from the ODE, to see whether it shrinks like h^2.
Writes lab/out/press_posthoc_h<h>_L<L>.json."""
import sys, os, json, math, time, numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot_core as kc
torch.set_num_threads(2)
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
h = float(sys.argv[1]); L = int(sys.argv[2])
t0 = time.time()
X, Y, Z = kc.grid((L, L, L), h)
phi = kc.hedgehog(X, Y, Z, kc.radial_bvp(60.0))
hist = None
while True:
    phi, hist, done = kc.relax_lbfgs(phi, h, maxit=200, tol=1e-11, wall=60, hist=hist)
    print(f"iters {len(hist)} E/12pi2 {hist[-1]/(12*math.pi**2):.8f} done {done} t {time.time()-t0:.0f}s", flush=True)
    if done or time.time() - t0 > 380:
        break
with torch.no_grad():
    rr = torch.sqrt(X**2 + Y**2 + Z**2)
    # central
    g = kc.derivs(phi, h)
    Dm = torch.stack([torch.stack([(g[i] * g[j]).sum(0) for j in range(3)]) for i in range(3)])
    trD = Dm[0, 0] + Dm[1, 1] + Dm[2, 2]
    trD2 = torch.einsum("ik...,ki...->...", Dm, Dm)
    pc = (0.5 * (trD**2 - trD2) - trD) / 3
    # forward/backward, same as the lattice energy
    gp, gm = kc.derivs_fb(phi, h)
    e2 = 0.5 * sum((gp[i] * gp[i]).sum(0) + (gm[i] * gm[i]).sum(0) for i in range(3))
    etot = kc.static_density(phi, h)
    e4 = etot - e2
    pf = (e4 - e2) / 3
    E2, E4 = float(e2.sum() * h**3), float(e4.sum() * h**3)
    f = torch.arccos(phi[0].clamp(-1, 1))
R = rr.numpy().ravel(); fr = f.numpy().ravel()
m = (R > 0) & (R < 0.45)
a_lat = float(np.polyfit(R[m], np.pi - fr[m], 1)[0]) if m.sum() > 3 else float("nan")
O = np.load(os.path.join(OUT, "press_profiles.npz"))
p0, r0 = float(O["p"][0]), float(O["r0"])
edges = np.arange(0, 4.7, h / 2); idx = np.digitize(R, edges) - 1; nb = len(edges) - 1
sel = idx < nb
cnt = np.bincount(idx[sel], minlength=nb)
avg = lambda v: np.bincount(idx[sel], weights=v[sel], minlength=nb) / np.maximum(cnt, 1)
rb = avg(R); pO = avg(np.interp(R, O["r"], O["p"]))
core = (cnt > 0) & (rb <= r0 / 2)
dc = float((np.abs(avg(pc.numpy().ravel()) - pO)[core] / p0).max())
dfb = float((np.abs(avg(pf.numpy().ravel()) - pO)[core] / p0).max())
# ODE slope
Pb = kc.radial_bvp(60.0); n = len(Pb["fp"])
a_ode = -float(np.polyfit(Pb["r"][2:1 + n][3:60]**2, Pb["fp"][1:][3:60], 1)[1])
res = dict(h=h, L=L, E_over_12pi2=hist[-1] / (12 * math.pi**2), iterations=len(hist), a_lat=a_lat, a_ode=a_ode,
           a_rel=a_lat / a_ode - 1, p0_from_slope_rel=(a_lat**2 * (a_lat**2 - 1)) / (a_ode**2 * (a_ode**2 - 1)) - 1,
           core_dev_central=dc, core_dev_fb=dfb, virial_fb=(E4 - E2) / (E4 + E2), E2_fb=E2, E4_fb=E4,
           vonLaue_fb_rel=float(pf.sum() / pf.abs().sum()), seconds=time.time() - t0,
           note="POST-HOC, not pre-registered")
json.dump(res, open(os.path.join(OUT, f"press_posthoc_h{h:g}_L{L}.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
ok = cnt > 0
np.savez(os.path.join(OUT, f"press_posthoc_h{h:g}_L{L}.npz"), h=h, r=rb[ok], p_fb=avg(pf.numpy().ravel())[ok],
         p_c=avg(pc.numpy().ravel())[ok], p_ode=pO[ok])
