"""press_lattice.py - P4: pressure p(r) and shear s(r) measured on the relaxed lattice knot (B=1, h=0.2, L=48),
compared with the exact ODE profile (press_profiles.npz from press_ode.py).
Relaxation exactly as knot_k1 (L-BFGS on psi, phi = psi/|psi|, Dirichlet vacuum wall), 2 threads.
T_ij = 2 D_ij + 2 trD D_ij - 2 (D^2)_ij - delta_ij e,  D_ij from 4th-order central differences (kc.derivs).
Writes lab/out/press_lattice.json and lab/out/press_lattice.npz."""
import sys, os, json, math, time, numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot_core as kc
torch.set_num_threads(2)

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
h, L = 0.2, 48
t0 = time.time()
Pb = kc.radial_bvp(60.0)
X, Y, Z = kc.grid((L, L, L), h)
phi = kc.hedgehog(X, Y, Z, Pb)
E_start = float(kc.energy(phi, h))
hist = None
for chunk in range(6):
    phi, hist, done = kc.relax_lbfgs(phi, h, maxit=200, tol=1e-11, log=print, every=5, wall=60, hist=hist)
    print(f"chunk {chunk} iters {len(hist)} E/12pi2 {hist[-1]/(12*math.pi**2):.8f} done {done} t {time.time()-t0:.0f}s",
          flush=True)
    if done or time.time() - t0 > 400:
        break
E = float(kc.energy(phi, h))
b = kc.baryon_density(phi, h)
B = float(b.sum() * h**3)
cx = [float((b * G).sum() / b.sum()) for G in (X, Y, Z)]

with torch.no_grad():
    g = kc.derivs(phi, h)                                    # 3 x (4, L, L, L)
    Dm = torch.stack([torch.stack([(g[i] * g[j]).sum(0) for j in range(3)]) for i in range(3)])   # 3x3xLLL
    trD = Dm[0, 0] + Dm[1, 1] + Dm[2, 2]
    D2 = torch.einsum("ik...,kj...->ij...", Dm, Dm)
    trD2 = D2[0, 0] + D2[1, 1] + D2[2, 2]
    e2 = trD; e4 = 0.5 * (trD**2 - trD2); e = e2 + e4
    T = 2 * Dm + 2 * trD * Dm - 2 * D2
    for i in range(3):
        T[i, i] -= e
    p = (T[0, 0] + T[1, 1] + T[2, 2]) / 3
    x = [X - cx[0], Y - cx[1], Z - cx[2]]
    rr = torch.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    rs = torch.where(rr > 0, rr, torch.ones_like(rr))
    n = [xi / rs for xi in x]
    Trr = sum(T[i, j] * n[i] * n[j] for i in range(3) for j in range(3))
    s = 1.5 * (Trr - p)                                      # T_rr = p + 2s/3
    # quartic-term trace pointwise: e4 - (2 trD^2 - 2 trD2 - 3 e4) = 0
    tr4 = (e4 - (2 * trD**2 - 2 * trD2 - 3 * e4)).abs().max().item()
    Ip = float(p.sum() * h**3); Ipabs = float(p.abs().sum() * h**3)
    E2l = float(e2.sum() * h**3); E4l = float(e4.sum() * h**3)
    Ir2p = float((rr**2 * p).sum() * h**3); Ir2s = float((rr**2 * s).sum() * h**3)

# shell averages (bins of width h/2)
R = rr.numpy().ravel(); pv = p.numpy().ravel(); sv = s.numpy().ravel(); ev = e.numpy().ravel()
edges = np.arange(0, 4.7 + 1e-9, h / 2)
idx = np.digitize(R, edges) - 1
nb = len(edges) - 1
cnt = np.bincount(idx[idx < nb], minlength=nb)
avg = lambda v: np.bincount(idx[idx < nb], weights=v[idx < nb], minlength=nb) / np.maximum(cnt, 1)
rb = avg(R); pb = avg(pv); sb = avg(sv); eb = avg(ev)
# ODE evaluated at the SAME points, averaged the same way
O = np.load(os.path.join(OUT, "press_profiles.npz"))
pO = np.interp(R, O["r"], O["p"]); sO = np.interp(R, O["r"], O["s"]); eO = np.interp(R, O["r"], O["T00"])
pOb = avg(pO); sOb = avg(sO); eOb = avg(eO)
p0 = float(O["p"][0]); r0 = float(O["r0"])
ok = cnt > 0
core = ok & (rb <= r0 / 2)
dev = np.abs(pb - pOb)[core] / p0
# also: lattice p at its own zero crossing
sg = np.sign(pb[ok]); rk = rb[ok]
j = int(np.argmax(sg < 0)); r0_lat = float(np.interp(0, [pb[ok][j], pb[ok][j - 1]], [rk[j], rk[j - 1]]))
M = E
res = {"h": h, "L": L, "half_box": (L - 1) / 2 * h, "E": E, "E_over_12pi2": E / (12 * math.pi**2), "B": B,
       "E_start_over_12pi2": E_start / (12 * math.pi**2), "iterations": len(hist), "centroid": cx,
       "core_rmax": r0 / 2, "core_bins": int(core.sum()), "core_points": int(cnt[core].sum()),
       "core_max_dev_over_p0": float(dev.max()), "core_dev_per_bin": dev.tolist(), "core_r_bins": rb[core].tolist(),
       "p_lat_first_bin": float(pb[ok][0]), "p_ode_first_bin": float(pOb[ok][0]), "r_first_bin": float(rb[ok][0]),
       "r0_lattice": r0_lat, "r0_ode": r0,
       "vonLaue_lattice_rel": Ip / Ipabs, "E2_lat": E2l, "E4_lat": E4l, "virial_lat": (E4l - E2l) / (E2l + E4l),
       "D_e4_lattice_from_p": M * Ir2p / 4, "D_e4_lattice_from_s": -(4 / 15) * M * Ir2s / 4,
       "quartic_trace_pointwise_max": tr4, "seconds": time.time() - t0}
res["D_lattice_ANW_from_p"] = res["D_e4_lattice_from_p"] / 5.44**4
res["D_lattice_ANW_from_s"] = res["D_e4_lattice_from_s"] / 5.44**4
res["pass_a"] = bool(res["core_max_dev_over_p0"] < 0.05)
json.dump(res, open(os.path.join(OUT, "press_lattice.json"), "w"), indent=1)
np.savez(os.path.join(OUT, "press_lattice.npz"), r=rb[ok], p=pb[ok], s=sb[ok], T00=eb[ok], p_ode=pOb[ok],
         s_ode=sOb[ok], T00_ode=eOb[ok], cnt=cnt[ok])
print(json.dumps({k: v for k, v in res.items() if k not in ("core_dev_per_bin", "core_r_bins")}, indent=1))
