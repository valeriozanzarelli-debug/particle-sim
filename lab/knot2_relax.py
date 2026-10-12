"""knot2_relax.py - relax the B = 1 knot with a chosen scheme (L-BFGS, start = old relaxed field of knot_k1.py).
usage: python knot2_relax.py scheme h L     (resumable, chunks <= 300 s; field -> KNOT_SCRATCH/knot2_field_<tag>.npz)
"""
import sys, os, json, math, time, numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot2_core as k2

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
SCR = os.environ.get("KNOT_SCRATCH", OUT)
scheme, h, L = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
k2.set_scheme(scheme)
tag = f"{scheme}_h{h:g}_L{L}"
ck = os.path.join(SCR, f"knot2_relax_state_{tag}.pt")
X, Y, Z = k2.grid((L, L, L), h)
t0 = time.time()
if os.path.exists(ck):
    S = torch.load(ck); phi, hist0, E_start = S["phi"], S["hist"], S["E_start"]
else:
    src = os.path.join(SCR, f"knot_field_h{h:g}_L{L}.npz")
    if os.path.exists(src):
        phi = torch.tensor(np.load(src)["phi"], dtype=torch.float64)
    else:
        phi = k2.hedgehog(X, Y, Z, k2.radial_bvp())
    phi = k2.normalize(phi); hist0 = None; E_start = float(k2.energy(phi, h))
phi, hist, done = k2.relax_lbfgs(phi, h, maxit=400, tol=1e-11, log=print, every=5, wall=290, hist=hist0)
if not done and len(hist) < 400:
    torch.save(dict(phi=phi, hist=hist, E_start=E_start), ck); print("CHECKPOINT", len(hist), flush=True); sys.exit(3)
if os.path.exists(ck):
    os.remove(ck)
E = hist[-1]
b = k2.baryon_density(phi, h); B = float(b.sum() * h**3)
v = torch.stack([torch.zeros_like(phi[0]), -phi[2], phi[1], torch.zeros_like(phi[0])])
Lam = float((v * k2.mass_apply(phi, v, h)).sum() * h**3) * 2
R = torch.sqrt(X**2 + Y**2 + Z**2)
rB2 = float((b * R**2).sum() / b.sum())
Rn = R.numpy().ravel(); f = torch.arccos(phi[0].clamp(-1, 1)).numpy().ravel()
edges = np.arange(0, Rn.max() + h / 2, h / 2); idx = np.digitize(Rn, edges)
rb = np.array([Rn[idx == i].mean() for i in range(1, len(edges)) if (idx == i).any()])
fb = np.array([f[idx == i].mean() for i in range(1, len(edges)) if (idx == i).any()])
res = dict(scheme=scheme, h=h, L=L, half_box=(L - 1) / 2 * h, R_eff=(L + 1) / 2 * h, iterations=len(hist), E=E,
           E_over_12pi2=E / (12 * math.pi**2), E_start_over_12pi2=E_start / (12 * math.pi**2), B=B, Lambda=Lam, rB2=rB2,
           r_half=float(np.interp(np.pi / 2, fb[::-1], rb[::-1])), seconds=time.time() - t0)
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(OUT, f"knot2_relax_{tag}.json"), "w"), indent=1)
np.savez(os.path.join(SCR, f"knot2_field_{tag}.npz"), phi=phi.numpy().astype(np.float32), h=h, rb=rb, fb=fb)
