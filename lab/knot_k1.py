"""knot_k1.py - K1: one mass, size-independent. Relax a B=1 hedgehog on an L^3 lattice of balls.
usage: python knot_k1.py h L [crude]
Writes lab/out/knot_k1_<tag>.json and the relaxed field to the scratch dir (KNOT_SCRATCH) for K2/K4/K5.
"""
import sys, os, json, math, time, numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot_core as kc

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
SCR = os.environ.get("KNOT_SCRATCH", OUT)
h = float(sys.argv[1]); L = int(sys.argv[2]); crude = len(sys.argv) > 3 and sys.argv[3] == "crude"
tag = f"h{h:g}_L{L}" + ("_crude" if crude else "")
t0 = time.time()
P = kc.radial_profile()
X, Y, Z = kc.grid((L, L, L), h)
if crude:
    Q = dict(r=P["r"], f=np.pi * np.exp(-P["r"] / 0.6))
    phi = kc.hedgehog(X, Y, Z, Q)
else:
    phi = kc.hedgehog(X, Y, Z, P)
ck = os.path.join(SCR, f"knot_k1_state_{tag}.pt")
hist0 = None
if os.path.exists(ck):
    S = torch.load(ck); phi, hist0 = S["phi"], S["hist"]
E_start = float(kc.energy(phi, h)) if hist0 is None else S["E_start"]
# NOTE: the pre-registered arrested Newton flow (dt = 0.15 h) went unstable in the first production attempt (h = 0.2,
# L = 48: converged to 1.205782 at it 250, then a stiff lattice mode blew up). Replaced by L-BFGS on psi, phi = psi/|psi|.
phi, hist, done = kc.relax_lbfgs(phi, h, maxit=200, tol=1e-11, log=print, every=5, wall=300, hist=hist0)
if not done and len(hist) < 200:
    torch.save(dict(phi=phi, hist=hist, E_start=E_start), ck); print("CHECKPOINT", len(hist)); sys.exit(3)
if os.path.exists(ck):
    os.remove(ck)

E = hist[-1]
b = kc.baryon_density(phi, h)
B = float(b.sum() * h**3)
# residual force
p = phi.clone().requires_grad_(True)
gr = torch.autograd.grad(kc.energy(p, h), p)[0]
Fmax = float(kc.tangent(phi, gr).abs().max() / h**3)
# moment of inertia: rigid isorotation about isospin axis 3, w = 1: v = (0, -phi2, phi1, 0); T = v.pi/2 = Lam/2
v = torch.stack([torch.zeros_like(phi[0]), -phi[2], phi[1], torch.zeros_like(phi[0])])
pi = kc.momentum_from_velocity(phi, v, h)
T = float(0.5 * (v * pi).sum() * h**3)
Lam = 2 * T
# radial profile f(r) = arccos(sigma) (bins of h/2) for later boosts
R = torch.sqrt(X**2 + Y**2 + Z**2).numpy().ravel(); f = torch.arccos(phi[0].clamp(-1, 1)).numpy().ravel()
edges = np.arange(0, R.max() + h / 2, h / 2)
idx = np.digitize(R, edges)
rb = np.array([R[idx == i].mean() for i in range(1, len(edges)) if (idx == i).any()])
fb = np.array([f[idx == i].mean() for i in range(1, len(edges)) if (idx == i).any()])
# shape: second moments of baryon and energy density
ed = kc.static_density(phi, h)
mom = lambda w: dict(r2=float((w * (X**2 + Y**2 + Z**2)).sum() / w.sum()),
                     aniso=float(max((w * A**2).sum() for A in (X, Y, Z)) / min((w * A**2).sum() for A in (X, Y, Z)) - 1))
res = dict(h=h, L=L, half_box=(L - 1) / 2 * h, crude=crude, iterations=len(hist), E=E, E_over_12pi2=E / (12 * math.pi**2),
           E_start_over_12pi2=E_start / (12 * math.pi**2), B=B, Fmax=Fmax, Lambda=Lam, r_half_over_h=P["r_half"] / h,
           baryon=mom(b), energy=mom(ed), seconds=time.time() - t0,
           history=[hist[i] / (12 * math.pi**2) for i in range(0, len(hist), 25)])
print(json.dumps({k: v for k, v in res.items() if k != "history"}, indent=1))
os.makedirs(OUT, exist_ok=True)
json.dump(res, open(os.path.join(OUT, f"knot_k1_{tag}.json"), "w"), indent=1)
os.makedirs(SCR, exist_ok=True)
np.savez(os.path.join(SCR, f"knot_field_{tag}.npz"), phi=phi.numpy().astype(np.float32), h=h, rb=rb, fb=fb)
