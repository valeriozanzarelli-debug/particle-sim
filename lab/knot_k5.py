"""knot_k5.py - K5: draw only the knot (isosurface of the baryon density), count primitives vs balls."""
import os, sys, json, numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot_core as kc

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
SCR = os.environ.get("KNOT_SCRATCH", OUT)
F = np.load(os.path.join(SCR, "knot_field_h0.2_L64.npz"))
phi = torch.tensor(F["phi"], dtype=torch.float64); h = float(F["h"]); N = phi.shape[1]
b = kc.baryon_density(phi, h).numpy()
lev = 0.1 * b.max()
org = -(N - 1) / 2 * h
T = kc.isosurface(b, lev, h, (org, org, org))
V = T.reshape(-1, 3)
c = V.mean(0)
# principal radii from the triangle-area weighted vertex cloud
area = 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1)
cen = T.mean(1)
C = (area[:, None, None] * np.einsum("ti,tj->tij", cen - c, cen - c)).sum(0) / area.sum()
ev = np.sqrt(np.linalg.eigvalsh(C))
balls = N**3
# baryon number enclosed by the drawn surface
Bin = float(b[b > lev].sum() * h**3)
res = dict(level=float(lev), triangles=int(len(T)), balls=balls, ratio=len(T) / balls, centroid=c.tolist(),
           centroid_err_over_h=float(np.linalg.norm(c) / h), principal_radii=ev.tolist(), aniso=float(ev.max() / ev.min() - 1),
           area=float(area.sum()), B_inside=Bin)
res["pass"] = dict(K5a=res["ratio"] < 0.01, K5b=res["centroid_err_over_h"] < 1 and res["aniso"] < 0.05)
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(OUT, "knot_k5.json"), "w"), indent=1)
np.save(os.path.join(SCR, "knot_k5_tris.npy"), T.astype(np.float32))
