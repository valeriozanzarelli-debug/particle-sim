"""nuclei_dev.py - DEV checks before the pre-registration (declared there): angular integrals of the rational maps,
rational-map profile energies vs HMS 1998, and shape measures of the UNRELAXED rational-map initial data on L=48, h=0.2."""
import math, json, os, time, numpy as np, torch
import nuclei_core as nc, knot_core as kc, knot2_core as k2
torch.set_num_threads(2)
k2.set_scheme("c4")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
L, h = 48, 0.2
X, Y, Z = kc.grid((L, L, L), h)
res = {}
for B in (1, 2, 3, 4):
    t0 = time.time()
    b, I = nc.angular_integrals(B, 400)
    r, f, E = nc.rm_profile(B, I, n=3000)
    phi = nc.rm_field(B, X, Y, Z, r, f)
    El = float(k2.energy(phi, h))
    w = kc.baryon_density(phi, h).numpy()
    tot, c, M, In = nc.moments(w, h)
    q = nc.multipoles(np.clip(w, 0, None), h, 4)
    d = dict(B_ang=b, I=I, E_ansatz_over=E / (12 * math.pi**2 * B), E_lat_over=El / (12 * math.pi**2 * B), B_lat=tot,
             M_eig=np.linalg.eigvalsh(M).tolist(), q=q, center_over_max=float(w[L//2-1:L//2+1, L//2-1:L//2+1, L//2-1:L//2+1].mean() / w.max()),
             asym={k: nc.asym(w, h, c, R) for k, R in dict(Cz90=nc.rot([0, 0, 1], math.pi/2), Cz45=nc.rot([0, 0, 1], math.pi/4),
                   Cz30=nc.rot([0, 0, 1], math.pi/6), Cx90=nc.rot([1, 0, 0], math.pi/2), C2x=nc.rot([1, 0, 0], math.pi),
                   C3=nc.rot([1, 1, 1], 2*math.pi/3), inv=-np.eye(3), Pz=np.diag([1., 1, -1])).items()}, sec=time.time() - t0)
    res[B] = d
    print(B, json.dumps(d), flush=True)
json.dump(res, open(os.path.join(OUT, "nuclei_dev.json"), "w"), indent=1)
