"""knot2_s3.py - S3: moment of inertia Lambda on the lattice with the analytic 1/R tail correction (pre-registered),
and |mu_p/mu_n| from it.  -> lab/out/knot2_s3.json
Linear tail: outside the core each pion component is harmonic; in a cube with vacuum walls pi_a = C (D_a - H_a),
D_a = x_a/r^3, H_a harmonic with H_a = D_a on the walls. Missing inertia dLambda = (4/3) C^2 Q / R_eff with
Q = R [ int_out |D|^2 + int_cube (2 D.H - |H|^2) ]  (sphere: 4 pi 9/5).
"""
import sys, os, json, math, numpy as np, torch
from scipy.fft import dstn, idstn
from scipy.interpolate import RegularGridInterpolator
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot2_core as k2

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
SCR = os.environ.get("KNOT_SCRATCH", OUT)
LAM_ODE = 106.76962149426222; C_ODE = 2.159644468019555
res = {}


def unit_cube_H(N):
    """H_1 on the unit cube [-1,1]^3 (interior N^3 grid, spacing 2/(N+1)); H_1 = D_1 = x/r^3 on the walls."""
    d = 2 / (N + 1); x = -1 + d * np.arange(1, N + 1)
    D1 = lambda X, Y, Z: X / (X * X + Y * Y + Z * Z) ** 1.5
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    b = np.zeros((N, N, N))
    # boundary neighbours of the 7-point Laplacian
    for ax in range(3):
        for side, xb in ((0, -1.0), (-1, 1.0)):
            sl = [slice(None)] * 3; sl[ax] = side
            C = [X[tuple(sl)], Y[tuple(sl)], Z[tuple(sl)]]; C[ax] = np.full_like(C[ax], xb)
            b[tuple(sl)] -= D1(*C) / d**2
    lam1 = (2 * np.cos(np.pi * np.arange(1, N + 1) / (N + 1)) - 2) / d**2
    lam = lam1[:, None, None] + lam1[None, :, None] + lam1[None, None, :]
    u = idstn(dstn(b, type=1) / lam, type=1)
    return x, u, X, Y, Z


def Q_cube(N):
    x, H, X, Y, Z = unit_cube_H(N)
    d = x[1] - x[0]
    r = np.sqrt(X * X + Y * Y + Z * Z); r[r == 0] = np.inf
    D1 = X / r**3
    inner = 3 * ((2 * D1 * H - H * H).sum() * d**3)          # 3 components by cubic symmetry
    return inner, x, H


# outside integral: int dOmega max|n_i|  (R = 1)
rng = np.random.default_rng(1); n = rng.normal(size=(4_000_000, 3)); n /= np.linalg.norm(n, axis=1)[:, None]
out_int = 4 * math.pi * float(np.abs(n).max(1).mean())
qs = {}
for N in (63, 127):
    inner, x, H = Q_cube(N); qs[N] = out_int + inner
Q = qs[127] + (qs[127] - qs[63]) / 3            # 2nd-order Richardson in grid spacing
x, H = Q_cube(127)[1:]
res["Q"] = dict(out=out_int, Q63=qs[63], Q127=qs[127], Q=Q, Q_sphere=4 * math.pi * 9 / 5, kappa_vs_sphere=Q / (4 * math.pi * 9 / 5))
print("Q", res["Q"], flush=True)
H1 = RegularGridInterpolator((x, x, x), H, bounds_error=False, fill_value=None)


def H_vec(Xp, Yp, Zp, R):
    """H_a(x; R) = H_a(x/R; 1)/R^2, components by cyclic permutation of H_1."""
    pts = lambda A, B, C: np.stack([A.ravel(), B.ravel(), C.ravel()], 1) / R
    h1 = H1(pts(Xp, Yp, Zp)).reshape(Xp.shape) / R**2
    h2 = H1(pts(Yp, Zp, Xp)).reshape(Xp.shape) / R**2
    h3 = H1(pts(Zp, Xp, Yp)).reshape(Xp.shape) / R**2
    return [h1, h2, h3]


def corrected(scheme_tag, field_file, h, L, Lam, r_fit=(2.5, 3.5)):
    F = np.load(field_file); phi = F["phi"].astype(np.float64)
    X, Y, Z = [a.numpy() for a in k2.grid((L, L, L), h)]
    Reff = (L + 1) / 2 * h
    r = np.sqrt(X * X + Y * Y + Z * Z)
    sel = (r >= r_fit[0]) & (r <= r_fit[1])
    Hs = H_vec(X[sel], Y[sel], Z[sel], Reff)
    Ds = [A[sel] / r[sel] ** 3 for A in (X, Y, Z)]
    num = sum((phi[a + 1][sel] * (Ds[a] - Hs[a])).sum() for a in range(3))
    den = sum(((Ds[a] - Hs[a]) ** 2).sum() for a in range(3))
    C = num / den
    resid = math.sqrt(sum(((phi[a + 1][sel] - C * (Ds[a] - Hs[a])) ** 2).sum() for a in range(3)) /
                      sum((phi[a + 1][sel] ** 2).sum() for a in range(3)))
    dL = 4 / 3 * C * C * Q / Reff
    return dict(tag=scheme_tag, h=h, L=L, R_eff=Reff, Lambda_lat=Lam, C_fit=C, fit_rel_residual=resid, dLambda=dL,
                Lambda_corr=Lam + dL, dev_from_ODE=(Lam + dL) / LAM_ODE - 1, dLambda_with_C_ODE=4 / 3 * C_ODE**2 * Q / Reff)


rows = {}
for h, L in ((0.2, 48), (0.2, 64), (0.15, 64)):
    j = json.load(open(os.path.join(OUT, f"knot2_relax_c4_h{h:g}_L{L}.json")))
    rows[f"c4_h{h:g}_L{L}"] = corrected("c4", os.path.join(SCR, f"knot2_field_c4_h{h:g}_L{L}.npz"), h, L, j["Lambda"]) | dict(E=j["E"], rB2=j["rB2"])
for h, L in ((0.2, 48), (0.2, 64), (0.15, 64)):
    j = json.load(open(os.path.join(OUT, f"knot_k1_h{h:g}_L{L}.json")))
    rows[f"old_h{h:g}_L{L}"] = corrected("old", os.path.join(SCR, f"knot_field_h{h:g}_L{L}.npz"), h, L, j["Lambda"]) | dict(E=j["E"], rB2=j["baryon"]["r2"])
for k, v in rows.items():
    print(k, {a: round(b, 4) if isinstance(b, float) else b for a, b in v.items()}, flush=True)
res["rows"] = rows
prim = rows["c4_h0.15_L64"]
res["primary"] = dict(Lambda_corr=prim["Lambda_corr"], dev=prim["dev_from_ODE"], pass_=abs(prim["dev_from_ODE"]) < 0.02)
rich = lambda a, b: b + (b - a) * 0.15**2 / (0.2**2 - 0.15**2)
res["secondary_old_richardson"] = dict(value=rich(rows["old_h0.2_L48"]["Lambda_corr"], rows["old_h0.15_L64"]["Lambda_corr"]))
res["secondary_old_richardson"]["dev"] = res["secondary_old_richardson"]["value"] / LAM_ODE - 1
# fit Lambda(R) = Lam_inf - a/R - b/R^3, c4 h = 0.2 L = 32, 48, 64
pts = []
for L in (32, 48, 64):
    j = json.load(open(os.path.join(OUT, f"knot2_relax_c4_h0.2_L{L}.json"))); pts.append(((L + 1) / 2 * 0.2, j["Lambda"]))
Rr = np.array([p[0] for p in pts]); Lm = np.array([p[1] for p in pts])
sol = np.linalg.solve(np.c_[np.ones(3), -1 / Rr, -1 / Rr**3], Lm)
res["secondary_fit_c4_h0.2"] = dict(R=Rr.tolist(), Lambda=Lm.tolist(), Lam_inf=sol[0], a=sol[1], b=sol[2],
                                     a_theory=4 / 3 * rows["c4_h0.2_L64"]["C_fit"] ** 2 * Q, dev=sol[0] / LAM_ODE - 1)

# |mu_p/mu_n| (ANW, one input M_Delta/M_N fixes e)
MN, MD = 938.92, 1232.0


def anw(E, Lam, rB2):
    lam = 3 / (2 * (MD - MN)); M = MN - 3 / (8 * lam)
    e = (E * Lam / (4 * M * lam)) ** 0.25; Fpi = 4 * M * e / E
    r2 = rB2 * 4 / (e * Fpi) ** 2
    mu_p = r2 * MN * (MD - MN) / 9 + MN / (2 * (MD - MN)); mu_n = r2 * MN * (MD - MN) / 9 - MN / (2 * (MD - MN))
    return dict(e=e, F_pi=Fpi, mu_p=mu_p, mu_n=mu_n, mu_ratio=abs(mu_p / mu_n))


P = k2.radial_bvp()
res["mu"] = dict(lattice_primary=anw(prim["E"], prim["Lambda_corr"], prim["rB2"]), ode=anw(P["E"], P["Lam"], P["rB2"]),
                 lattice_uncorrected=anw(prim["E"], prim["Lambda_lat"], prim["rB2"]), experiment=2.7928 / 1.9130)
res["mu"]["pass_"] = abs(res["mu"]["lattice_primary"]["mu_ratio"] / res["mu"]["experiment"] - 1) < 0.03
res["mu"]["dev_exp"] = res["mu"]["lattice_primary"]["mu_ratio"] / res["mu"]["experiment"] - 1
res["mu"]["dev_ode"] = res["mu"]["lattice_primary"]["mu_ratio"] / res["mu"]["ode"]["mu_ratio"] - 1
cv = lambda o: bool(o) if isinstance(o, (bool, np.bool_)) else float(o)
print(json.dumps({k: v for k, v in res.items() if k != "rows"}, indent=1, default=cv))
json.dump(res, open(os.path.join(OUT, "knot2_s3.json"), "w"), indent=1, default=cv)
