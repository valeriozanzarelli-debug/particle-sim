"""kk_shape_spectra.py - TEST 1: rest-mass spectra of the zig-zag ball on hidden SHAPES, continuum checks,
circle x shape product, and charge conservation in each case. Criteria: lab/out/kk_shape_prereg.json.
Writes lab/out/kk_shape_spectra.json and kk_shape_modes.npz (mode patterns for the figure)."""
import json, time, sys
import numpy as np
from scipy.special import jn_zeros
from scipy.sparse.linalg import eigsh
import kk_shape_core as S

t0 = time.time()
out = {}
rng = np.random.default_rng(S.SEED)
SH = [("disk", None), ("ellipse", 1.5), ("ellipse", 2.0), ("ellipse", 3.0), ("triangle", None), ("stadium", 1.0)]


def name(k, p):
    return k if p is None else f"{k}_{p:g}"


def theta_map(kind, R, par=None, w=0):
    """w = 0: sharp mass wall (pre-registered). w > 0: smooth wall over w sites (POST-HOC variant, added after
    seeing that the sharp wall mixes the physical Dirac cone with the lattice doubler cone at (pi/2, pi/2))"""
    mask, info = S.shape_mask(kind, np.pi * R**2, par, margin=4 + w)
    return S.soft_wall(mask, w), mask, info


# ---------------- 1a. disk: walk vs continuum infinite-mass Dirac (Berry-Mondragon) ----------------
cont = np.array([k for k, l, s in S.dirac_disk_levels(1.0, 9.0)])
out["disk_check"] = {}
for R in (10, 20, 40):
    th, mask, info = theta_map("disk", R)
    Reff = np.sqrt(info["sites"] / np.pi)
    E = S.rest_masses_sparse(th, k=16)
    mid, spl, worst = S.pair_up(E)
    x = mid[:6] * Reff
    err = (x - cont[:6]) / cont[:6]
    out["disk_check"][R] = dict(R_eff=Reff, x_walk=x.tolist(), x_cont=cont[:6].tolist(), rel_err=err.tolist(),
                                pair_split_x=(spl[:6] * Reff).tolist(), worst_split_over_gap=worst)
    print(f"[disk R={R}] E*R_eff (pair midpoints) {np.round(x, 3)}  continuum {np.round(cont[:6], 3)}  "
          f"max|err| {np.max(np.abs(err)):.3f}  doubler splitting*R {np.round(spl[:3] * Reff, 3)}  ({time.time() - t0:.0f}s)")
    sys.stdout.flush()
errs = [np.max(np.abs(out["disk_check"][R]["rel_err"])) for R in (10, 20, 40)]
out["disk_check"]["PASS"] = bool(errs[1] < 0.05 and errs[0] > errs[1] > errs[2])
x20 = np.array(out["disk_check"][20]["x_walk"]); x40 = np.array(out["disk_check"][40]["x_walk"])
xr = 2 * x40 - x20
out["disk_check"]["richardson_posthoc"] = dict(x=xr.tolist(), rel_err=((xr - cont[:6]) / cont[:6]).tolist())
print("  POST-HOC Richardson 2 x(40) - x(20):", np.round(xr, 3), "rel err", np.round((xr - cont[:6]) / cont[:6], 4))
print("  disk PASS (err<5% at R=20 and decreasing):", out["disk_check"]["PASS"], np.round(errs, 4))

# ---------------- 1b. all shapes: spectra at R = 20 and 40 (convergence), mode patterns ----------------
out["shapes"] = {}
modes = {}
for kind, par in SH:
    nm = name(kind, par)
    res = {}
    for R in (20, 40):
        th, mask, info = theta_map(kind, R, par)
        Reff = np.sqrt(info["sites"] / np.pi)
        if R == 20:
            E, V, sites = S.rest_masses_sparse(th, k=16, vectors=True)
        else:
            E = S.rest_masses_sparse(th, k=16)
        mid, spl, worst = S.pair_up(E)
        res[R] = dict(x=(mid * Reff).tolist(), split_x=(spl * Reff).tolist(), worst=worst)
        if R == 20:
            o = np.argsort(E)
            # first 4 levels (each a doubler pair)
            dens = [S.mode_density(V[:, o[2 * i]], sites, *th.shape) + S.mode_density(V[:, o[2 * i + 1]], sites, *th.shape)
                    for i in range(4)]   # density of the doubler PAIR (removes the physical x doubler interference)
            modes[nm] = dict(mask=mask, dens=np.array(dens), x=(mid[:4] * Reff))
    x20, x40 = np.array(res[20]["x"][:6]), np.array(res[40]["x"][:6])
    conv = float(np.max(np.abs(x40 - x20) / x40))
    m = np.array(res[40]["x"][:3])
    out["shapes"][nm] = dict(R20=res[20], R40=res[40], conv_20_40=conv, m2_m1=float(m[1] / m[0]),
                             m3_m1=float(m[2] / m[0]), koide=S.koide(m))
    print(f"[{nm:10s}] E*R_eff (R=40) {np.round(x40, 3)}  change 20->40 {conv:.3f}  m2/m1 {m[1] / m[0]:.3f} "
          f"m3/m1 {m[2] / m[0]:.3f} Q {S.koide(m):.4f}  ({time.time() - t0:.0f}s)")
    sys.stdout.flush()

# ---------------- 1c. Laplacian (Dirichlet, 5-point) vs exact: disk and equilateral triangle ----------------
out["laplacian_check"] = {}
R = 30
mask, info = S.shape_mask("disk", np.pi * R**2)
L, _ = S.laplacian(mask)
lam = np.sort(eigsh(L.tocsc(), k=8, sigma=0, which="LM")[0])
Reff = np.sqrt(info["sites"] / np.pi)
ex = sorted([(jn_zeros(l, 3)[m] / Reff)**2 for l in range(0, 6) for m in range(3) for _ in range(1 if l == 0 else 2)])[:8]
e_d = (lam[:6] - np.array(ex[:6])) / np.array(ex[:6])
mask, info = S.shape_mask("triangle", np.pi * R**2)
L, _ = S.laplacian(mask)
lamt = np.sort(eigsh(L.tocsc(), k=8, sigma=0, which="LM")[0])
s_eff = np.sqrt(4 * info["sites"] / np.sqrt(3))
ext = sorted([16 * np.pi**2 / (9 * s_eff**2) * (m * m + m * n + n * n) for m in range(1, 8) for n in range(1, 8)])[:8]
e_t = (lamt[:6] - np.array(ext[:6])) / np.array(ext[:6])
lam60 = {}
for kk in ("disk", "triangle"):
    m60, i60 = S.shape_mask(kk, np.pi * 60**2)
    L60, _ = S.laplacian(m60)
    lam60[kk] = np.sort(eigsh(L60.tocsc(), k=8, sigma=0, which="LM")[0])
R60 = np.sqrt(S.shape_mask("disk", np.pi * 60**2)[1]["sites"] / np.pi)
ex60 = np.array(ex[:6]) * Reff**2 / R60**2
e_d60 = (lam60["disk"][:6] - ex60) / ex60
s60 = np.sqrt(4 * S.shape_mask("triangle", np.pi * 60**2)[1]["sites"] / np.sqrt(3))
ext60 = np.array(ext[:6]) * s_eff**2 / s60**2
e_t60 = (lam60["triangle"][:6] - ext60) / ext60
print(f"[Laplacian R=60] disk rel err {np.round(e_d60, 4)}  triangle {np.round(e_t60, 4)}")
out["laplacian_check"] = dict(disk_rel_err_R60=e_d60.tolist(), triangle_rel_err_R60=e_t60.tolist(), disk_rel_err=e_d.tolist(), triangle_rel_err=e_t.tolist(),
                              disk_ratio_l2_l1=float(lam[1] / lam[0]), triangle_ratio=float(lamt[1] / lamt[0]),
                              PASS=bool(np.max(np.abs(e_d)) < 0.03 and np.max(np.abs(e_t)) < 0.03))
print(f"[Laplacian R=30] disk rel err {np.round(e_d, 4)}  triangle rel err {np.round(e_t, 4)}  PASS {out['laplacian_check']['PASS']}")

# ---------------- 1d. circle x shape: charge sectors ----------------
# without the M block (it does not commute with the circle hop) levels come in clusters of 4
# (2 x residual doubler pair); compare cluster means: E_n^2 =? k0^2 + E_{n=0}^2
out["product"] = {}
N0 = 64
for kind, par in (("disk", None), ("stadium", 1.0), ("triangle", None)):
    th, mask, info = theta_map(kind, 20, par)
    cl = lambda E: np.sort(E)[:8].reshape(2, 4).mean(1)
    E0 = cl(S.rest_masses_sparse(th, k=8, n=0, block=False))
    row = {}
    for n in (1, 2):
        k0 = 2 * np.pi * n / N0
        En = np.sort(S.rest_masses_sparse(th, k=8, n=n, N0=N0, block=False))
        pred = np.sqrt(k0**2 + E0**2)
        row[n] = dict(k0=k0, E=cl(En).tolist(), lowest=float(En[0]), pred=pred.tolist(),
                      rel_err=((cl(En) - pred) / pred).tolist())
    out["product"][name(kind, par)] = dict(E_n0=E0.tolist(), sectors=row)
    print(f"[circle(N0={N0}) x {name(kind, par)}] n=0 {np.round(E0, 4)}  n=1 {np.round(row[1]['E'], 4)} vs "
          f"sqrt(k0^2+E0^2) {np.round(row[1]['pred'], 4)}  rel err n=1 {np.round(row[1]['rel_err'], 4)} "
          f"n=2 {np.round(row[2]['rel_err'], 4)}")
    sys.stdout.flush()
# lattice caveat seen while debugging (N0 = 24, k0 a = 0.26): charged levels fall BELOW k0 and are boundary-
# enhanced; recorded here as a number
th, mask, info = theta_map("disk", 20)
out["product_caveat_N0_24_n1_lowest"] = float(np.sort(S.rest_masses_sparse(th, k=4, n=1, N0=24, block=False))[0])
print(f"  caveat: N0=24, n=1 lowest level {out['product_caveat_N0_24_n1_lowest']:.4f} < k0 = {2*np.pi/24:.4f} (lattice wall effect)")

# ---------------- 1e. charge conservation ----------------
out["charge"] = {}
# (i) hidden space = ONE loop (charge = momentum around it); shape only through curvature
N, T, th0 = 96, 400, 0.1
s = np.arange(N)
loops = [("circle", None), ("ellipse", 2.0), ("trefoil", 0.2), ("stadium", 1.0)]
for kind, par in loops:
    kap, _ = S.loop_sample(kind, par, N)
    kap_lat = kap / N                                        # curvature in 1/sites
    for lab, thv, V in (("mass=curvature", th0 * kap / np.mean(kap), None),
                        ("daCosta", np.full(N, th0), -kap_lat**2 / (8 * th0))):
        U = S.loop_U(thv, V)
        psi = np.concatenate([np.exp(2j * np.pi * s / N), np.exp(2j * np.pi * s / N)]) / np.sqrt(2 * N)
        nv = np.rint(np.fft.fftfreq(N) * N)
        def q(p):
            F = np.abs(np.fft.fft(p.reshape(2, N), axis=1))**2
            w = F.sum(0); return float(np.sum(nv * w) / w.sum())
        q0 = q(psi); dr = 0.0
        for t in range(T):
            psi = U @ psi; dr = max(dr, abs(q(psi) - q0))
        out["charge"][f"loop {kind} {lab}"] = dr
        print(f"[charge, loop only] {kind:8s} {lab:15s} max drift of n over {T} ticks: {dr:.3e}")
# (ii) hidden space = 2-D shape only; candidate charge J_z (sharp wall = pre-registered; soft wall w=8 = post-hoc)
for w in (0, 8):
  for kind, par in SH[:1] + SH[2:3] + SH[4:]:
    th, mask, info = theta_map(kind, 20, par, w)
    n1, n2 = th.shape
    c1, c2 = info["c1"], info["c2"]
    X1, X2 = np.meshgrid(np.arange(n1) - c1, np.arange(n2) - c2, indexing="ij")
    if kind == "triangle":                                   # centroid of the lattice triangle
        c1, c2 = c1 + X1[mask].mean(), c2 + X2[mask].mean()
        X1, X2 = np.meshgrid(np.arange(n1) - c1, np.arange(n2) - c2, indexing="ij")
    r = np.hypot(X1, X2); ph = np.arctan2(X2, X1)
    f = np.exp(-(r - 8)**2 / (2 * 3.0**2)) * np.exp(1j * ph) * mask
    psi = np.zeros((4, 1, n1, n2), complex)
    for c in range(4):
        psi[c, 0] = f / 2                                    # spinor (1,1) in L and in R: sigma_x = +1
    W = S.ShapeWalk(th, 1)
    J0 = S.Jz(psi, c1, c2); dr = 0.0
    for t in range(T):
        psi = W.tick(psi); dr = max(dr, abs(S.Jz(psi, c1, c2) - J0))
    out["charge"][f"plane w={w} {name(kind, par)}"] = dict(J0=J0, drift=dr)
    print(f"[charge, shape only, wall w={w}] {name(kind, par):10s} J_z(0) = {J0:.3f}  max drift {dr:.3f}")
    sys.stdout.flush()
# (iii) circle (8 sites) x shape, plus a random potential that depends on the shape coordinates only
for kind, par in SH[:1] + SH[2:3] + SH[4:]:
    th, mask, info = theta_map(kind, 12, par)
    n1, n2 = th.shape; N0 = 8
    V = (0.2 * rng.random((n1, n2)))[None]
    psi = (rng.normal(size=(4, 1, n1, n2)) + 1j * rng.normal(size=(4, 1, n1, n2))) * mask
    psi = psi * np.exp(2j * np.pi * np.arange(N0) / N0)[None, :, None, None]
    psi /= np.linalg.norm(psi)
    W = S.ShapeWalk(th, N0, V=V)
    q0 = S.circle_charge(psi); dr = 0.0
    for t in range(T):
        psi = W.tick(psi); dr = max(dr, abs(S.circle_charge(psi) - q0))
    out["charge"][f"product {name(kind, par)}"] = dr
    print(f"[charge, circle x shape] {name(kind, par):10s} n(0) = {q0:.3f}  max drift {dr:.2e}")

out["seconds"] = time.time() - t0
with open("out/kk_shape_spectra.json", "w") as fh:
    json.dump(out, fh, indent=1, default=float)
np.savez_compressed("out/kk_shape_modes.npz", **{f"{k}__mask": v["mask"] for k, v in modes.items()},
                    **{f"{k}__dens": v["dens"] for k, v in modes.items()}, **{f"{k}__x": v["x"] for k, v in modes.items()})
print(f"done {time.time() - t0:.0f}s")
