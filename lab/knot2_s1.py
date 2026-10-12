"""knot2_s1.py - S1: lattice dispersion of small pion waves (measured with the dynamics code), stencil mismatch on the
knot background, Cherenkov thresholds and the knot spectrum above them.  -> lab/out/knot2_s1.json
"""
import sys, os, json, math, time, numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot2_core as k2
from scipy.optimize import brentq, curve_fit

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
SCR = os.environ.get("KNOT_SCRATCH", OUT)
h, N, eps = 0.2, 64, 1e-3
t0 = time.time()
res = dict(h=h, N=N, eps=eps)


def measure_omega(scheme, m, diag):
    k2.set_scheme(scheme); k2.PERIODIC = True
    shape = (N, N, 2) if diag else (N, 2, 2)
    X, Y, Z = k2.grid(shape, h)
    kk = 2 * math.pi * m / (N * h)
    x0 = float(X.min())   # grid is centred on half-integers: phase from the first point so m = N/2 is not zero
    arg = kk * (X - x0 + Y - x0) if diag else kk * (X - x0)
    c = torch.cos(arg)
    phi = torch.stack([torch.sqrt(1 - (eps * c) ** 2), eps * c, torch.zeros_like(c), torch.zeros_like(c)])
    pi = torch.zeros_like(phi)
    kmag = kk * (math.sqrt(2) if diag else 1)
    w_an = (math.sqrt(2) if diag else 1) * k2.omega_1d(kk * h, scheme) / h
    dt = min(0.25 * h, 0.1 / w_an)          # omega dt <= 0.1: RK4 phase error ~ 1e-6
    T = 1.5 * 2 * math.pi / w_an
    n = int(math.ceil(T / dt))
    ts, A = [0.0], [1.0]
    for i in range(n):
        phi, pi = k2.rk4(phi, pi, h, dt)
        ts.append((i + 1) * dt); A.append(float((phi[1] * c).sum() / (c * c).sum() / eps))
    ts, A = np.array(ts), np.array(A)
    (A0, w), _ = curve_fit(lambda t, a, w: a * np.cos(w * t), ts, A, p0=[1.0, w_an])
    k2.PERIODIC = False
    return kmag, abs(w), w_an


for scheme in ("old", "c4"):
    for diag in (False,):          # diagonal dropped for time (declared): analytic omega^2 = sum_i omega_1d(k_i)^2 only
        rows = []
        for m in (range(1, N // 2 + 1) if not diag else range(2, N // 2 + 1, 2)):
            kmag, w, w_an = measure_omega(scheme, m, diag)
            rows.append((m, kmag, w, w_an))
        R = np.array(rows); kh = R[:, 1] * h; w = R[:, 2]
        vp = w / R[:, 1]
        vg = np.gradient(w, R[:, 1])
        out = dict(kh=kh.tolist(), omega=w.tolist(), omega_analytic=R[:, 3].tolist(), v_phase=vp.tolist(), v_group=vg.tolist(),
                   max_rel_dev_from_analytic=float(np.max(np.abs(w / R[:, 3] - 1))))
        for lab, lim in (("pi4", math.pi / 4 + 1e-9), ("pi2", math.pi / 2 + 1e-9), ("full", 10)):
            sel = (kh <= lim) if not diag else (kh / math.sqrt(2) <= lim)   # per-component k h for the diagonal
            iv = np.argmax(np.abs(vp[sel] - 1)); ig = np.argmax(np.abs(vg[sel] - 1))
            out[f"max_vphase_dev_{lab}"] = [float(abs(vp[sel][iv] - 1)), float(kh[sel][iv])]
            out[f"max_vgroup_dev_{lab}"] = [float(abs(vg[sel][ig] - 1)), float(kh[sel][ig])]
        out["code_check_pass"] = out["max_rel_dev_from_analytic"] < 0.005
        res[f"S1a_{scheme}_{'diag' if diag else 'axis'}"] = out
        print(scheme, "diag" if diag else "axis", {k: v for k, v in out.items() if k.startswith("max") or k == "code_check_pass"},
              f"{time.time()-t0:.0f}s", flush=True)

# ---- S1b: local light speed along x on the old rest knot (x-gradient-only perturbations)
F = np.load(os.path.join(SCR, "knot_field_h0.2_L48.npz"))
phi = torch.tensor(F["phi"], dtype=torch.float64)
k2.set_scheme("old")
gp, gm = k2.derivs_fb(phi, h)
P = phi.numpy()
I4 = np.eye(4)
K = np.zeros(P.shape[1:] + (4, 4)); M = np.zeros_like(K)
a = 1.0
for j in (1, 2):
    for g in (gp[j].numpy(), gm[j].numpy()):
        gg = np.moveaxis(g, 0, -1)
        K += 0.5 * ((gg * gg).sum(-1)[..., None, None] * I4 - gg[..., :, None] * gg[..., None, :])
        a = a + 0.5 * (gg * gg).sum(-1)
    gb = np.moveaxis(0.5 * (gp[j] + gm[j]).numpy(), 0, -1)
    ph = np.moveaxis(P, 0, -1)
    gt = gb - ph * (ph * gb).sum(-1, keepdims=True)
    M -= gt[..., :, None] * gt[..., None, :]
K += I4; M += a[..., None, None] * I4
ph = np.moveaxis(P, 0, -1); Pr = I4 - ph[..., :, None] * ph[..., None, :]; pp = ph[..., :, None] * ph[..., None, :]
Kp = Pr @ K @ Pr + pp; Mp = Pr @ M @ Pr + pp
Lc = np.linalg.cholesky(Mp); Li = np.linalg.inv(Lc)
ev = np.linalg.eigvalsh(Li @ Kp @ np.swapaxes(Li, -1, -2))
c2min = ev[..., 0]
X, Y, Z = k2.grid(phi.shape[1:], h)
Rr = torch.sqrt(X**2 + Y**2 + Z**2).numpy()
i = np.unravel_index(np.argmin(c2min), c2min.shape)
# energy-weighted mean deficit
e = k2.static_density(phi, h).numpy()
res["S1b_old_rest_knot"] = dict(c_eff_min=float(np.sqrt(c2min.min())), at_r=float(Rr[i]),
                                c_eff_energy_weighted=float((np.sqrt(c2min) * e).sum() / e.sum()),
                                note="c2 and c4: c_eff = 1 identically (kinetic and static built with the same stencil)")
print("S1b", res["S1b_old_rest_knot"], flush=True)

# ---- S1c/d: Cherenkov thresholds and knot spectrum above them
th = {}
for scheme in ("old", "c4"):
    for v in (0.3, 0.6, 0.8):
        f = lambda x: k2.omega_1d(x, scheme) / x - v
        th[f"{scheme}_v{v}"] = brentq(f, 1e-3, math.pi) if f(math.pi) < 0 else None
res["S1c_cherenkov_kh"] = th
spec = {}
for hh, L, NX in ((0.2, 48, 112), (0.15, 64, 150)):
    F = np.load(os.path.join(SCR, f"knot_field_h{hh:g}_L{L}.npz"))
    rb, fb = F["rb"], F["fb"]
    prof = dict(r=np.concatenate([[0.0], rb, [rb[-1] + hh, 1e3]]), f=np.concatenate([[np.pi], fb, [0.0, 0.0]]))
    X, Y, Z = k2.grid((NX, L, L), hh)
    for v in (0.6, 0.8):
        g = 1 / math.sqrt(1 - v * v)
        q = k2.hedgehog(X, Y, Z, prof, gx=g).numpy()
        Fq = np.abs(np.fft.rfft(q[1:], axis=1)) ** 2
        pw = Fq.sum(axis=(0, 2, 3)); kx = np.fft.rfftfreq(NX) * 2 * math.pi
        for scheme in ("old", "c4"):
            kc_ = th[f"{scheme}_v0.8"]
            spec[f"h{hh:g}_v{v}_{scheme}_above_kh{kc_:.2f}"] = float(pw[kx >= kc_].sum() / pw.sum())
        spec[f"h{hh:g}_v{v}_above_kh_pi2"] = float(pw[kx >= math.pi / 2].sum() / pw.sum())
res["S1d_knot_spectrum_fraction"] = spec
print("S1c", th, "\nS1d", spec)
res["seconds"] = time.time() - t0
json.dump(res, open(os.path.join(OUT, "knot2_s1.json"), "w"), indent=1)
