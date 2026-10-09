"""dirac_surf2d.py - 2+1 D surfer tests of the Dirac surfer rule (split from dirac_free.py; same tolerances).
  S1b |v| <= c;  S2c/S2d equivariance (chi^2 and ensemble mean);  central surfers: loop radius and speed;
  S3 2D drift = group velocity;  G2 surfers of a spin-up packet rotate with Omega = 1/(2 m sigma^2).
Grid 256 x 256, L = 72 lambda_C (dx = 0.28), sigma = 4; main run T = 24 hbar/mc^2 (7.6 ZBW periods).
Run: timeout 590 python lab/dirac_surf2d.py
"""
import json, os, time
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2 as chi2dist
from dirac_core import Dirac, gaussian_packet, SEED, LAMBDA_C, T_UNIT, DT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TOL = dict(zbw_freq_rel=0.02, zbw_amp_rel=0.05, zbw_proj_abs=1e-6, vmax=1 + 1e-12,
           quantile_1d=1e-3, chi2_p=1e-3, mean_sigma=4.0, drift_rel=0.01, g_rel=0.01,
           omega_rel=0.05)
OUT = os.path.join(os.path.dirname(__file__), "out")
res, verdict = {"2d": {}}, []
t_start = time.time()


def check(name, ok, msg):
    verdict.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def zbw_fit(t, x, tmax):
    sel = t <= tmax
    f = lambda t, c0, c1, A, w, ph: c0 + c1 * t + A * np.sin(w * t + ph)
    best = None
    for ph0 in (0, 1.5, 3, 4.5):
        try:
            p, _ = curve_fit(f, t[sel], x[sel], p0=[x[0], 0.0, 0.5, 2.0, ph0], maxfev=20000)
        except RuntimeError:
            continue
        r = np.sum((f(t[sel], *p) - x[sel])**2)
        if best is None or r < best[1]:
            best = (p, r)
    p = best[0]
    return abs(p[2]), abs(p[3]), p


rng = np.random.default_rng(SEED + 1)
D2 = Dirac(256, 72.0, 2)
sig2 = 4.0
X, Y = D2.x
ts2 = np.arange(0, 24 + 1e-9, 0.05)
psi_m = gaussian_packet(D2, sig2, [0, 0], [1, 1])
pkm = D2.fft(psi_m)
mean_tr = {"mixed_k0": (np.array([D2.mean(D2.evolve(pkm, t), X) for t in ts2]),
                        np.array([D2.mean(D2.evolve(pkm, t), Y) for t in ts2]))}
# surfers, mixed packet at rest
N2 = 3000
psi_m = gaussian_packet(D2, sig2, [0, 0], [1, 1])
pos0 = D2.birth(psi_m, N2, rng)
t0 = time.time()
times2, traj2, vmax2 = D2.run_surfers(psi_m, pos0, 24.0, record_every=1)
print(f"2D surfers: {N2} x {len(times2) - 1} steps in {time.time() - t0:.1f} s", flush=True)
check("S1b 2D |v| <= c", vmax2 <= TOL["vmax"], f"max |v| = {vmax2:.12f}")


def chi2_test(psi, pos, nb=12, rmax=3.0):
    """bins in x and y over [-rmax sigma, rmax sigma]^2 plus one overflow bin"""
    rho = D2.density_current(psi)[0] * D2.dx**2
    # (fix after first run: edges aligned with grid-cell boundaries; with arbitrary edges the expected mass of a
    #  cell cut by an edge was given entirely to one bin, which alone made chi2 fail at t = 0, before any motion)
    edges = (np.round(np.linspace(-rmax * sig2, rmax * sig2, nb + 1) / D2.dx) + 0.5) * D2.dx
    ix = np.digitize(X, edges) - 1; iy = np.digitize(Y, edges) - 1
    inside = (ix >= 0) & (ix < nb) & (iy >= 0) & (iy < nb)
    P = np.zeros((nb, nb)); np.add.at(P, (ix[inside], iy[inside]), rho[inside])
    P = np.append(P.ravel(), rho.sum() - P.sum())
    jx = np.digitize(pos[0], edges) - 1; jy = np.digitize(pos[1], edges) - 1
    ins = (jx >= 0) & (jx < nb) & (jy >= 0) & (jy < nb)
    Nh = np.zeros((nb, nb)); np.add.at(Nh, (jx[ins], jy[ins]), 1)
    Nh = np.append(Nh.ravel(), (~ins).sum())
    E = P / P.sum() * pos.shape[1]
    keep = E >= 5
    c2 = np.sum((Nh[keep] - E[keep])**2 / E[keep]) + (np.sum(Nh[~keep]) - np.sum(E[~keep]))**2 / max(np.sum(E[~keep]), 1e-9)
    dof = keep.sum()
    return c2, dof, chi2dist.sf(c2, dof)


eq = {}
for tchk in (0.0, 12.0, 24.0):
    i = int(round(tchk / DT))
    c2, dof, pv = chi2_test(D2.evolve(pkm, tchk), traj2[i])
    eq[tchk] = (c2, dof, pv)
    check(f"S2c 2D equivariance chi2 t={tchk}", pv > TOL["chi2_p"], f"chi2 = {c2:.1f} / {dof} dof, p = {pv:.3f}")
res["2d"]["equivariance"] = {str(k): v for k, v in eq.items()}
# ensemble mean vs field mean (ZBW circle of the mean)
mx, my = mean_tr["mixed_k0"]
idx = np.round(ts2 / DT).astype(int)
ex, ey = traj2[idx, 0].mean(axis=1), traj2[idx, 1].mean(axis=1)
se2 = np.std(traj2[-1, 0]) / np.sqrt(N2)
dev2 = max(np.max(np.abs(ex - mx)), np.max(np.abs(ey - my)))
check("S2d 2D ensemble mean = field <r>(t)", dev2 < TOL["mean_sigma"] * se2, f"max dev {dev2:.4f}, SE {se2:.4f}")
# central surfers: loop radius and speed
rr0 = np.hypot(pos0[0], pos0[1])
cen = np.argsort(rr0)[:20]
sp = np.sqrt(np.sum(np.diff(traj2[:, :, cen], axis=0)**2, axis=1)) / DT
loops = []
for k, i in enumerate(cen):
    xs, ys = traj2[:, 0, i], traj2[:, 1, i]
    tt = times2 <= np.pi * 2
    Ax_, wx_, _ = zbw_fit(times2, xs, 2 * np.pi)
    loops.append((Ax_, wx_, float(np.mean(sp[times2[1:] <= 2 * np.pi, k]))))
loops = np.array(loops)
res["2d"]["central_surfers"] = dict(loop_radius_mean=float(loops[:, 0].mean()), omega_mean=float(loops[:, 1].mean()),
                                    speed_mean_first_2periods=float(loops[:, 2].mean()))
print(f"20 central surfers, first 2 ZBW periods: loop radius {loops[:, 0].mean():.4f} lambda_C "
      f"({loops[:, 0].mean() * LAMBDA_C:.3e} m), omega {loops[:, 1].mean():.4f}, mean speed {loops[:, 2].mean():.4f} c", flush=True)

# drift 2D: positive and mixed with k0 = (0.5, 0)
psi = gaussian_packet(D2, sig2, [0.5, 0], [1, 1])
pkk = D2.fft(psi)
Hp = D2.apply_H(pkk); pplus = 0.5 * (pkk + Hp / D2.E); pminus = pkk - pplus
vpred_mix = (np.sum(np.abs(pplus)**2 * D2.hx / D2.E) - np.sum(np.abs(pminus)**2 * D2.hx / D2.E)) / np.sum(np.abs(pkk)**2)
psi_p = D2.ifft(pplus); psi_p /= np.sqrt(D2.norm(psi_p))
ppk = D2.fft(psi_p)
vpred_pos = np.sum(np.abs(ppk)**2 * D2.hx / D2.E) / np.sum(np.abs(ppk)**2)
traj_k = {}
for nm, ps, pred in (("positive", psi_p, vpred_pos), ("mixed", psi, vpred_mix)):
    p0 = D2.birth(ps, 1000, rng)
    T = 20.0
    tt, tr, vm = D2.run_surfers(ps, p0, T, record_every=5)
    traj_k[nm] = (tt, tr)
    d = (tr[-1, 0] - tr[0, 0]) / T
    tol = TOL["drift_rel"] * abs(pred) + 4 * d.std() / np.sqrt(1000) + (0.5 / T if nm == "mixed" else 0)
    res["2d"][f"drift_{nm}"] = dict(drift=d.mean(), pred=pred, vmax=vm)
    check(f"S3 2D drift k0=0.5 {nm}", abs(d.mean() - pred) < tol and vm <= TOL["vmax"],
          f"surfer drift {d.mean():.5f} c, predicted {pred:.5f} c (v_g(0.5) = {0.5 / np.sqrt(1.25):.5f})")

# G2: rigid rotation of surfers in the spin-up packet (sigma = 4)
ps = D2.project(gaussian_packet(D2, sig2, [0, 0], [1, 0]), +1); ps /= np.sqrt(D2.norm(ps))
p0 = D2.birth(ps, 1000, rng)
tt, tr, vm = D2.run_surfers(ps, p0, 20.0, record_every=10)
ang = np.unwrap(np.arctan2(tr[:, 1], tr[:, 0]), axis=0)
rr = np.hypot(tr[0, 0], tr[0, 1])
sel = (rr > 0.5 * sig2) & (rr < 1.5 * sig2)
Om = np.median((ang[-1, sel] - ang[0, sel]) / tt[-1])
Om_pred = 1 / (2 * sig2**2)
res["2d"]["spin_rotation"] = dict(Omega=Om, Omega_pred_pauli=Om_pred, vmax=vm,
                                  vphi_at_sigma=Om * sig2, vphi_at_sigma_SI=Om * sig2 * 299792458.0)
check("G2 surfers rotate with Omega = 1/(2 m sigma^2)", abs(Om / Om_pred - 1) < TOL["omega_rel"],
      f"Omega = {Om:.5f} (pred {Om_pred:.5f}) mc^2/hbar; speed at r = sigma: {Om * sig2:.4f} c")
# post-hoc (added after the first run, NOT a test): the frozen prediction forgot that the packet spreads,
# sigma(t)^2 = sigma^2 (1 + (t / 2 m sigma^2)^2), and Omega(t) = 1 / (2 m sigma(t)^2):
Om_spread = np.arctan(tt[-1] / (2 * sig2**2)) / tt[-1]
res["2d"]["spin_rotation"]["Omega_pred_with_spreading_posthoc"] = Om_spread
print(f"   post-hoc: average Omega with spreading = {Om_spread:.5f}; measured / that = {Om / Om_spread:.4f}", flush=True)

# figure 2: Zitterbewegung loop (2D) + surfers
fig, ax = plt.subplots(1, 3, figsize=(13, 4.2), dpi=85)
mx, my = mean_tr["mixed_k0"]
ax[0].plot(mx[ts2 <= 8], my[ts2 <= 8], lw=1.2)
ax[0].set_aspect("equal"); ax[0].set_title("<r>(t), pacchetto misto fermo\n(cerchio di raggio ~hbar/2mc)")
ax[0].set_xlabel("x [hbar/mc]"); ax[0].set_ylabel("y [hbar/mc]")
for k, i in enumerate(cen[:6]):
    s = times2 <= 8
    ax[1].plot(traj2[s, 0, i], traj2[s, 1, i], lw=0.8)
ax[1].set_aspect("equal"); ax[1].set_title("6 surfisti centrali, t < 8 hbar/mc^2\n(anelli a ~c, raggio ~0.5 hbar/mc)")
ax[1].set_xlabel("x [hbar/mc]")
tt, tr = traj_k["mixed"]
r0 = np.hypot(tr[0, 0], tr[0, 1])
for i in np.argsort(r0)[:6]:
    ax[2].plot(tr[:, 0, i], tr[:, 1, i], lw=0.8)
ax[2].set_title("pacchetto misto con k0 = 0.5 mc:\nanelli + deriva (eliche in x,y,t)")
ax[2].set_xlabel("x [hbar/mc]")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "dirac_zbw_loop.png")); plt.close(fig)

res["verdict"] = verdict
res["runtime_s"] = time.time() - t_start
with open(os.path.join(OUT, "dirac_surf2d.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
print(f"=== {sum(v for _, v in verdict)}/{len(verdict)} PASS, runtime {time.time() - t_start:.1f} s ===", flush=True)
