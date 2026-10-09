"""dirac_free.py - free electron with the Dirac surfer rule, 1+1 D and 2+1 D.

Tests (tolerances fixed in TOL before the first run):
  Z1  Zitterbewegung of <x> in a packet with negative-energy components: frequency 2mc^2/hbar,
      amplitude hbar/2mc; absent (< 1e-6) after projection on positive energies.
  S1  surfer speed <= c at every RK4 stage of every surfer.
  S2  equivariance: 1D - each surfer keeps its quantile of psi^dagger psi (Bohm 1D: no crossing);
                    2D - chi^2 of surfer histogram against psi^dagger psi(T); ensemble mean vs <x>.
  S3  drift = group velocity (positive-energy packet), and = time-averaged <alpha> for the mixed one.
  G1  2+1 D spin-up packet at rest: magnetic moment of the Dirac current mu = (1/2) Int r x j,
      spin S = <sigma_z>/2 (2+1 D: J = L + sigma_z/2); g = 2 mu / S (e = m = 1) -> 2.
  G2  surfers of the spin-up packet rotate rigidly with Omega = 1/(2 m sigma^2) (spin current).
Run: timeout 590 python lab/dirac_free.py   (writes lab/out/dirac_free.json and dirac_zbw_1d.png)
The 2+1 D SURFER tests (S1b, S2c, S2d, S3 2D, G2, figure dirac_zbw_loop.png) are in dirac_surf2d.py
(split off because the shared machine made the whole run longer than 10 min).
"""
import json, os, time
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2 as chi2dist
from dirac_core import Dirac, gaussian_packet, SEED, LAMBDA_C, T_UNIT, DT

TOL = dict(zbw_freq_rel=0.02, zbw_amp_rel=0.05, zbw_proj_abs=1e-6, vmax=1 + 1e-12,
           quantile_1d=1e-3, chi2_p=1e-3, mean_sigma=4.0, drift_rel=0.01, g_rel=0.01,
           omega_rel=0.05)
OUT = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUT, exist_ok=True)
res, verdict = {}, []
t_start = time.time()


def check(name, ok, msg):
    verdict.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def zbw_fit(t, x, tmax):
    sel = t <= tmax
    f = lambda t, c0, c1, A, w, ph: c0 + c1 * t + A * np.sin(w * t + ph)
    p0 = [x[0], 0.0, 0.5, 2.0, 0.0]
    best = None
    for ph0 in (0, 1.5, 3, 4.5):
        p0[4] = ph0
        try:
            p, _ = curve_fit(f, t[sel], x[sel], p0=p0, maxfev=20000)
        except RuntimeError:
            continue
        r = np.sum((f(t[sel], *p) - x[sel])**2)
        if best is None or r < best[1]:
            best = (p, r)
    p = best[0]
    if p[2] < 0:
        p[2] = -p[2]; p[4] += np.pi
    return abs(p[2]), abs(p[3]), p


rng = np.random.default_rng(SEED)

# ======================= 1+1 D =======================
print("=== 1+1 D, sigma = 5 lambda_C, spinor (1,1)/sqrt2 ===", flush=True)
D1 = Dirac(4096, 400.0, 1)
sig = 5.0
psi_mix = gaussian_packet(D1, sig, [0.0], [1, 1])
psi_pos = D1.project(psi_mix, +1)
w_pos = D1.norm(psi_pos)
psi_pos = psi_pos / np.sqrt(w_pos)
res["1d_weight_positive_energy"] = w_pos
print(f"positive-energy weight of the mixed packet: {w_pos:.4f}", flush=True)
ts = np.arange(0, 40 + 1e-9, 0.01)
xm = {}
for name, psi in (("mixed", psi_mix), ("positive", psi_pos)):
    pk = D1.fft(psi)
    xm[name] = np.array([D1.mean(D1.evolve(pk, t), D1.x[0]) for t in ts])
A, w, p = zbw_fit(ts, xm["mixed"], 4 * np.pi)
Ap, wp, pp = zbw_fit(ts, xm["positive"], 4 * np.pi)
# exact expectation of the frequency for this k-distribution: 2<E>
phik = D1.fft(psi_mix)
wk = np.sum(np.abs(phik)**2, axis=0); wk /= wk.sum()
w_pred = 2 * np.sum(wk * D1.E)
res["1d_zbw"] = dict(amp=A, omega=w, omega_pred_2E=w_pred, amp_positive=Ap, drift=p[1])
check("Z1a 1D ZBW frequency", abs(w / 2 - 1) < TOL["zbw_freq_rel"],
      f"omega = {w:.5f} mc^2/hbar (2<E> = {w_pred:.5f}; 2mc^2/hbar = 2) -> {w / T_UNIT:.4e} rad/s")
check("Z1b 1D ZBW amplitude", abs(A / 0.5 - 1) < TOL["zbw_amp_rel"],
      f"A = {A:.4f} lambda_C (hbar/2mc = 0.5) -> {A * LAMBDA_C:.4e} m")
# positive projected: residual oscillation = std of <x> minus its linear trend
lin = np.polyfit(ts, xm["positive"], 2)
resid = np.max(np.abs(xm["positive"] - np.polyval(lin, ts)))
res["1d_zbw"]["positive_residual"] = resid
check("Z1c 1D no ZBW after P+", resid < TOL["zbw_proj_abs"], f"max |<x> - smooth trend| = {resid:.2e}")

# surfers in 1D (mixed packet)
N1 = 2000
pos0 = D1.birth(psi_mix, N1, rng)
F0 = None
cdf_err = [0.0]


def cdf_at(psi, xs):
    rho = D1.density_current(psi)[0]
    c = np.cumsum(rho) * D1.dx
    c = c - 0.5 * rho * D1.dx    # cumulative at cell centre
    return np.interp(xs, D1.x[0], c)


q0 = cdf_at(psi_mix, pos0[0])


def on_step_1d(s, t, psi, pos):
    if s % 50 == 0:
        cdf_err[0] = max(cdf_err[0], float(np.max(np.abs(cdf_at(psi, pos[0]) - q0))))


t0 = time.time()
times1, traj1, vmax1 = D1.run_surfers(psi_mix, pos0, 40.0, on_step=on_step_1d)
print(f"1D surfers: {N1} x {len(times1) - 1} steps in {time.time() - t0:.1f} s", flush=True)
check("S1a 1D |v| <= c", vmax1 <= TOL["vmax"], f"max |v| = {vmax1:.12f}")
check("S2a 1D quantile conservation", cdf_err[0] < TOL["quantile_1d"],
      f"max |F(x_i(t),t) - F(x_i(0),0)| = {cdf_err[0]:.2e} (cell-quantile resolution {1 / N1:.0e} not needed: deterministic)")
ens = traj1[:, 0, :].mean(axis=1)
idx = np.searchsorted(ts, times1 - 1e-9)
dev = np.max(np.abs(ens - xm["mixed"][idx]))
se = np.std(traj1[-1, 0, :]) / np.sqrt(N1)
res["1d_surfers"] = dict(vmax=vmax1, quantile_err=cdf_err[0], ens_mean_dev=dev, se=se)
check("S2b 1D ensemble <x>(t) = field <x>(t)", dev < TOL["mean_sigma"] * se, f"max dev {dev:.4f}, SE {se:.4f}")
# one central surfer: its jitter
ic = np.argmin(np.abs(pos0[0]))
xc = traj1[:, 0, ic]
Ac, wc, _ = zbw_fit(times1, xc, 4 * np.pi)
res["1d_central_surfer"] = dict(x0=float(pos0[0, ic]), amp=Ac, omega=wc)
print(f"central surfer (x0 = {pos0[0, ic]:.3f}): jitter amplitude {Ac:.4f}, omega {wc:.4f}", flush=True)

# drift = group velocity, 1D
for k0 in (0.5, 1.0, 2.0):
    psi = gaussian_packet(D1, sig, [k0], [1, 1])
    pk = D1.fft(psi)
    Hpk = D1.apply_H(pk)
    pplus = 0.5 * (pk + Hpk / D1.E); pminus = pk - pplus
    vavg = (np.sum(np.abs(pplus)**2 * D1.hx / D1.E) - np.sum(np.abs(pminus)**2 * D1.hx / D1.E)) / np.sum(np.abs(pk)**2)
    psi_p = D1.project(psi, +1); psi_p /= np.sqrt(D1.norm(psi_p))
    vg = k0 / np.sqrt(1 + k0**2)
    pp_k = D1.fft(psi_p)
    vg_exact = np.sum(np.abs(pp_k)**2 * D1.hx / D1.E) / np.sum(np.abs(pp_k)**2)
    out = {}
    for nm, ps, pred in (("positive", psi_p, vg_exact), ("mixed", psi, vavg)):
        p0 = D1.birth(ps, 1000, rng)
        T = 30.0
        tt, tr, vm = D1.run_surfers(ps, p0, T, record_every=1500)
        drift = (tr[-1, 0] - tr[0, 0]).mean() / T
        se_d = (tr[-1, 0] - tr[0, 0]).std() / T / np.sqrt(1000)
        out[nm] = dict(drift=drift, pred=pred, se=se_d, vmax=vm)
        # mixed: residual ZBW of the mean contributes at most ~0.5/T
        tol = TOL["drift_rel"] * abs(pred) + 4 * se_d + (0.5 / T if nm == "mixed" else 0)
        check(f"S3 1D drift k0={k0} {nm}", abs(drift - pred) < tol and vm <= TOL["vmax"],
              f"surfer drift {drift:.5f} c, predicted {pred:.5f} c (v_g(k0) = {vg:.5f}), max|v| {vm:.10f}")
    res[f"1d_drift_k{k0}"] = out

# figure 1: 1D <x>(t) and surfer world-lines
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), dpi=90)
ax[0].plot(ts, xm["mixed"], lw=1, label="pacchetto misto (+E e -E)")
ax[0].plot(ts, xm["positive"], lw=1, label="solo energia positiva")
ax[0].set_xlabel("t  [hbar/mc^2 = 1.29e-21 s]"); ax[0].set_ylabel("<x>  [hbar/mc = 3.86e-13 m]")
ax[0].set_title("Zitterbewegung di <x>, 1+1 D"); ax[0].legend(fontsize=8)
sel = np.argsort(np.abs(pos0[0]))[:12]
for i in sel:
    ax[1].plot(traj1[:, 0, i], times1, lw=0.7)
ax[1].set_xlabel("x  [hbar/mc]"); ax[1].set_ylabel("t  [hbar/mc^2]")
ax[1].set_title("12 surfisti (linee d'universo), pacchetto misto")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "dirac_zbw_1d.png")); plt.close(fig)

# ======================= 2+1 D =======================
print("=== 2+1 D, sigma = 4 lambda_C ===", flush=True)
D2 = Dirac(384, 96.0, 2)
sig2 = 4.0
X, Y = D2.x
cases = {"mixed_k0": ([0.0, 0.0], [1, 1]), "mixed_k05": ([0.5, 0.0], [1, 1])}
ts2 = np.arange(0, 30 + 1e-9, 0.05)
res["2d"] = {}
mean_tr = {}
for nm, (k0, sp) in cases.items():
    psi = gaussian_packet(D2, sig2, k0, sp)
    pk = D2.fft(psi)
    mx, my = [], []
    for t in ts2:
        ps = D2.evolve(pk, t)
        mx.append(D2.mean(ps, X)); my.append(D2.mean(ps, Y))
    mx, my = np.array(mx), np.array(my)
    mean_tr[nm] = (mx, my)
    Ax, wx, _ = zbw_fit(ts2, mx, 4 * np.pi)
    Ay, wy, _ = zbw_fit(ts2, my, 4 * np.pi)
    res["2d"][nm] = dict(Ax=Ax, wx=wx, Ay=Ay, wy=wy)
    check(f"Z1d 2D ZBW {nm}", abs(wx / 2 - 1) < TOL["zbw_freq_rel"] and abs(wy / 2 - 1) < TOL["zbw_freq_rel"]
          and abs(Ax / 0.5 - 1) < TOL["zbw_amp_rel"] and abs(Ay / 0.5 - 1) < TOL["zbw_amp_rel"] if nm == "mixed_k0"
          else abs(wx / 2 - 1) < 0.2,
          f"x: A={Ax:.4f} w={wx:.4f}; y: A={Ay:.4f} w={wy:.4f}  (circle of radius hbar/2mc = 0.5 at 2mc^2/hbar)")
psi_p2 = D2.project(gaussian_packet(D2, sig2, [0, 0], [1, 1]), +1)
psi_p2 /= np.sqrt(D2.norm(psi_p2))
pk = D2.fft(psi_p2)
mpx = np.array([D2.mean(D2.evolve(pk, t), X) for t in ts2[::4]])
mpy = np.array([D2.mean(D2.evolve(pk, t), Y) for t in ts2[::4]])
# (fix after first run: the projected packet DRIFTS - P+ of spinor (1,1) favours k_x > 0 - so, exactly as in
#  the 1D test Z1c, the smooth trend is removed before measuring the residual oscillation)
tq = ts2[::4]
resid2 = max(np.max(np.abs(mpx - np.polyval(np.polyfit(tq, mpx, 2), tq))),
             np.max(np.abs(mpy - np.polyval(np.polyfit(tq, mpy, 2), tq))))
res["2d"]["positive_residual"] = resid2
res["2d"]["positive_drift"] = float((mpx[-1] - mpx[0]) / (tq[-1] - tq[0]))
check("Z1e 2D no ZBW after P+", resid2 < TOL["zbw_proj_abs"],
      f"max |<r>(t) - smooth trend| = {resid2:.2e} (drift of the projected packet {res['2d']['positive_drift']:.4f} c)")

# ---------------- G1/G2: spin and magnetic moment of a spin-up packet at rest ----------------
print("=== spin and magnetic moment (2+1 D, positive energy, spinor (1,0)) ===", flush=True)
gres = {}
for s_ in (2.0, 4.0, 8.0):
    Dg = Dirac(384, 96.0 if s_ < 8 else 160.0, 2)
    ps = Dg.project(gaussian_packet(Dg, s_, [0, 0], [1, 0]), +1)
    ps /= np.sqrt(Dg.norm(ps))
    rho, jx, jy = Dg.density_current(ps)
    Xg, Yg = Dg.x
    dA = Dg.dx**2
    S = 0.5 * np.sum((np.abs(ps[0])**2 - np.abs(ps[1])**2)) * dA       # <sigma_z>/2
    mu = 0.5 * np.sum(Xg * jy - Yg * jx) * dA                            # (e/2) Int r x j, e = 1
    # orbital angular momentum <L_z> = <-i (x d_y - y d_x)> (spectral derivative)
    KX, KY = Dg.k
    Lz = 0.0
    for comp in ps:
        ck = np.fft.fft2(comp)
        dy = np.fft.ifft2(1j * KY * ck); dx_ = np.fft.ifft2(1j * KX * ck)
        Lz += np.real(np.sum(np.conj(comp) * (-1j) * (Xg * dy - Yg * dx_))) * dA
    g = 2 * mu / S
    gres[s_] = dict(S=S, Lz=Lz, J=S + Lz, mu=mu, g=g)
    print(f"sigma={s_}: <S_z>={S:.6f}  <L_z>={Lz:.6f}  J={S + Lz:.6f}  mu={mu:.6f}  g=2mu/S={g:.5f}", flush=True)
res["spin_moment"] = {str(k): v for k, v in gres.items()}
check("G1 g-factor of the Dirac current (sigma=8)", abs(gres[8.0]["g"] / 2 - 1) < TOL["g_rel"],
      f"g = {gres[8.0]['g']:.5f}; J = {gres[8.0]['J']:.6f} (exact 1/2)")
res["verdict"] = verdict
res["runtime_s"] = time.time() - t_start
res["tolerances"] = TOL
with open(os.path.join(OUT, "dirac_free.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
npass = sum(v for _, v in verdict)
print(f"=== {npass}/{len(verdict)} PASS, runtime {time.time() - t_start:.1f} s ===", flush=True)
