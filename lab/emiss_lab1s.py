"""emiss_lab1s.py - why the Laboratory atom keeps oscillating between "1s" and "2p" after the light has gone
(data of the main thread: /mnt/project-files/ripresa/dati/atomo_accoppiato/serie.csv, run of accoppiato.rs).

1. Replays EXACTLY the Laboratory's listening of its 1s (zz.rs Atomo::ascolta(0.885, 0.915, 600): Hann window over
   600 ticks, every tick, candidate frequencies every 0.0025, the candidate nearest the fine spectral peak) from the
   same start (the radial 1s with a = 0.5).
2. Listens the clean levels of the 3-D rule near 1s (0.8996, 0.8863, 0.8404) with a long Gaussian window and
   decomposes the Laboratory "1s" on them.
3. Predicts, with no free parameter, the beat of p_exc = 1 - |<1s_lab|psi(t)>|^2 for the atom left alone
   (period 2 pi / (w_a - w_b), depth from the weights) and checks it against a free run of the rule and against
   the Laboratory series.
4. Describes the extra levels: radial profile (where they live) and their weight on the four "colours"
   (sub-lattices (i+j+k) mod 2 and the staggered pattern) - is it the first cell (NUCLEO_1S) or the doubling?
Run: python lab/emiss_lab1s.py   ->  lab/out/emiss_lab1s.json
"""
import json, os, time
import numpy as np
from emiss_core import *
from emiss_levels import listen, rayleigh

t0 = time.time()
n, a = 48, 0.5
Za = 60 * ALPHA
D = Dirac3D(n, a)
res = {}

# 1. the Laboratory's 1s
om_r, rg, g, f = radial_state(-1, 1, 0.5, Za, rmax=40.0)
psi0 = dirac_3d_from_radial(rg, g, f, -1, 0.5, D.x)
lo, hi, ticks = 0.885, 0.915, 600
cand = lo + 0.0025 * np.arange(int(round((hi - lo) / 0.0025)) + 1)
acc = np.zeros((len(cand),) + psi0.shape, complex)
ser = []
p = psi0.copy()
for t in range(ticks):
    w = np.sin(np.pi * (t + 0.5) / ticks)**2
    tm = t * a
    ser.append((tm, ovl(psi0, p)))
    for k, om in enumerate(cand):
        acc[k] += (w * np.exp(1j * om * tm)) * p
    D.tick(p, t)
tt = np.array([s[0] for s in ser]); cc = np.array([s[1] for s in ser])
fine = np.arange(lo, hi + 1e-9, 0.0002)
h = np.sin(np.pi * tt / (ticks * a))**2
S = np.abs(np.exp(1j * np.outer(fine, tt)) @ (h * cc))**2
kb = np.argmax(S[1:-1]) + 1
om_lab = fine[kb] + 0.0002 * 0.5 * (S[kb - 1] - S[kb + 1]) / (S[kb - 1] - 2 * S[kb] + S[kb + 1])
kc = np.argmin(np.abs(cand - om_lab))
lab = acc[kc] / np.sqrt(norm2(acc[kc]))
res["lab_1s"] = dict(omega_peak=om_lab, candidate=cand[kc], rayleigh=rayleigh(D, lab)[0])
print(f"Laboratory 1s: fine peak {om_lab:.5f} (manifest 0.89896), candidate used {cand[kc]:.4f}  ({time.time() - t0:.0f} s)", flush=True)

# 2. clean levels (start: the Laboratory state itself, it contains all three)
levels = {}
for name, w0 in (("a", 0.8996), ("b", 0.8863), ("c", 0.8404)):
    f_ = os.path.join(CK, f"lab1s_level_{name}.npy")
    if os.path.exists(f_):
        phi = np.load(f_)
    else:
        phi = listen(D, psi0, w0, 1200.0)
        np.save(f_, phi)
    om, rs = rayleigh(D, phi)
    levels[name] = (phi, om)
    print(f"level {name}: omega {om:.6f}, residual {rs:.1e}  ({time.time() - t0:.0f} s)", flush=True)
W = {k: abs(ovl(v[0], lab))**2 for k, v in levels.items()}
Wstart = {k: abs(ovl(v[0], psi0))**2 for k, v in levels.items()}
res["weights_lab"] = W
res["weights_radial_start"] = Wstart
res["omegas"] = {k: v[1] for k, v in levels.items()}
print("weights of the Laboratory 1s:", {k: round(v, 4) for k, v in W.items()}, " of the radial start:", {k: round(v, 4) for k, v in Wstart.items()})

# 3. the beat predicted (no parameter) vs a free run of the rule
T = 1600.0
tt = np.arange(0, T + 1e-9, 8.0)
amp = sum(W[k] * np.exp(-1j * levels[k][1] * tt) for k in W)
p_pred = 1 - np.abs(amp)**2
p = lab.copy()
p_run = []
for t in range(int(T / a) + 1):
    if t % 16 == 0:
        p_run.append(1 - abs(ovl(lab, p))**2)
    D.tick(p, t)
p_run = np.array(p_run)
period_pred = 2 * np.pi / (levels["a"][1] - levels["b"][1])
res["free"] = dict(t=tt.tolist(), p_pred=p_pred.tolist(), p_run=p_run.tolist(), period_pred=period_pred,
                   max_abs_diff=float(np.max(np.abs(p_pred - p_run[:len(tt)]))))
print(f"free atom: beat period predicted {period_pred:.1f}; p_exc predicted range {p_pred.min():.3f}-{p_pred.max():.3f}; "
      f"free run range {p_run.min():.3f}-{p_run.max():.3f}; max |pred - run| {res['free']['max_abs_diff']:.2e}", flush=True)
# the Laboratory series
d = np.genfromtxt("/mnt/project-files/ripresa/dati/atomo_accoppiato/serie.csv", delimiter=",", names=True)
m = d["t"] > 400
x = d["p_exc"][m] - d["p_exc"][m].mean()
om = np.linspace(0.002, 0.05, 4000)
Sx = np.abs(np.exp(1j * np.outer(om, d["t"][m])) @ (x * np.hanning(len(x))))
res["lab_series"] = dict(beat_omega=float(om[np.argmax(Sx)]), beat_period=float(2 * np.pi / om[np.argmax(Sx)]),
                         p_exc_t4=float(d["p_exc"][0]))
print(f"Laboratory series: p_exc beat at w = {om[np.argmax(Sx)]:.5f} (period {2 * np.pi / om[np.argmax(Sx)]:.0f}); "
      f"predicted w_a - w_b = {levels['a'][1] - levels['b'][1]:.5f}", flush=True)

# 4. where the extra levels live
X, Y, Zc = np.meshgrid(D.x, D.x, D.x, indexing="ij")
r = np.sqrt(X**2 + Y**2 + Zc**2)
I = np.arange(n)
par = (I[:, None, None] + I[None, :, None] + I[None, None, :]) % 2
desc = {}
for k, (phi, om) in levels.items():
    rho = np.sum(np.abs(phi)**2, 0)
    # staggered content: weight of the Fourier components near the corners (pi,...)
    F = np.fft.fftn(phi, axes=(1, 2, 3))
    kk = np.abs(np.fft.fftfreq(n))
    corner = np.zeros((2, 2, 2))
    for ix in (0, 1):
        for iy in (0, 1):
            for iz in (0, 1):
                mx = (kk >= 0.25) if ix else (kk < 0.25)
                my = (kk >= 0.25) if iy else (kk < 0.25)
                mz = (kk >= 0.25) if iz else (kk < 0.25)
                sel = mx[:, None, None] & my[None, :, None] & mz[None, None, :]
                corner[ix, iy, iz] = np.sum(np.abs(F[:, sel])**2)
    corner /= corner.sum()
    desc[k] = dict(omega=om, r_mean=float(np.sum(rho * r)), inside_1cell=float(np.sum(rho[r < 0.5])),
                   inside_r1=float(np.sum(rho[r < 1.0])),
                   k_corners={f"{ix}{iy}{iz}": float(corner[ix, iy, iz]) for ix in (0, 1) for iy in (0, 1) for iz in (0, 1)})
    print(f"level {k} ({om:.4f}): <r> = {desc[k]['r_mean']:.3f}, P(r<0.5) = {desc[k]['inside_1cell']:.3f}, "
          f"P(r<1) = {desc[k]['inside_r1']:.3f}, k-space corners {{000: {corner[0,0,0]:.3f}, 100+010+001: "
          f"{corner[1,0,0]+corner[0,1,0]+corner[0,0,1]:.3f}, 110+101+011: {corner[1,1,0]+corner[1,0,1]+corner[0,1,1]:.3f}, 111: {corner[1,1,1]:.3f}}}", flush=True)
res["levels_where"] = desc
json.dump(res, open(os.path.join(OUT, "emiss_lab1s.json"), "w"), indent=1, default=float)
print(f"done ({time.time() - t0:.0f} s)")
