"""zz_atom3d.py - T3 (3-D check): the FULL 3-D zig-zag walk (cubic lattice, 4-component walkers) binds an
electron to a point nucleus through the Coulomb phase alone.

Physical hydrogen (Z alpha = 1/137) does not fit a 3-D lattice that also resolves hbar/mc (see
zz_hydrogen.py). The Dirac-Coulomb problem depends on Z alpha only, so here the same rule runs for a
hydrogen-like NEODYMIUM ion, Z = 60 (Z alpha = 0.4378, Bohr radius lambda_C/(Z alpha) = 2.28 lambda_C),
on lattices a = 0.4, 0.3, 0.24 lambda_C in a periodic box of side 24 lambda_C (nucleus between 8 sites).
Inputs: m_e (theta = a/lambda_C), alpha, Z = 60. Nothing else.

How the ground state is found (the "listening" rule, as in zz_listen.py): start from an electron blob
(positive-energy part of a spin-up exponential of the Bohr size, projected with the free walk's own
eigenvectors), run the walk, and add up the balls' amplitudes with a lock-in at the blob's mean frequency
under a Gaussian time window (std 27 hbar/mc^2): states away from that frequency average out. The energy is
the frequency of the filtered state (Rayleigh quotient of two ticks).

Tests (TOL fixed before the first run):
  A1 1s binding energy extrapolated to a -> 0 (fit in a^2) vs Sommerfeld-Dirac formula: < 1 %.
  A2 the surfer of the 3-D state circulates about the spin axis: density-weighted v_phi / (Z alpha c sin th)
     = 1 within 5 % (Dirac: |v| = c sqrt(1 - gamma^2) sin th = Z alpha c sin th exactly); radial and
     z components < 5 % of Z alpha c.
Run: timeout 590 python lab/zz_atom3d.py
"""
import json, os, time
import numpy as np
from zz_core import bloch_batch, dirac_formula, dirac_velocity, ALPHA, MC2_EV, SEED
from zz_fast import FastDirac3D

TOL = dict(E=0.01, v=0.05)
OUT = os.path.join(os.path.dirname(__file__), "out")
res, verdict = {}, []
t_start = time.time()
Z = 60
ZA = Z * ALPHA
HALF = 12.0
S_T = 27.0                      # lock-in window std (time, hbar/mc^2)


def check(name, ok, msg):
    verdict.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def run(a):
    n = int(round(2 * HALF / a))
    x = (np.arange(n) - n / 2 + 0.5) * a
    X, Y, Zc = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Zc**2)
    V = -ZA / r
    W = FastDirac3D(n, a, a, V)
    # initial electron blob: exp(-r Z alpha) * (1,0,1,0)/sqrt2, positive-energy part (free walk eigenvectors)
    psi = np.exp(-r * ZA)[None] * (np.array([1, 0, 1, 0], complex) / np.sqrt(2))[:, None, None, None]
    pk = np.fft.fftn(psi, axes=(1, 2, 3)).reshape(4, -1)
    kk = 2 * np.pi * np.fft.fftfreq(n)
    K = np.stack([q.ravel() for q in np.meshgrid(kk, kk, kk, indexing="ij")], -1)
    out = np.zeros_like(pk)
    for b0 in range(0, K.shape[0], 200000):
        sl = slice(b0, b0 + 200000)
        w, v = np.linalg.eig(bloch_batch(K[sl], a, 0.5))
        om = -np.angle(w)
        o = np.argsort(-om, axis=1)[:, :2]                      # the two positive frequencies (spin)
        Q, _ = np.linalg.qr(np.take_along_axis(v, o[:, None, :], axis=2))   # orthonormal basis of the pair
        c = np.einsum("nab,na->nb", Q.conj(), pk[:, sl].T)
        out[:, sl] = np.einsum("nab,nb->na", Q, c).T
    psi = np.fft.ifftn(out.reshape(4, n, n, n), axes=(1, 2, 3))
    psi /= np.sqrt(np.sum(np.abs(psi)**2))

    def rayleigh(p):
        q = W.tick(W.tick(p.copy()))
        lam = np.vdot(p, q) / np.vdot(p, p)
        return -np.angle(lam) / (2 * a)
    W.tick_count = 0
    om0 = rayleigh(psi)
    W.tick_count = 0
    nt = int(round(4.5 * S_T / a))
    if nt % 2:
        nt += 1
    c0, s = nt / 2, S_T / a
    acc = np.zeros_like(psi)
    t0 = time.time()
    for t in range(nt + 1):
        acc += np.exp(-0.5 * ((t - c0) / s)**2) * np.exp(1j * om0 * a * t) * psi
        W.tick(psi)
        if t == 5:
            print(f"  a = {a}: n = {n}, {nt} ticks, ~{(time.time() - t0) / 6 * nt:.0f} s", flush=True)
    acc /= np.sqrt(np.sum(np.abs(acc)**2))
    W.tick_count = 0
    om1 = rayleigh(acc)
    # surfer velocity field of the filtered state
    v, rho = dirac_velocity(acc)
    rc = np.sqrt(X**2 + Y**2)
    rc[rc == 0] = 1e-12
    sinth = rc / r
    vphi = (-Y * v[0] + X * v[1]) / rc
    vrad = (X * v[0] + Y * v[1] + Zc * v[2]) / r
    wgt = rho * (r < 3 / ZA)
    vphi_ratio = float(np.sum(wgt * vphi / (ZA * np.maximum(sinth, 0.2)) * (sinth > 0.2)) / np.sum(wgt * (sinth > 0.2)))
    return dict(a=a, n=n, ticks=nt, omega_start=om0, omega=om1, binding=1 - om1, binding_start=1 - om0,
                vphi_over_Zac=vphi_ratio, vrad_rms_over_Zac=float(np.sqrt(np.sum(wgt * vrad**2) / wgt.sum()) / ZA),
                vz_rms_over_Zac=float(np.sqrt(np.sum(wgt * v[2]**2) / wgt.sum()) / ZA), seconds=time.time() - t0), acc


E_dirac = 1 - dirac_formula(1, -1, ZA)
print(f"Z = {Z}, Z alpha = {ZA:.5f}; Dirac 1s binding (check only) = {E_dirac:.6f} mc^2 = {E_dirac * MC2_EV / 1e3:.3f} keV", flush=True)
res["runs"] = []
last = None
for a in (0.4, 0.3, 0.24):
    r_, acc = run(a)
    res["runs"].append(r_)
    last = acc
    print(f"  a = {a}: binding {r_['binding']:.6f} mc^2 = {r_['binding'] * MC2_EV / 1e3:.3f} keV (start blob {r_['binding_start']:.5f});"
          f" v_phi/(Z alpha c sin) = {r_['vphi_over_Zac']:.4f}, v_r rms {r_['vrad_rms_over_Zac']:.3f}, v_z rms {r_['vz_rms_over_Zac']:.3f}"
          f"  ({r_['seconds']:.0f} s)", flush=True)
aa = np.array([q["a"] for q in res["runs"]])
bb = np.array([q["binding"] for q in res["runs"]])
fit = np.polyfit(aa**2, bb, 2)
fit1 = np.polyfit(aa[1:]**2, bb[1:], 1)
ex, ex1 = fit[-1], fit1[-1]
res["extrapolated"] = dict(binding_quadratic=ex, binding_linear_two_finest=ex1, dirac=E_dirac)
print(f"  a -> 0: {ex:.6f} mc^2 (quadratic in a^2), {ex1:.6f} (linear, two finest); Dirac {E_dirac:.6f}", flush=True)
check("A1 3-D walk 1s binding (Z = 60) vs Dirac", abs(ex / E_dirac - 1) < TOL["E"], f"{100 * (ex / E_dirac - 1):+.3f} % (linear: {100 * (ex1 / E_dirac - 1):+.3f} %)")
rl = res["runs"][-1]
check("A2 3-D surfer circulates at Z alpha c sin(theta)", abs(rl["vphi_over_Zac"] - 1) < TOL["v"] and rl["vrad_rms_over_Zac"] < TOL["v"] and rl["vz_rms_over_Zac"] < TOL["v"],
      f"v_phi/(Z alpha c sin th) = {rl['vphi_over_Zac']:.4f}, v_r rms/(Z alpha c) = {rl['vrad_rms_over_Zac']:.4f}, v_z {rl['vz_rms_over_Zac']:.4f} (a = 0.24)")
np.save(os.path.join(OUT, "zz_atom3d_density.npy"), np.sum(np.abs(last)**2, axis=0)[:, :, last.shape[3] // 2].astype(np.float32))
res["verdict"] = verdict
res["seconds"] = time.time() - t_start
json.dump(res, open(os.path.join(OUT, "zz_atom3d.json"), "w"), indent=1, default=float)
npass = sum(v for _, v in verdict)
print(f"=== {npass} PASS / {len(verdict) - npass} FAIL  ({time.time() - t_start:.0f} s)", flush=True)
