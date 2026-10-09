"""zz_listen.py - T3 (dynamics): the hydrogen levels found by LISTENING to the balls (the ball rule only).

Same radial walk as zz_hydrogen.py (kappa = -1: s states), coarse lattice theta = 1.0 (so that the
binding phase per tick, theta x 2.7e-5, is not too small; the first run at theta = 0.5 with 3e6 ticks, and
a second at theta = 1.0 with 1.5e6 ticks, did not finish in 590 s on the shared machine), box 20 Bohr radii. The electron starts as a
Gaussian blob g = r exp(-r^2 / (2 (2 a0)^2)), f = 0: a mixture of 1s, 2s, 3s, ... and free states.
  Step 1 (spectroscopy): run 1 000 000 ticks and record how much of the starting blob is still there,
     C(n) = <chi0|chi_n> (demodulated by the rest frequency). Its spectrum shows a peak at every frequency
     where the atom rings.
  Step 2 (lock-in): run again and add up the balls' amplitudes with the phase e^{i omega n} of the 1s and 2s
     peaks under Gaussian time windows (three widths for 1s, one for 2s). The sum IS a stationary state; its frequency
     (Rayleigh quotient of one tick) is its energy.
  Reference: the same walk listened to for an infinite time = resolvent / shift-invert eigenvector
     (RadialWalk.eigen), as used for all numbers in zz_hydrogen.py.
Tests (TOL fixed before the first run):
  S1 the spectrum shows separate peaks for 1s and 2s whose positions agree with the resolvent within the
     spectral resolution 2 pi / N_ticks (per tick).
  S2 1s energy from the longest lock-in vs resolvent: relative difference of the binding < 1e-4; and the
     error falls as the listening window grows.
  S3 2s energy from the longest lock-in vs resolvent: < 1e-2 (2s and 3s are only ~ 1.5 resolutions apart).
Run: timeout 590 python lab/zz_listen.py
"""
import json, os, time
import numpy as np
from zz_core import RadialWalk, ALPHA, MC2_EV
from zz_fast import radial_listen

TOL = dict(peak_res=1.0, s1=1e-4, s2=1e-2)
OUT = os.path.join(os.path.dirname(__file__), "out")
res, verdict = {}, []
t_start = time.time()
A0 = 1 / ALPHA
TH = 1.0
N1 = 1_000_000
N2 = 1_000_000


def check(name, ok, msg):
    verdict.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


rw = RadialWalk(-1, TH, 20 * A0)
r = rw.r
g = r * np.exp(-r**2 / (2 * (2 * A0)**2))
f = 0 * g
chi0 = np.array([(g - 1j * f) / np.sqrt(2), (g + 1j * f) / np.sqrt(2)])
chi0 /= np.linalg.norm(chi0)
# sponge on the outer 20 % of the box: free (ionized) parts of the blob leave instead of bouncing back
# (run 2 showed equally spaced false peaks, 5.06 eV apart: the free part bouncing in the box)
xs = np.clip((r - 0.8 * r[-1]) / (0.2 * r[-1]), 0, 1)
DAMP = np.exp(-0.02 * xs**2)
print(f"radial walk theta = {TH}, N = {rw.N} sites, {N1} ticks for the spectrum", flush=True)
t0 = time.time()
C, _, _ = radial_listen(rw.M, rw.refl0, rw.refl1, chi0, N1, np.zeros(0), np.zeros(0), np.zeros(0), TH, True, DAMP)
print(f"  step 1 done in {time.time() - t0:.0f} s", flush=True)
# spectrum: positive frequency nu per tick = theta x binding (in m c^2)
npad = 1 << 24
win = np.hanning(N1)
S = np.abs(np.fft.fft(C * win, npad))**2
nu = 2 * np.pi * np.fft.fftfreq(npad)
band = (nu > 0.2e-5 * TH) & (nu < 1.5 * 0.5 * ALPHA**2 * TH * 1.3)
idx = np.nonzero(band)[0]
Sb = S[idx]
peaks = [idx[i] for i in range(1, len(idx) - 1) if Sb[i] > Sb[i - 1] and Sb[i] >= Sb[i + 1] and Sb[i] > 1e-3 * Sb.max()]


def refine(i):
    y0, y1, y2 = np.log(S[i - 1]), np.log(S[i]), np.log(S[i + 1])
    d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
    return nu[i] + d * (nu[1] - nu[0])


pk = sorted([refine(i) for i in peaks], reverse=True)
pk_eV = [p / TH * MC2_EV for p in pk]
print("  peaks (binding energies the atom rings at, eV): " + ", ".join(f"{e:.4f}" for e in pk_eV[:6]), flush=True)
res["peaks_eV"] = pk_eV[:8]
# reference: infinite listening (resolvent) at the same theta
ref = {}
for nm, k in (("1s", 0), ("2s", 1)):
    om, u, resid = rw.eigen(1 - pk[k] / TH)
    ref[nm] = 1 - om
print(f"  resolvent: 1s {ref['1s'] * MC2_EV:.6f} eV, 2s {ref['2s'] * MC2_EV:.6f} eV (this coarse lattice theta = 1.0)", flush=True)
reso = 2 * np.pi / N1 / TH * MC2_EV
d1 = abs(pk_eV[0] - ref["1s"] * MC2_EV)
d2 = abs(pk_eV[1] - ref["2s"] * MC2_EV)
res["resolvent_eV"] = {k: v * MC2_EV for k, v in ref.items()}
res["spectral_resolution_eV"] = reso
check("S1 spectrum peaks at 1s and 2s", d1 < TOL["peak_res"] * reso and d2 < TOL["peak_res"] * reso,
      f"|peak - resolvent| = {d1:.2e} / {d2:.2e} eV, resolution 2 pi/(N theta) = {reso:.2e} eV")

# step 2: lock-in at the two peaks, three window widths each
widths = np.array([N2 / 6, N2 / 12, N2 / 24])
nus = np.array([pk[0], pk[0], pk[0], pk[1]])
cent = np.full(4, N2 / 2)
wd = np.array([widths[0], widths[1], widths[2], widths[0]])
t0 = time.time()
_, acc, _ = radial_listen(rw.M, rw.refl0, rw.refl1, chi0, N2, nus, cent, wd, TH, False, DAMP)
print(f"  step 2 done in {time.time() - t0:.0f} s", flush=True)
Wm = rw.matrix()
res["lockin"] = {}
for k in range(4):
    nm = "1s" if k < 3 else "2s"
    v = acc[k].T.reshape(-1)
    v = v / np.linalg.norm(v)
    lam = np.vdot(v, Wm @ v)
    b = 1 - (-np.angle(lam) / TH)
    rel = abs(b / ref[nm] - 1)
    res["lockin"][f"{nm} window {wd[k]:.0f} ticks"] = dict(binding_eV=b * MC2_EV, rel_diff=rel)
    print(f"  lock-in {nm}, window std {wd[k]:.0f} ticks: binding {b * MC2_EV:.6f} eV, relative difference to resolvent {rel:.1e}", flush=True)
r1 = [res["lockin"][f"1s window {w:.0f} ticks"]["rel_diff"] for w in widths]
r2 = res["lockin"][f"2s window {widths[0]:.0f} ticks"]["rel_diff"]
check("S2 1s by listening = 1s by resolvent", r1[0] < TOL["s1"] and r1[0] <= r1[1] <= r1[2], f"{r1[0]:.1e} (shorter windows {r1[1]:.1e}, {r1[2]:.1e})")
check("S3 2s by listening", r2 < TOL["s2"], f"{r2:.1e}")
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 3.2))
    sel = (nu > 0) & (nu < 16 / MC2_EV * TH)
    ax.semilogy(nu[sel] / TH * MC2_EV, S[sel] / S[sel].max())
    for nm in ("1s", "2s"):
        ax.axvline(ref[nm] * MC2_EV, color="k", ls=":", lw=0.8)
    ax.set_xlabel("binding energy (eV) - walk at theta = 1.0"); ax.set_ylabel("spectrum of C(n)")
    ax.set_title("listening to the hydrogen walk"); ax.set_ylim(1e-8, 2)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "zz_listen.png"), dpi=90)
except Exception as ex:
    print("figure skipped:", ex)
res["verdict"] = verdict
res["seconds"] = time.time() - t_start
json.dump(res, open(os.path.join(OUT, "zz_listen.json"), "w"), indent=1, default=float)
npass = sum(v for _, v in verdict)
print(f"=== {npass} PASS / {len(verdict) - npass} FAIL  ({time.time() - t_start:.0f} s)", flush=True)
