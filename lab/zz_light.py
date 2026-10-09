"""zz_light.py - T1: LIGHT in the zig-zag rule (spin-1 walkers, theta_m = 0, no charge).

The walker carries the Riemann-Silberstein vector F = E + i c B (Bialynicki-Birula 1994/1996).
One tick = translate along x,y,z (order alternated), the S_j = +1 part one cell forward, the -1 part one
cell back, the 0 part stays. No other ingredient. Units: lattice spacing a = 1, tick = a/c = 1.

Tests (tolerances in TOL, fixed before the first run):
  L1 dispersion of the lattice vs omega = c|k|: relative error vs lambda/a along axis / face / body
     diagonal (Bloch analysis of the two-tick walk); PASS if <= 2 % at lambda = 8a and the error falls
     like (a/lambda)^2 (ratio 8a -> 16a in [3, 5]).
  L2 propagation at c: centroid of the energy density of a 3-D packet (lambda = 8a and 16a) moving along
     the three directions; PASS if |v/c - v_group(Bloch)/c| < 1e-2 and |v/c - 1| < 2 % at 8a, 0.5 % at 16a.
  L3 transversality: a packet with div F = 0 (physical light). Energy fraction in the walk's
     non-moving (frequency 0) mode, and |k.F| / (|k||F|) after the run. PASS if static fraction < 1 %
     at lambda = 8a and falls by >= 8x at 16a.
  L4 spin 1: Berry phase of the moving (helicity +1) mode carried around loops of directions
     (great circle -> 2 pi = +1; cone of half angle 60 deg -> 2 pi (1 - cos 60) = pi). Tolerance 1e-2 rad.
     (Measured for real light in a coiled fibre: Tomita & Chiao, PRL 57, 937 (1986).)
  L5 (moved to zz_slit.py) Young double slit (the same geometry as docs/REGOLA_BASE.md: lambda = 1, slits 2 wide 8 apart,
     Gaussian beam waist 7 on the slits, screen at L = 160): time-averaged Poynting flux of the walk on
     the screen vs the exact 2-D Rayleigh-Sommerfeld flux (lab/optics.py). PASS if coarse correlation
     >= 0.98 and bright-fringe rms error <= 8 % (same thresholds as the photon rule) at lambda = 8a.
Run: timeout 590 python lab/zz_light.py
"""
import json, os, sys, time
import numpy as np
from zz_core import Walk, bloch_two_ticks, bloch_batch, quasi_energies, BASIS_ONE, SPIN1, SEED

sys.path.insert(0, os.path.dirname(__file__))
TOL = dict(disp8=0.02, disp_ratio=(3.0, 5.0), vg_match=1e-2, v8=0.02, v16=0.005,
           static8=0.01, static_drop=8.0, berry=1e-2, corr=0.98, fringe_rms=0.08)
OUT = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUT, exist_ok=True)
res, verdict = {}, []
t_start = time.time()
DIRS = {"axis (1,0,0)": np.array([1.0, 0, 0]), "face (1,1,0)": np.array([1.0, 1, 0]) / np.sqrt(2),
        "body (1,1,1)": np.array([1.0, 1, 1]) / np.sqrt(3)}


def check(name, ok, msg):
    verdict.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def light_modes(ka):
    """frequencies (per tick) and eigenvectors of the two-tick light walk; index of the static mode"""
    om, v = quasi_energies(ka, 0.0, 1)
    i0 = int(np.argmin(np.abs(om)))
    return om, v, i0


# ---------------- L1 dispersion ----------------
print("=== L1 dispersion of the light walk (Bloch analysis) ===", flush=True)
res["L1"] = {}
lams = [4, 6, 8, 12, 16, 32]
errs = {}
for dn, d in DIRS.items():
    row = []
    for lam in lams:
        k = 2 * np.pi / lam
        om, v, i0 = light_modes(k * d)
        wmax = np.max(np.abs(om))
        row.append(wmax / k - 1)
    errs[dn] = row
    print(f"  {dn:14s} " + "  ".join(f"l/a={l:2d}: {100 * e:+.4f}%" for l, e in zip(lams, row)), flush=True)
res["L1"]["lambda_over_a"] = lams
res["L1"]["rel_err"] = errs
e8 = max(abs(errs[d][lams.index(8)]) for d in DIRS)
e16 = max(abs(errs[d][lams.index(16)]) for d in DIRS)
ratio = e8 / e16 if e16 > 0 else np.inf
check("L1 dispersion error at lambda = 8a", e8 <= TOL["disp8"], f"max |omega/ck - 1| = {100 * e8:.3f} %")
check("L1 error ~ (a/lambda)^2", TOL["disp_ratio"][0] <= ratio <= TOL["disp_ratio"][1], f"ratio 8a/16a = {ratio:.2f} (2nd order: 4)")


# ---------------- L2 + L3 packets in 3-D ----------------
def transverse_packet(W, sigma, k0vec):
    """helicity +1 (physical, div F = 0) packet built in k-space: F(k) = A(k) e_+(k_hat)"""
    n = W.n
    x = W.x1
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    c0 = -15.0 * k0vec / np.linalg.norm(k0vec)
    env = np.exp(-((X - c0[0])**2 + (Y - c0[1])**2 + (Z - c0[2])**2) / (4 * sigma**2)
                 + 1j * (k0vec[0] * X + k0vec[1] * Y + k0vec[2] * Z))
    Ek = np.fft.fftn(env)
    kk = 2 * np.pi * np.fft.fftfreq(n)
    KX, KY, KZ = np.meshgrid(kk, kk, kk, indexing="ij")
    K = np.stack([KX, KY, KZ], -1)
    kn = np.linalg.norm(K, axis=-1)
    kn[kn == 0] = 1e-30
    Sk = np.einsum("...i,ijk->...jk", K / kn[..., None], np.array(SPIN1))
    w, v = np.linalg.eigh(Sk.reshape(-1, 3, 3))
    e_plus = v[:, :, 2]                               # eigenvalue +1 of S.k_hat: moves along +k_hat
    # gauge-free: project ONE fixed polarization (helicity + at the central k0) with P+(k) = e+ e+^dag
    w0, v0 = np.linalg.eigh(np.einsum("i,ijk->jk", k0vec / np.linalg.norm(k0vec), np.array(SPIN1)))
    eps0 = v0[:, 2]
    proj = e_plus * (e_plus.conj() @ eps0)[:, None]
    Fk = Ek[..., None] * proj.reshape(n, n, n, 3)
    F = np.moveaxis(np.fft.ifftn(Fk, axes=(0, 1, 2)), -1, 0)
    return F / np.sqrt(np.sum(np.abs(F)**2)), (X, Y, Z)


def static_fraction(F):
    """energy fraction in the walk's frequency-0 (non-moving) Bloch mode, and |k.F|/(|k||F|)"""
    n = F.shape[1]
    Fk = np.fft.fftn(F, axes=(1, 2, 3)).reshape(3, -1).T
    kk = 2 * np.pi * np.fft.fftfreq(n)
    KX, KY, KZ = np.meshgrid(kk, kk, kk, indexing="ij")
    K = np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], -1)
    wgt = np.sum(np.abs(Fk)**2, axis=1)
    sel = wgt > 1e-12 * wgt.max()
    tot = wgt[sel].sum()
    Wb = bloch_batch(K[sel], 0.0, 1)
    w, v = np.linalg.eig(Wb)
    i0 = np.argmin(np.abs(np.angle(w)), axis=1)                    # the frequency-0 mode
    v0 = v[np.arange(len(i0)), :, i0]
    v0 /= np.linalg.norm(v0, axis=1)[:, None]
    stat = np.sum(np.abs(np.sum(v0.conj() * Fk[sel], axis=1))**2)
    kn = np.linalg.norm(K[sel], axis=1)
    kn[kn == 0] = 1
    div = np.sqrt(np.sum(np.abs(np.sum(K[sel] * Fk[sel], axis=1))**2 / kn**2) / tot)
    return stat / tot, div


print("=== L2/L3 3-D packets (n = 64, sigma = 5a, 30 ticks) ===", flush=True)
res["L2"], res["L3"] = {}, {}
n, sigma, T = 64, 5.0, 30
for lam in (8, 16):
    k = 2 * np.pi / lam
    for dn, d in DIRS.items():
        W = Walk("light", n, 1.0)
        F, (X, Y, Z) = transverse_packet(W, sigma, k * d)
        F0 = F.copy()
        f0, div0 = static_fraction(F) if dn == "body (1,1,1)" else (None, None)

        def centroid(F):
            rho = np.sum(np.abs(F)**2, axis=0)
            return np.array([np.sum(rho * X), np.sum(rho * Y), np.sum(rho * Z)]) / rho.sum()
        r0 = centroid(F)
        for _ in range(T):
            F = W.tick(F)
        r1 = centroid(F)
        v = np.dot(r1 - r0, d) / T
        vperp = np.linalg.norm((r1 - r0) - np.dot(r1 - r0, d) * d) / T
        # group velocity of the moving mode at k0 from the Bloch dispersion (central difference)
        h = 1e-4
        g = []
        for q in range(3):
            e = np.zeros(3); e[q] = h
            g.append((np.max(light_modes(k * d + e)[0]) - np.max(light_modes(k * d - e)[0])) / (2 * h))
        vg = np.dot(np.array(g), d)
        # added after the second run (information): what this PACKET should do. A 5a-wide packet contains
        # many directions, so even in the continuum its centroid moves slower than c along d.
        Fk0 = np.fft.fftn(F0, axes=(1, 2, 3)).reshape(3, -1)
        wk = np.sum(np.abs(Fk0)**2, axis=0)
        sk = wk > 1e-10 * wk.max()
        kk1 = 2 * np.pi * np.fft.fftfreq(n)
        KK = np.stack([q.ravel() for q in np.meshgrid(kk1, kk1, kk1, indexing="ij")], -1)[sk]
        kn_ = np.linalg.norm(KK, axis=1); kn_[kn_ == 0] = 1
        v_cont = float(np.sum(wk[sk] * (KK @ d) / kn_) / wk[sk].sum())
        def wmax(Kq):
            return np.max(-np.angle(np.linalg.eigvals(bloch_batch(Kq, 0.0, 1))) / 2, axis=1)
        gw = np.zeros(len(KK))
        for q in range(3):
            e = np.zeros(3); e[q] = 1e-4
            gw += d[q] * (wmax(KK + e) - wmax(KK - e)) / 2e-4
        v_walkpred = float(np.sum(wk[sk] * gw) / wk[sk].sum())
        print(f"    (added) packet prediction: continuum Maxwell {v_cont:.5f} c, walk dispersion {v_walkpred:.5f} c, measured {v:.5f} c", flush=True)
        res["L2"][f"{lam}a {dn}"] = dict(v=v, v_perp=vperp, v_group_bloch=vg, v_packet_continuum=v_cont, v_packet_walk=v_walkpred)
        tol = TOL["v8"] if lam == 8 else TOL["v16"]
        check(f"L2 speed lambda={lam}a {dn}", abs(v - vg) < TOL["vg_match"] and abs(v - 1) < tol,
              f"centroid speed {v:.5f} c (Bloch group velocity {vg:.5f} c, sideways drift {vperp:.1e} c)")
        if f0 is not None:
            f1, div1 = static_fraction(F)
            res["L3"][f"{lam}a"] = dict(static_fraction=f0, div_initial=div0, div_after=div1, static_after=f1)
            print(f"  transversality lambda={lam}a (body diagonal): static-mode energy {f0:.2e} (after {T} ticks {f1:.2e}),"
                  f" |k.F|/|k||F| {div0:.1e} -> {div1:.2e}", flush=True)
s8, s16 = res["L3"]["8a"]["static_fraction"], res["L3"]["16a"]["static_fraction"]
check("L3 transversality at 8a", s8 < TOL["static8"], f"static (longitudinal-like) energy fraction {s8:.2e}")
check("L3 transversality improves with lambda/a", s8 / max(s16, 1e-300) >= TOL["static_drop"], f"8a/16a ratio {s8 / max(s16, 1e-300):.1f}")

# ---------------- L4 Berry phase (spin 1) ----------------
print("=== L4 Berry phase of the moving light mode around loops of directions ===", flush=True)


def berry(kmag, cone_deg, spin, theta=0.0, npts=400, branch="plus"):
    """discrete Berry phase of the upper (helicity +) moving band around a cone about z"""
    vs = []
    th = np.radians(cone_deg)
    for p in np.linspace(0, 2 * np.pi, npts, endpoint=False):
        kv = kmag * np.array([np.sin(th) * np.cos(p), np.sin(th) * np.sin(p), np.cos(th)])
        om, v = quasi_energies(kv, theta, spin)
        vs.append(v[:, -1] / np.linalg.norm(v[:, -1]))   # highest frequency = moving along +k (helicity +)
    prod = 1.0 + 0j
    for i in range(npts):
        prod *= np.vdot(vs[i], vs[(i + 1) % npts])
    return -np.angle(prod)


res["L4"] = {}
for cone, expect in ((90, 2 * np.pi), (60, np.pi)):
    g = berry(2 * np.pi / 16, cone, 1)
    dev = np.angle(np.exp(1j * (g - expect)))
    res["L4"][f"cone {cone}"] = dict(berry=g, expected=expect)
    check(f"L4 light Berry phase, cone {cone} deg", abs(dev) < TOL["berry"],
          f"gamma = {g:+.5f} rad (mod 2pi), expected spin 1 x solid angle = {expect:.5f} -> exp(i gamma) = {np.exp(1j * g).real:+.4f}")
# added after the first run (information, no test): lattice convergence of the 60-degree loop
res["L4"]["cone 60 convergence"] = {}
for lam in (8, 16, 32, 64):
    g = berry(2 * np.pi / lam, 60, 1)
    dev = float(np.angle(np.exp(1j * (g - np.pi))))
    res["L4"]["cone 60 convergence"][lam] = dev
    print(f"  cone 60 deg, lambda = {lam}a: gamma - pi = {dev:+.5f} rad", flush=True)


# L5 (double slit) moved to lab/zz_slit.py after the first run (numpy too slow here: 753 s estimated).

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from optics import XS
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    for dn in DIRS:
        ax[0].loglog(lams, np.abs(errs[dn]) + 1e-16, "o-", label=dn)
    ax[0].loglog(lams, 2.0 * (np.array(lams, float) / 8)**-2 * e8, "k:", label="~ (a/lambda)^2")
    ax[0].set_xlabel("lambda / a"); ax[0].set_ylabel("|omega/(ck) - 1|"); ax[0].legend(fontsize=7)
    ax[0].set_title("light walk: lattice dispersion")
    ax[1].axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "zz_light.png"), dpi=90)
except Exception as ex:
    print("figure skipped:", ex)

res["verdict"] = verdict
res["seconds"] = time.time() - t_start
with open(os.path.join(OUT, "zz_light.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
npass = sum(v for _, v in verdict)
print(f"=== {npass} PASS / {len(verdict) - npass} FAIL  ({time.time() - t_start:.0f} s)", flush=True)
