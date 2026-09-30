"""TASK 3 — the SAME rule (photon_rule.py, estimator E2, lattice fan, hat window, flux-ratio rule at
objects) on the other experiments.

  free     : free Gaussian beam, width 2·std(x) of the photon balls vs w0·sqrt(1+(z/zR)^2)
  which    : double slit with a polarizer on each slit (parallel H,H / orthogonal H,V; light enters at 45°)
  malus    : polarizer at angle θ on a free beam
             M1 absorbing polarizer  -> flux-ratio rule at the object (local coin of the photon ball)
             M2 polarizing beam splitter + absorber on the rejected port -> deterministic Bohmian branching
  build    : single-photon build-up (clicks one by one: correlation after 50, 200, 1000, 6000 photons)

Run: python lab/validate.py free | which | malus1 | malus2 | build      (each < 4 min on 8 cores)
"""
import sys, os, time, numpy as np
from multiprocessing import Pool
from optics import *
from photon_rule import (Fan, lattice_fan, field_lattice, guide, double_slit, edge_rule, K)

HERE = os.path.dirname(os.path.abspath(__file__))
SMOOTH = dict(ds=0.1, delta=0.5, dth=0.01)          # sources without sharp edges


def flux_z(U, Uz):
    return np.imag(np.conj(U) * Uz).sum(-1)


# ─────────────────────────────── free Gaussian beam ─────────────────────────────────
W0F = 3.0
ZR = np.pi * W0F ** 2 / LAM


def source_fan(w0, jones=(1.0,), p=SMOOTH):
    J = np.asarray(jones, complex)
    return lattice_fan(-3.5 * w0, 3.5 * w0, p["ds"], lambda s: np.exp(-s ** 2 / w0 ** 2), dth=p["dth"],
                       jones=[J])


def _free(args):
    x0, seed = args
    fan = source_fan(W0F)
    rng = np.random.default_rng(seed)
    keep, r = edge_rule(x0, fan, SMOOTH["delta"], np.exp(-x0 ** 2 / W0F ** 2), "E2", rng)
    X = x0[keep]; z = SMOOTH["delta"]; out = []
    for zt in (ZR, 2 * ZR, 5 * ZR):
        X, _ = guide(X, fan, z, zt, "E2")
        z = zt; out.append(X.copy())
    return out, (~keep).sum()


def run_free(pool, n=4000):
    rng = np.random.default_rng(11)
    x0 = rng.normal(0, W0F / 2, n)                       # photon ball emitted with |U|^2
    res = pool.map(_free, [(c, 100 + i) for i, c in enumerate(np.array_split(x0, 32))])
    lost = sum(r[1] for r in res)
    print(f"free beam: w0={W0F}, zR={ZR:.2f}, {n} photons, lost at the source by the flux rule: {lost}", flush=True)
    ok = True
    # exact RS second moment for comparison (nonparaxial correction)
    s = np.linspace(-3.5 * W0F, 3.5 * W0F, 1401); w = (np.exp(-s ** 2 / W0F ** 2) * (s[1] - s[0]))[:, None].astype(complex)
    for i, (zt, name) in enumerate(((ZR, "zR"), (2 * ZR, "2zR"), (5 * ZR, "5zR"))):
        X = np.concatenate([r[0][i] for r in res])
        wm = 2 * X.std()
        wth = W0F * np.sqrt(1 + (zt / ZR) ** 2)
        xg = np.linspace(-6 * wth, 6 * wth, 3001)
        U, Ux, Uz = rs_field(xg, zt, s, w); jz = flux(U, Ux, Uz)[1]
        wrs = 2 * np.sqrt((jz * xg ** 2).sum() / jz.sum())
        err = wm / wth - 1
        ok &= abs(err) < 0.03
        print(f"  z={name:4s} ({zt:6.1f}): width measured {wm:6.3f}  Gaussian-beam theory {wth:6.3f}  "
              f"(exact RS {wrs:6.3f})  error {100*err:+.2f}%  (stat. ±{100/np.sqrt(2*X.size):.1f}%)  "
              f"{'PASS' if abs(err) < 0.03 else 'FAIL'}", flush=True)
    return ok


# ─────────────────────────────── which path ─────────────────────────────────────────
def run_which(pool):
    Hj, Vj = np.array([1, 0], complex), np.array([0, 1], complex)
    Din = np.array([1, 1], complex) / np.sqrt(2)
    ok = True
    for name, f in (("parallel H,H", (Hj, Hj)), ("orthogonal H,V", (Hj, Vj))):
        # polarizer f on each slit: the balls through slit i carry f_i <f_i|Din>
        jones = [fi * np.vdot(fi, Din) for fi in f]
        ref = screen_reference("flux", jones=jones)
        xl, ab, dt = double_slit(dict(jones=jones, nph=12000), pool)   # half is absorbed by the slit polarizers
        np.save(os.path.join(HERE, "out", f"which_{name.split()[0]}.npy"), xl)
        V = visibility(xl); Vth = visibility(ref, True)
        b = ideal_baseline(ref, xl.size, 100)
        print(f"  (ideal sampler with {xl.size} clicks: corr {b[0]:.4f}+-{b[1]:.4f}, fringe rms {100*b[2]:.1f}%)", flush=True)
        c = report(f"which-path {name}", xl, ref, f"absorbed {100*ab:.1f}% [{dt:.0f}s]")
        good = (V < 0.15) if "orth" in name else (V > 0.9)
        ok &= good
        print(f"  visibility measured {V:.3f}  theory {Vth:.3f}  -> {'PASS' if good else 'FAIL'}", flush=True)
    return ok


# ─────────────────────────────── Malus ──────────────────────────────────────────────
ZP, ZEND_M1, ZEND_M2, ALPHA = 5.0, 40.0, 200.0, 0.25
GUIDE_M = dict(dzmax=float(os.environ.get("MALUS_DZMAX", 5.0)), grow=float(os.environ.get("MALUS_GROW", 0.1)))


def refan(fan_in, zp, transform, half=12.0, p=SMOOTH):
    """Re-emission at an object on the plane zp: each ball that hits it becomes a new origin whose
    amplitude = local field estimate at the hit point (from its neighbours) x object transmission."""
    n = int(round(2 * half / p["ds"]))
    s = -half + (np.arange(n) + 0.5) * 2 * half / n
    U, _, _ = field_lattice(s, zp, fan_in, "E2")
    a = transform(s, U) * (2 * half / n)
    return Fan(s, a, zp, dth=p["dth"])


def polarizer_T(theta):
    ax = np.array([np.cos(theta), np.sin(theta)], complex)
    return lambda s, U: (U @ ax.conj())[:, None] * ax[None, :]


def pbs_T(theta, alpha=ALPHA):
    ax = np.array([np.cos(theta), np.sin(theta)], complex)
    bx = np.array([-np.sin(theta), np.cos(theta)], complex)
    g = 1 / np.sqrt(np.cos(alpha))                    # thin-prism transmission: keeps the z-flux
    def T(s, U):
        return g * ((U @ ax.conj())[:, None] * ax[None, :] * np.exp(1j * K * np.sin(alpha) * s)[:, None]
                    + (U @ bx.conj())[:, None] * bx[None, :] * np.exp(-1j * K * np.sin(alpha) * s)[:, None])
    return T


def _malus(args):
    x0, theta, kind, seed = args
    rng = np.random.default_rng(seed)
    H = (1.0, 0.0)
    fan1 = source_fan(W0F, H)
    keep, _ = edge_rule(x0, fan1, SMOOTH["delta"], np.exp(-x0 ** 2 / W0F ** 2), "E2", rng)
    X = x0[keep]; lost_src = (~keep).sum()
    X, u = guide(X, fan1, SMOOTH["delta"], ZP, "E2")
    # flux arriving at the object, read locally
    U, Ux, Uz = field_lattice(X, ZP, fan1, "E2"); jin = flux_z(U, Uz)
    fan2 = refan(fan1, ZP, polarizer_T(theta) if kind == "M1" else pbs_T(theta))
    # first step behind the object: move along the OUTGOING current (implicit Euler, 2 iterations);
    # flying straight here would shift a tilted branch by delta*alpha and bias the branching.
    from photon_rule import slope as _slope
    Xd = X + SMOOTH["delta"] * u
    for _ in range(2):
        Xd = X + SMOOTH["delta"] * _slope(Xd, ZP + SMOOTH["delta"], fan2, "E2", u)
    U2, _, Uz2 = field_lattice(Xd, ZP + SMOOTH["delta"], fan2, "E2"); jout = flux_z(U2, Uz2)
    r = jout / jin
    if kind == "M2":
        # lossless (unitary) element: |T psi|^2 >= |psi|^2 everywhere -> the rule never removes the photon.
        # (Evaluating j_out at a finite delta behind a TILTING element shifts the beam sideways by
        #  delta*alpha and makes r fluctuate +-15% for purely numerical reasons; see REGOLA_FOTONE.md.)
        survive = np.ones(X.size, bool)
    else:
        survive = rng.random(X.size) < np.minimum(1, r)
    X = Xd[survive]
    X, _ = guide(X, fan2, ZP + SMOOTH["delta"], ZEND_M1 if kind == "M1" else ZEND_M2, "E2", **GUIDE_M)
    ax = np.array([np.cos(theta), np.sin(theta)], complex)
    if kind == "M1":
        passed = X
    else:
        passed = X[X > 0]                               # rejected branch goes to x < 0 -> absorber
    # polarization of the local wave at the photon after passing
    Ue, _, _ = field_lattice(passed, ZEND_M1 if kind == "M1" else ZEND_M2, fan2, "E2")
    pol = (np.abs(Ue @ ax.conj()) ** 2 / np.maximum((np.abs(Ue) ** 2).sum(-1), 1e-300))
    return passed.size, x0.size, lost_src, (~survive).sum(), pol.sum(), r.min() if r.size else 1, r.max() if r.size else 1


def run_malus(pool, n=3000, kinds=("M1", "M2")):
    ok = True
    for kind in kinds:
        print(f"Malus {kind}: {'absorbing polarizer, flux-ratio rule' if kind=='M1' else 'PBS + absorber, deterministic branching'}", flush=True)
        for deg in [float(a) for a in os.environ.get("MALUS_ANGLES", "0,22.5,45,67.5,90").split(",")]:
            th = np.deg2rad(deg); t = time.time()
            rng = np.random.default_rng(int(deg * 10) + (0 if kind == "M1" else 7))
            x0 = rng.normal(0, W0F / 2, n)
            res = pool.map(_malus, [(c, th, kind, 1000 * i + int(deg)) for i, c in enumerate(np.array_split(x0, 16))])
            npass = sum(r[0] for r in res); ntot = sum(r[1] for r in res)
            lost = sum(r[2] for r in res)
            pol = sum(r[4] for r in res) / max(npass, 1)
            P = npass / ntot; Pth = np.cos(th) ** 2
            sig = np.sqrt(max(Pth * (1 - Pth), 1 / ntot) / ntot)
            good = abs(P - Pth) < max(3 * sig, 0.01)
            ok &= good
            print(f"  θ={deg:5.1f}°  P(pass) {P:.4f}  cos²θ {Pth:.4f}  (±{sig:.4f})  lost at source {lost}  "
                  f"r range [{min(r[5] for r in res):.3f},{max(r[6] for r in res):.3f}]  "
                  f"<|<a|ψ>|²> after = {pol:.4f}  {'PASS' if good else 'FAIL'}  [{time.time()-t:.0f}s]", flush=True)
    return ok


# ─────────────────────────────── build-up ───────────────────────────────────────────
def run_build(pool):
    ref = screen_reference("flux")
    xl, ab, dt = double_slit(dict(nph=6000, seed=21), pool)
    rng = np.random.default_rng(5); xl = rng.permutation(xl)      # photons arrive one at a time
    for n in (50, 200, 1000, 3000, len(xl)):
        base = ideal_baseline(ref, n, 100)
        print(f"  after {n:5d} photons: corr {coarse_corr(xl[:n], ref):.3f}   (ideal sampler {base[0]:.3f}±{base[1]:.3f})", flush=True)
    np.save(os.path.join(HERE, "out", "buildup_clicks.npy"), xl)
    return True


if __name__ == "__main__":
    os.makedirs(os.path.join(HERE, "out"), exist_ok=True)
    what = sys.argv[1:] or ["free", "which", "malus", "build"]
    with Pool(8) as pool:
        for w in what:
            t = time.time()
            ok = {"free": run_free, "which": run_which, "malus": run_malus, "build": run_build,
                  "malus1": lambda p: run_malus(p, kinds=("M1",)), "malus2": lambda p: run_malus(p, kinds=("M2",))}[w](pool)
            print(f"== {w}: {'PASS' if ok else 'FAIL'}  [{time.time()-t:.0f}s]", flush=True)
