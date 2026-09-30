"""TASK 1 — why did guidance by the EXACT field give only 0.95?

We guide the photon ball with the exact field and switch on, one at a time, the
ingredients that the old diagnosis (exact_guide.py) got wrong:

  law     'bohm-paraxial' : dx/dz = Im(U* dxU) / (k |U|^2)          (old)
          'flux'          : dx/dz = j_x / j_z  (photon moves at c along the current)
  field   'legacy'        : e^{ikr}/sqrt(r)*(1+z/r)/2, dxU from phase only (old)
          'rs'            : exact 2-D Rayleigh-Sommerfeld (Hankel), exact gradient
  start   'inc'           : |U_inc|^2 inside the slits at z0 (old)
          'jz'            : flux j_z(x, z0) of the field used for guidance (equivariant)
  ref     'paraxial'      : paraxial Fresnel (old reference)
          'flux'          : exact j_z on the screen (what a flat detector counts)

Run:  python lab/exact_guide2.py          (steps 5 and 7 = the two good ones, ~3.5 min)
      python lab/exact_guide2.py 0 3     (any steps; one or two per command keeps it < 5 min)
"""
import numpy as np, time, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
from multiprocessing import Pool
from optics import *

NPH = 6000


def legacy_field(x, z, s, w):
    r = np.sqrt((x[:, None] - s[None, :]) ** 2 + z ** 2)
    e = np.exp(1j * K * r) / np.sqrt(r) * (0.5 * (1 + z / r)) * w[:, 0][None, :]
    U = e.sum(1)
    Ux = (e * (1j * K * (x[:, None] - s[None, :]) / r)).sum(1)
    Uz = (e * (1j * K * z / r)).sum(1)
    return U[:, None], Ux[:, None], Uz[:, None]


def slope(x, z, cfg, s, w):
    f = legacy_field if cfg["field"] == "legacy" else rs_field
    U, Ux, Uz = f(x, z, s, w)
    jx, jz = flux(U, Ux, Uz)
    if cfg["law"] == "flux":
        return jx / jz
    return jx / (K * (np.abs(U) ** 2).sum(-1))


def integrate(args):
    """RK2/RK4 in z. Step dz(z) = min(dz, max(dzmin, grow*z)): small near the slits
    (near-field structure on the scale of z), up to dz in the far field."""
    x, cfg = args
    near = sources(cfg["m"]); far = sources(min(cfg["m"], 100))   # ds = 0.02 is plenty once z > 1
    z = cfg["z0"]
    while z < L - 1e-9:
        s, w = near if z < 1.0 else far
        dz = min(cfg["dz"], max(cfg.get("dzmin", cfg["dz"]), cfg.get("grow", 1e9) * z), L - z)
        if cfg["rk"] == 2:
            k1 = slope(x, z, cfg, s, w)
            x = x + dz * slope(x + 0.5 * dz * k1, z + 0.5 * dz, cfg, s, w)
        else:
            k1 = slope(x, z, cfg, s, w)
            k2 = slope(x + 0.5 * dz * k1, z + 0.5 * dz, cfg, s, w)
            k3 = slope(x + 0.5 * dz * k2, z + 0.5 * dz, cfg, s, w)
            k4 = slope(x + dz * k3, z + dz, cfg, s, w)
            x = x + dz / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        z += dz
    return x


def start(cfg, rng):
    if cfg["start"] == "inc":
        x = rng.normal(0, W0 / 2, 60 * NPH)
        return x[in_slit(x)][:NPH]
    s, w = sources(cfg["m"])
    if cfg["start"] == "thin":
        # photon enters with |U_inc|^2; it survives the edge with prob min(1, j_z/(k|U_inc|^2)) at z0
        x = rng.normal(0, W0 / 2, 80 * NPH)
        x = x[in_slit(x)][:3 * NPH]
        U, Ux, Uz = rs_field(x, cfg["z0"], s, w)
        r = flux(U, Ux, Uz)[1] / (K * inc_amp(x) ** 2)
        keep = rng.random(x.size) < np.minimum(1, r)
        print(f"    edge rule: survival {keep.mean():.4f} (mean r {r.mean():.4f}, r>1 for {100*(r>1).mean():.0f}% of entries)", flush=True)
        return x[keep][:NPH]
    xg = np.arange(-7, 7, 0.002)
    f = legacy_field if cfg["field"] == "legacy" else rs_field
    U, Ux, Uz = f(xg, cfg["z0"], s, w)
    p = np.clip(flux(U, Ux, Uz)[1], 0, None)
    c = np.cumsum(p); c /= c[-1]
    return np.interp(rng.random(NPH), c, xg)


def run(cfg, pool, refs):
    rng = np.random.default_rng(cfg.get("seed", 1))
    x0 = start(cfg, rng)
    t = time.time()
    parts = pool.map(integrate, [(c, cfg) for c in np.array_split(x0, 32)])
    xl = np.concatenate(parts)
    tag = " ".join(f"{k}={v}" for k, v in cfg.items() if k != "seed")
    out = {}
    for rn in cfg["refs"]:
        out[rn] = report(f"  [{tag}] vs {rn:8s}", xl, refs[rn], f"[{time.time()-t:.0f}s]")
    return xl, out


if __name__ == "__main__":
    refs = {"paraxial": screen_reference("paraxial"), "flux": screen_reference("flux")}
    b = ideal_baseline(refs["flux"], NPH, 200)
    print(f"ideal sampler, {NPH} clicks: corr {b[0]:.4f}+-{b[1]:.4f}, fringe err rms {100*b[2]:.1f}% max {100*b[3]:.1f}%", flush=True)
    base = dict(law="bohm-paraxial", field="legacy", start="inc", z0=0.05, dz=0.5, rk=2, m=600,
                refs=("paraxial", "flux"))
    AD = dict(rk=4, dz=2.0, dzmin=0.01, grow=0.1)      # adaptive RK4
    steps = [
        ("0 old diagnosis (exact_guide.py)", {}),
        ("1 old field+start, flux law dx/dz = jx/jz", dict(law="flux", m=300)),
        ("2 exact RS field, flux law, old start |inc|^2 at z0=0.05", dict(law="flux", field="rs", m=300, **AD)),
        ("3 exact RS field, flux law, start j_z(z0=0.05)", dict(law="flux", field="rs", start="jz", m=300, **AD)),
        ("4 exact RS, flux, start j_z(z0=0.5)", dict(law="flux", field="rs", start="jz", z0=0.5, m=100, **AD)),
        ("5 exact RS, flux, start j_z(z0=2)", dict(law="flux", field="rs", start="jz", z0=2.0, m=100, **AD)),
        ("6 as 5, half steps (convergence)", dict(law="flux", field="rs", start="jz", z0=2.0, m=100,
                                                 rk=4, dz=1.0, dzmin=0.005, grow=0.05)),
        ("7 exact RS, flux, |inc|^2 start + local edge rule at z0=0.05", dict(law="flux", field="rs", start="thin", m=300, **AD)),
        ("8 as 5 but OLD law dx/dz = jx/(k|U|^2)", dict(field="rs", start="jz", z0=2.0, m=100, **AD)),
    ]
    only = [int(a) for a in sys.argv[1:]] or [5, 7]
    with Pool(8) as pool:
        for i, (name, over) in enumerate(steps):
            if only and i not in only:
                continue
            print(name, flush=True)
            cfg = dict(base); cfg.update(over)
            xl, _ = run(cfg, pool, refs)
            os.makedirs(os.path.join(HERE, "out"), exist_ok=True)
            np.save(os.path.join(HERE, "out", f"exact_guide2_step{i}.npy"), xl)
