"""Convergence study of the number of wave balls (the rule's only resolution parameters: η and ε).

  python lab/base_conv.py wave       field of the balls vs the exact RS integral, η = 0.5, 1, 2, 4 (2-D and 3-D)
  python lab/base_conv.py photons    double slit with photons, η = 0.5, 1, 2, 4 (ε = 1e-3)
  python lab/base_conv.py eps        double slit with photons, ε = 1e-2, 1e-3, 1e-4 (η = 2)

Nothing here changes the frozen rule: each case is the rule with one number changed, to show that above the
Nyquist count the results stop changing and below it they degrade.
"""
import sys
import time
import numpy as np
from base_rule import *
import optics as OP
from base_exp2d import ds_scene, run_double
from base_exp3d import rs3_points, polar_grid


def rule_with(**kw):
    r = dict(RULE)
    r.update(kw)
    return r


def wave():
    rows = []
    s0, w0 = OP.sources(2000)
    zs = (0.02, 0.1, 0.5, 2.0, 10.0, 160.0)
    refs = {}
    for z in zs:
        x = np.linspace(-6, 6, 241) if z < 20 else np.linspace(-60, 60, 481)
        U, Ux, Uz = OP.rs_field(x, z, s0, w0)
        refs[z] = (x, U[:, 0], np.imag(np.conj(U[:, 0]) * Uz[:, 0]))
    log("2-D double slit: max|U_balls − U_RS| / max|U_RS|  (and rms j_z error / max j_z) at distance z behind the slits")
    log("   eta   origins(l=0)  " + "  ".join(f"z={z:<6}" for z in zs))
    for eta in (0.5, 1.0, 2.0, 4.0):
        sc, _ = ds_scene(rule=rule_with(eta=eta))
        n0 = len(sc.origins(1)[0])
        errs = []
        for z in zs:
            x, Ur, jr = refs[z]
            U, Ux, Uz = sc.field(np.column_stack([x, np.full(x.size, z)]))
            jz = np.imag(np.conj(U[:, 0]) * Uz[:, 0])
            errs.append((np.abs(U[:, 0] - Ur).max() / np.abs(Ur).max(), np.sqrt(np.mean((jz - jr) ** 2)) / jr.max()))
        log(f"   {eta:<5} {n0:6d}        " + "  ".join(f"{a:.1e}/{b:.1e}" for a, b in errs))
        rows.append(dict(case=f"2-D eta={eta}", origins_slits=n0, err=errs))
    # 3-D: circular aperture D=6, plane wave: field on the screen (z=50) and near (z=0.05)
    log("3-D circular aperture D=6: relative error of |E_x|² and of j_z vs exact vector RS (fine polar quadrature)")
    xs, ys, wq = polar_grid(lambda r: np.ones_like(r), 0, 3.0, 600, 512)
    for z, rr in ((0.05, np.linspace(0, 3.5, 36)), (2.0, np.linspace(0, 4, 41)), (50.0, np.linspace(0, 25, 51))):
        P = np.column_stack([rr, np.zeros_like(rr), np.full(rr.size, z)])
        jr, Er = rs3_points(xs, ys, wq, P)
        line = []
        for eta in (0.5, 1.0, 2.0, 4.0):
            sc = Scene3D(rule_with(eta=eta))
            sc.add(Obj3("ap", "source", 0.0, mask=("circ_in", 0, 0, 3.0), src=dict(w0=np.inf)))
            j, rho, E = sc.current(P, 0)
            line.append((eta, len(sc.origins(0)[0]), np.abs(E[:, 0] - Er).max() / np.abs(Er).max(),
                         np.sqrt(np.mean((j[:, 2] - jr) ** 2)) / np.abs(jr).max()))
        log(f"   z={z}: " + "   ".join(f"η={e}: {n} origins, |E| err {a:.1e}, j_z err {b:.1e}" for e, n, a, b in line))
        rows.append(dict(case=f"3-D z={z}", data=line))
    save_results("conv_wave", rows, {})


def photons(kind="eta", vals=None):
    rows = []
    ref = OP.screen_reference("flux")
    if vals is None:
        vals = (0.5, 1.0, 2.0, 4.0) if kind == "eta" else (1e-2, 1e-3, 1e-4)
    for v in vals:
        r = rule_with(eta=v) if kind == "eta" else rule_with(eps=v)
        t = time.time()
        sc, res, xl = run_double(12000, seed=7, rule=r, verbose=False)
        c = OP.coarse_corr(xl, ref)
        e, rms, mx = OP.fringe_errors(xl, ref)
        chi = OP.chi2(xl, ref)
        ab = 1 - len(xl) / res["n"]
        bs = balls_summary(sc)
        log(f"  {kind}={v}: lmax {sc.d['lmax']}, balls/photon (level 0) {bs['balls_level0']}, all levels {bs['balls_all_levels']:.3g}; "
            f"{len(xl)} clicks: corr {c:.4f}  chi2/dof {chi:.2f}  fringe rms {100 * rms:.1f}% max {100 * mx:.1f}%  "
            f"per fringe {np.round(100 * e).astype(int).tolist()}  absorbed {ab:.4f}  [{time.time() - t:.0f}s]")
        rows.append(dict(**{kind: v}, clicks=len(xl), corr=c, chi2=chi, rms=rms, max=mx, absorbed=ab, balls=bs,
                         runtime=time.time() - t))
    b = OP.ideal_baseline(ref, int(np.mean([q["clicks"] for q in rows])), 200)
    log(f"  ideal sampler with the same number of clicks: corr {b[0]:.4f}±{b[1]:.4f}, fringe rms {100 * b[2]:.1f}%; exact absorbed 0.7806")
    save_results(f"conv_{kind}" + ("" if len(vals) in (3, 4) else "_" + "_".join(str(v) for v in vals)), rows, dict(ideal=b))


if __name__ == "__main__":
    w = sys.argv[1] if len(sys.argv) > 1 else "wave"
    t = time.time()
    vals = tuple(float(a) for a in sys.argv[2:]) or None
    {"wave": wave, "photons": lambda: photons("eta", vals), "eps": lambda: photons("eps", vals)}[w]()
    log(f"== conv {w} [{time.time() - t:.0f}s]")
