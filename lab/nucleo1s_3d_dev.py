"""nucleo1s_3d_dev.py - DEVELOPMENT runs of phase 2 (declared as such in docs/NUCLEO_1S.md, before the 3-D pre-registration):
the Laboratorio start state (radial 1s with the same a, on the cubic lattice, N = 48, a = 0.5, Z = 60) in each colour,
listened and harmonic-inverted.  Usage: python nucleo1s_3d_dev.py lab [nt]"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_3d import Colour3D, radial_state, spinor_field, fdm, fdm_stable, dirac_omega, pencil  # noqa
from zz_core import ALPHA

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
Za = 60 * ALPHA


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


T0 = time.time()
if sys.argv[1] == "lab":
    nt = int(sys.argv[2]) if len(sys.argv) > 2 else 2400
    n, a = 48, 0.5
    om_r, r, g, f = radial_state(-1, 1, a, Za)
    log("radial 1s (same a):", om_r)
    res = {}
    for q in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)):
        col = Colour3D(n, a, Za, q=q)
        E, O = col.field(spinor_field(-1, r, g, f, a))
        C, _, _ = col.listen(E, O, nt, log=log)
        roots = fdm_stable(C, a, 0.75, 1.02, dmin=1e-3)
        log("colour", q, [(round(x["omega"], 5), round(x["weight"], 4), x["ok"]) for x in roots])
        res[str(q)] = dict(C_re=C.real.tolist(), C_im=C.imag.tolist(), roots=roots)
    res["dirac_1s"] = dirac_omega(1, -1, Za)
    json.dump(res, open(os.path.join(OUT, "nucleo1s_3d_dev_lab.json"), "w"))

if sys.argv[1] == "states":
    # filtered states of the four lines (colour 0) and what they look like
    nt = 4000
    n, a = 48, 0.5
    om_r, r, g, f = radial_state(-1, 1, a, Za)
    col = Colour3D(n, a, Za)
    E, O = col.field(spinor_field(-1, r, g, f, a))
    oms = [0.84037, 0.88627, 0.8993, 0.8996]
    C, _, acc = col.listen(E, O, nt, filt=oms, log=log)
    info = {}
    for om, (ae, ao) in zip(oms, acc):
        nr = np.sqrt(np.sum(np.abs(ae)**2) + np.sum(np.abs(ao)**2))
        ae, ao = ae / nr, ao / nr
        rho_e = np.sum(np.abs(ae)**2, -1); rho_o = np.sum(np.abs(ao)**2, -1)
        rr = np.concatenate([col.rE.ravel(), col.rO.ravel()]); rho = np.concatenate([rho_e.ravel(), rho_o.ravel()])
        # nearest two sites
        k2 = np.argsort(rr)[:2]
        # radial shells
        sh = {}
        for R in (0.5, 1.0, 2.0, 4.0):
            sh[R] = float(rho[rr < R].sum())
        mean_r = float((rho * rr).sum())
        # angular: weight along the colour axis vs perpendicular, for r in 0.5..3
        Xe, Ye, Ze = col.coords(0); Xo, Yo, Zo = col.coords(1)
        X = np.concatenate([Xe.ravel(), Xo.ravel()]); Y = np.concatenate([Ye.ravel(), Yo.ravel()]); Z = np.concatenate([Ze.ravel(), Zo.ravel()])
        near = np.argsort(rr)[:2]
        ax = np.array([X[near[0]], Y[near[0]], Z[near[0]]]); ax /= np.linalg.norm(ax)
        cosax = (X * ax[0] + Y * ax[1] + Z * ax[2]) / np.maximum(rr, 1e-9)
        m = (rr > 0.5) & (rr < 3)
        P2 = float(np.sum(rho[m] * 0.5 * (3 * cosax[m]**2 - 1)) / np.sum(rho[m]))
        # Rayleigh frequency of the filtered state
        E2, O2 = np.empty_like(ae), np.empty_like(ao)
        col.ticks = 0
        e1, o1 = col.tick(ae.copy(), ao.copy(), E2, O2)
        E3, O3 = np.empty_like(ae), np.empty_like(ao)
        e2, o2 = col.tick(e1.copy(), o1.copy(), E3, O3)
        lam = (np.vdot(ae, e2) + np.vdot(ao, o2))
        info[om] = dict(rayleigh=float(-np.angle(lam) / (2 * a)), absl=float(abs(lam)), two_nearest=float(rho[k2].sum()),
                        shells=sh, mean_r=mean_r, P2_axis=P2, axis=ax.tolist())
        log(om, info[om])
    json.dump(info, open(os.path.join(OUT, "nucleo1s_3d_dev_states.json"), "w"), indent=1)

if sys.argv[1] == "scan":
    # Lab-like start (radial level of the same a) on lattices a = 0.4, 0.3, 0.24 (box 24), colour 0
    a = float(sys.argv[2]); nt = int(sys.argv[3])
    n = int(round(24 / a / 2)) * 2
    om_r, r, g, f = radial_state(-1, 1, a, Za)
    col = Colour3D(n, a, Za)
    E, O = col.field(spinor_field(-1, r, g, f, a))
    C, _, _ = col.listen(E, O, nt, log=log)
    roots = fdm_stable(C, a, 0.70, 1.02, dmin=1e-3)
    log("a", a, "n", n, "radial", om_r, [(round(x["omega"], 5), round(x["weight"], 4), x["ok"]) for x in roots])
    json.dump(dict(a=a, n=n, radial=om_r, roots=roots, C_re=C.real.tolist(), C_im=C.imag.tolist()),
              open(os.path.join(OUT, f"nucleo1s_3d_dev_scan_{a}.json"), "w"))

if sys.argv[1] == "taste":
    # the four lines: what are they? long filter at a = 0.5 (Lab lattice), colour 0
    nt = 8000
    n, a = 48, 0.5
    om_r, r, g, f = radial_state(-1, 1, a, Za)
    col = Colour3D(n, a, Za)
    E, O = col.field(spinor_field(-1, r, g, f, a))
    oms = [0.84037, 0.88627, 0.8995]
    C, _, acc = col.listen(E, O, nt, filt=oms, log=log)
    m = n // 2
    kk = np.fft.fftfreq(m) * m          # integer wave numbers
    K = np.abs(kk) > m / 4             # "high" half of the compact zone along an axis
    info = {}
    for om, (ae, ao) in zip(oms, acc):
        nr = np.sqrt(np.sum(np.abs(ae)**2) + np.sum(np.abs(ao)**2)); ae, ao = ae / nr, ao / nr
        d = {}
        d["E_weight"] = float(np.sum(np.abs(ae)**2))
        d["L_weight"] = float(np.sum(np.abs(ae[..., :2])**2) + np.sum(np.abs(ao[..., :2])**2))
        for nm, arr in (("E", ae), ("O", ao)):
            F = np.sum(np.abs(np.fft.fftn(arr, axes=(0, 1, 2)))**2, -1)
            F /= F.sum()
            pat = {}
            for hx in (0, 1):
                for hy in (0, 1):
                    for hz in (0, 1):
                        mk = (K[:, None, None] == hx) & (K[None, :, None] == hy) & (K[None, None, :] == hz)
                        pat[f"{hx}{hy}{hz}"] = round(float(F[mk].sum()), 4)
            d["kpattern_" + nm] = pat
        # overlap with start, Rayleigh with two ticks
        col.ticks = 0
        e1, o1 = col.tick(ae.copy(), ao.copy(), np.empty_like(ae), np.empty_like(ao))
        e2, o2 = col.tick(e1.copy(), o1.copy(), np.empty_like(ae), np.empty_like(ao))
        lam = np.vdot(ae, e2) + np.vdot(ao, o2)
        d["rayleigh"] = float(-np.angle(lam) / (2 * a)); d["abs_lam"] = float(abs(lam))
        rho = np.concatenate([np.sum(np.abs(ae)**2, -1).ravel(), np.sum(np.abs(ao)**2, -1).ravel()])
        rr = np.concatenate([col.rE.ravel(), col.rO.ravel()])
        d["two_nearest"] = float(rho[np.argsort(rr)[:2]].sum()); d["mean_r"] = float((rho * rr).sum())
        log(om, d)
        info[om] = d
    json.dump(info, open(os.path.join(OUT, "nucleo1s_3d_dev_taste.json"), "w"), indent=1)

if sys.argv[1] == "nuc":
    # nucleus agglomerate exploration: the two nearest sites of the colour get V0 = s0 * (-Za/r0), next shell s1
    a = float(sys.argv[2]); nt = int(sys.argv[3]); s0 = float(sys.argv[4]); s1 = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    n = int(round(24 / a / 2)) * 2
    r0 = np.sqrt(3) / 2 * a
    r1 = np.sqrt(11) / 2 * a

    def Vf(rr):
        V = -Za / rr
        V = np.where(np.abs(rr - r0) < 1e-6 * a, s0 * V, V)
        V = np.where(np.abs(rr - r1) < 1e-6 * a, s1 * V, V)
        return V
    om_r, r, g, f = radial_state(-1, 1, a, Za)
    col = Colour3D(n, a, Za, Vfun=Vf)
    E, O = col.field(spinor_field(-1, r, g, f, a))
    C, _, _ = col.listen(E, O, nt)
    roots = [x for x in pencil(C, a, 0.70, 0.98) if x[1] > 3e-3]
    log("a", a, "s0", s0, "s1", s1, [(round(x[0], 5), round(x[1], 3)) for x in roots])
