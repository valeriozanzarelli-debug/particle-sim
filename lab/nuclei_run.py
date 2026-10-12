"""nuclei_run.py JOB [--wall S] - relax one configuration (N1 rational maps, N2 clusters) with L-BFGS (c4 static energy of
knot2_core, read-only), resumable in chunks; when converged (or at the iteration cap) writes lab/out/nuclei_JOB.json
with all the pre-registered measures and the final baryon/energy densities (npz) for the figure.
Jobs: rm1 rm2 rm3 rm4 (rational maps, L=48, h=0.2), rm4_L64 (box check), pair25 pair30 (B=2 from two hedgehogs),
tet4 sq4 (B=4 from four hedgehogs)."""
import math, json, os, sys, time, argparse, numpy as np, torch
import nuclei_core as nc, knot_core as kc, knot2_core as k2
torch.set_num_threads(2)
k2.set_scheme("c4")
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
H = 0.2; MAXIT = 600; TOL = 1e-9
PI2 = 12 * math.pi**2


def job_spec(job):
    if job.startswith("rm"):
        B = int(job[2]); L = 64 if job.endswith("L64") else 48
        return dict(kind="rm", B=B, L=L)
    if job.startswith("pair"):
        d = float(job[4:]) / 10
        return dict(kind="cluster", B=2, L=48, centers=[(-d / 2, 0, 0), (d / 2, 0, 0)], isos=[None, "z"], sep=d)
    if job == "tet4":
        a = 2.5 / (2 * math.sqrt(2))        # edge 2.5
        V = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
        return dict(kind="cluster", B=4, L=48, centers=[tuple(a * np.array(v)) for v in V], isos=[None, "x", "y", "z"], sep=2.5)
    if job == "sq4":
        s = 2.5 / 2
        return dict(kind="cluster", B=4, L=48, centers=[(s, s, 0), (-s, s, 0), (-s, -s, 0), (s, -s, 0)], isos=[None, "z", None, "z"], sep=2.5)
    raise SystemExit("unknown job " + job)


def initial(sp, log):
    L = sp["L"]; X, Y, Z = kc.grid((L, L, L), H)
    if sp["kind"] == "rm":
        b, I = nc.angular_integrals(sp["B"], 400)
        r, f, E = nc.rm_profile(sp["B"], I, n=3000)
        log(f"  rational map B={sp['B']}: angular B={b:.5f} I={I:.4f}  ansatz E/12pi^2B={E/(PI2*sp['B']):.5f}")
        return nc.rm_field(sp["B"], X, Y, Z, r, f), dict(I=I, E_ansatz_over=E / (PI2 * sp["B"]))
    prof = kc.radial_bvp()
    isos = [None if a is None else nc.rot_pi(a) for a in sp["isos"]]
    phi = None
    for c, R in zip(sp["centers"], isos):
        u = kc.hedgehog(X, Y, Z, prof, center=c, iso=R, taper=(3.0, 3.6))
        phi = u if phi is None else kc.qmul(phi, u)
    return kc.normalize(phi), {}


def measures(phi, h, B):
    w = kc.baryon_density(phi, h).numpy()
    e = k2.static_density(phi, h).numpy()
    E = float(e.sum() * h**3)
    tot, c, M, In = nc.moments(w, h)
    wp = np.clip(w, 0, None)
    _, cE, ME, _ = nc.moments(e, h)
    ev, evec = np.linalg.eigh(M)
    q = nc.multipoles(wp, h, 6)
    pk = nc.peaks(wp, h, 0.5)
    n = w.shape[0]; m = n // 2
    cen = float(w[m - 1:m + 1, m - 1:m + 1, m - 1:m + 1].mean() / w.max())
    rots = dict(Cz90=nc.rot([0, 0, 1], math.pi / 2), Cz45=nc.rot([0, 0, 1], math.pi / 4), Cz30=nc.rot([0, 0, 1], math.pi / 6),
                Cx90=nc.rot([1, 0, 0], math.pi / 2), Cy90=nc.rot([0, 1, 0], math.pi / 2), C2x=nc.rot([1, 0, 0], math.pi),
                C2y=nc.rot([0, 1, 0], math.pi), C2z=nc.rot([0, 0, 1], math.pi), C3=nc.rot([1, 1, 1], 2 * math.pi / 3),
                inv=-np.eye(3), Pz=np.diag([1., 1, -1]))
    asg = {k: nc.asym(wp, h, c, R) for k, R in rots.items()}
    # principal-axis frame: the axis with the most distinct eigenvalue
    gaps = [abs(ev[0] - ev[1]) + abs(ev[0] - ev[2]), 0, abs(ev[2] - ev[0]) + abs(ev[2] - ev[1])]
    ia = int(np.argmax(gaps)); ax = evec[:, ia]
    asp = {f"Cax{d}": nc.asym(wp, h, c, nc.rot(ax, math.radians(d))) for d in (30, 45, 90, 180)}
    Rm = np.eye(3) - 2 * np.outer(ax, ax)
    asp["Pax"] = nc.asym(wp, h, c, Rm)
    # ring profile: density on the circle of maximal density in the plane through the centroid normal to ax
    u = evec[:, (ia + 1) % 3]; v = np.cross(ax, u)
    rr = np.linspace(0, 3, 61); th = np.linspace(0, 2 * np.pi, 72, endpoint=False)
    def samp(P):
        idx = P / h + (n - 1) / 2
        from scipy import ndimage
        return ndimage.map_coordinates(wp, idx, order=3)
    rad = np.array([samp(np.stack([c + r_ * (math.cos(t) * u + math.sin(t) * v) for t in th]).T).mean() for r_ in rr])
    r_ring = float(rr[np.argmax(rad)])
    ring = samp(np.stack([c + r_ring * (math.cos(t) * u + math.sin(t) * v) for t in th]).T)
    return dict(E=E, E_over=E / (PI2 * B), B_meas=tot, centroid=c.tolist(), rms_B=float(math.sqrt(np.trace(M))),
                rms_E=float(math.sqrt(np.trace(ME))), M_eig=ev.tolist(), M_vec=evec.T.tolist(), inertia_eig=np.linalg.eigvalsh(In).tolist(),
                q=q, n_peaks=int(len(pk)), peaks=pk.tolist(), center_over_max=cen, asym_grid=asg, axis=ax.tolist(),
                asym_axis=asp, r_ring=r_ring, ring_var=float((ring.max() - ring.min()) / ring.max()),
                ring_center_over_max=float(rad[0] / rad.max()))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("job"); ap.add_argument("--wall", type=float, default=480); ap.add_argument("--finalize", action="store_true")
    a = ap.parse_args()
    sp = job_spec(a.job)
    logf = open(os.path.join(OUT, "nuclei_logs", f"{a.job}.log"), "a")
    def log(s):
        print(s, flush=True); logf.write(s + "\n"); logf.flush()
    ck = os.path.join(OUT, f"nuclei_ckpt_{a.job}.pt")
    t0 = time.time()
    if os.path.exists(ck):
        st = torch.load(ck, weights_only=False); phi, hist, extra = st["phi"], st["hist"], st["extra"]
        log(f"resume {a.job}: {len(hist)} outer iterations so far, E/12pi^2B = {hist[-1]/(PI2*sp['B']):.6f}")
    else:
        phi, extra = initial(sp, log); hist = []
        m0 = measures(phi, H, sp["B"]); extra["start"] = m0
        log(f"start {a.job}: E/12pi^2B = {m0['E_over']:.6f}  B = {m0['B_meas']:.4f}")
    if st_done(ck):
        log("already done"); return
    if a.finalize:      # budget stop (declared): measure the current, not converged, field
        done = False; MAX = len(hist)
    else:
        MAX = MAXIT
    phi, hist, done = (phi, hist, False) if a.finalize else k2.relax_lbfgs(phi, H, maxit=MAXIT, tol=TOL, log=lambda s: log(s.replace("E/12pi^2", "E/12pi^2(total)")),
                                     every=10, wall=a.wall - (time.time() - t0) - 40, hist=hist)
    fin = done or len(hist) >= MAX
    torch.save(dict(phi=phi, hist=hist, extra=extra, done=fin), ck)
    log(f"chunk end: {len(hist)} its, E/12pi^2B = {hist[-1]/(PI2*sp['B']):.7f}, converged={done}, {time.time()-t0:.0f}s")
    if fin:
        p = phi.detach().requires_grad_(True); E = k2.energy(p, H); (g,) = torch.autograd.grad(E, p)
        gn = float(kc.tangent(phi, g).norm() / H**3 / math.sqrt(phi[0].numel()))
        m = measures(phi, H, sp["B"])
        spj = {k: (v if not isinstance(v, list) else [list(map(float, x)) if isinstance(x, tuple) else x for x in v]) for k, v in sp.items()}
        res = dict(job=a.job, spec=spj, h=H, converged=done, outer_iterations=len(hist), grad_rms=gn, final=m,
                   hist_E_over=[x / (PI2 * sp["B"]) for x in hist[::5]], **extra)
        json.dump(res, open(os.path.join(OUT, f"nuclei_{a.job}.json"), "w"), indent=1)
        np.savez_compressed(os.path.join(OUT, f"nuclei_fields_{a.job}.npz"), b=kc.baryon_density(phi, H).numpy().astype(np.float32),
                            e=k2.static_density(phi, H).numpy().astype(np.float32))
        log(f"FINAL {a.job}: E/12pi^2B = {m['E_over']:.6f}  B = {m['B_meas']:.4f}  rms_B = {m['rms_B']:.3f}  q = {m['q']}")


def st_done(ck):
    return os.path.exists(ck) and torch.load(ck, weights_only=False).get("done", False)


if __name__ == "__main__":
    main()
