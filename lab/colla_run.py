"""colla_run.py JOB [--wall S] - relax one COLLA configuration (resumable, checkpoints), then measure.
JOB = <model>_B<n>[start]   model in {an, ghs0, ghs50, ghs80, ghs95, sx50, sx200} (parameters in colla_units.json)
start: '' = relaxed massless N1 field rm<n> + 3 % smooth noise localised on the knot (seed 7); 'pair' = two hedgehogs d=2.5 (attractive channel);
       'tet' = four hedgehogs on a tetrahedron edge 2.5 (attractive channels), as nuclei_run.
Writes lab/out/colla_JOB.json and lab/out/colla_fields_JOB.npz (baryon density, final field in float32)."""
import math, json, os, sys, time, argparse, numpy as np, torch
import colla_core as cc, knot_core as kc, nuclei_core as nc
import nuclei_run as nr
torch.set_num_threads(int(os.environ.get('COLLA_THREADS', '2')))
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
H = 0.2; MAXIT = 400; TOL = 1e-7   # 1e-8 for an_B1..4 and ghs95_B1; 1e-7 from 22:25 (budget)


def spec(job):
    model, rest = job.split("_B")
    B = int(rest[0]); start = rest[1:]
    U = json.load(open(os.path.join(OUT, "colla_units.json")))
    kw = U["models"][model]
    L = 52 if model.startswith("ghs") else 48
    return dict(model=model, B=B, start=start, kw=kw, L=L)


def initial(sp):
    L = sp["L"]; X, Y, Z = kc.grid((L, L, L), H)
    if sp["start"] == "":
        phi = torch.load(os.path.join(OUT, f"nuclei_ckpt_rm{sp['B']}.pt"), weights_only=False)["phi"]
        if L > phi.shape[1]:
            p = (L - phi.shape[1]) // 2
            vac = torch.tensor([1.0, 0, 0, 0], dtype=phi.dtype).view(4, 1, 1, 1)
            phi = torch.nn.functional.pad((phi - vac)[None], (p, p, p, p, p, p))[0] + vac
        from scipy.ndimage import gaussian_filter
        rng = np.random.default_rng(7)
        nz = np.stack([gaussian_filter(rng.standard_normal(phi.shape[1:]), 2.5) for _ in range(4)])
        nz = torch.tensor(nz / np.abs(nz).max(), dtype=phi.dtype)
        w = (1 - phi[0]) / 2            # smooth noise localised on the knot (no high-frequency, no vacuum noise)
        return kc.normalize(phi + 0.03 * w * nz)
    prof = kc.radial_bvp()
    if sp["start"] == "pair":
        d = 2.5; cen = [(-d / 2, 0, 0), (d / 2, 0, 0)]; isos = [None, "z"]
    elif sp["start"] == "tet":
        a = 2.5 / (2 * math.sqrt(2)); V = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
        cen = [tuple(a * np.array(v)) for v in V]; isos = [None, "x", "y", "z"]
    else:
        raise SystemExit("bad start")
    phi = None
    for c, R in zip(cen, isos):
        u = kc.hedgehog(X, Y, Z, prof, center=c, iso=None if R is None else nc.rot_pi(R), taper=(3.0, 3.6))
        phi = u if phi is None else kc.qmul(phi, u)
    return kc.normalize(phi)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("job"); ap.add_argument("--wall", type=float, default=480)
    a = ap.parse_args()
    sp = spec(a.job); kw = sp["kw"]
    logf = open(os.path.join(OUT, "colla_logs", f"{a.job}.log"), "a")

    def log(s):
        print(s, flush=True); logf.write(s + "\n"); logf.flush()
    ck = os.path.join(OUT, f"colla_ckpt_{a.job}.pt")
    t0 = time.time()
    if os.path.exists(ck):
        st = torch.load(ck, weights_only=False)
        if st.get("done"):
            log("already done"); return
        phi, hist = st["phi"], st["hist"]
        log(f"resume {a.job}: {len(hist)} its, E = {hist[-1]:.6f}")
    else:
        phi = initial(sp); hist = []
        log(f"start {a.job} {kw}: E = {float(cc.energy(phi, H, **kw)):.6f}")
    phi, hist, done = cc.relax(phi, H, maxit=MAXIT, tol=TOL, log=log, every=10, wall=a.wall - (time.time() - t0) - 60,
                               hist=hist, **kw)
    fin = done or len(hist) >= MAXIT
    torch.save(dict(phi=phi, hist=hist, done=fin), ck)
    log(f"chunk end: {len(hist)} its, E = {hist[-1]:.8f}, converged={done}, {time.time()-t0:.0f}s")
    if fin:
        m = nr.measures(phi, H, sp["B"])
        E = float(cc.energy(phi, H, **kw))
        m["E"] = E; m["E_over"] = E / (cc.PI2 * sp["B"]); m["E_over_8pi2"] = E / (8 * math.pi**2)
        res = dict(job=a.job, spec=sp, h=H, converged=done, iterations=len(hist), final=m, hist=hist[::5])
        json.dump(res, open(os.path.join(OUT, f"colla_{a.job}.json"), "w"), indent=1)
        np.savez_compressed(os.path.join(OUT, f"colla_fields_{a.job}.npz"),
                            b=kc.baryon_density(phi, H).numpy().astype(np.float32))
        log(f"FINAL {a.job}: E = {E:.6f}  E/12pi^2B = {m['E_over']:.6f}  B = {m['B_meas']:.4f}  peaks = {m['n_peaks']}  "
            f"q = {[round(m['q'][l], 3) for l in (2, 3, 4)]}  ring c/max = {m['ring_center_over_max']:.3f}")


if __name__ == "__main__":
    main()
