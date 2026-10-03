"""MANY photons: the three "zoom levels" of light, with measured switching criteria.

LEVEL 1 -> 2  shared wave: n photons in the SAME mode (laser = coherent state, also thermal light) are a
              product state φ^⊗n; the n-photon rule reduces to n independent photon balls reading ONE set
              of wave balls (multi_core self-test).  Cost ∝ wave balls + n.  Test: laser pulses on the
              double slit of REGOLA_FOTONE.md (same single-photon rule, photon_rule.py): clicks per pulse
              Poisson, no correlation between photons of a pulse, pattern = single-photon pattern.
LEVEL 2 -> 3  wave only (Maxwell): drop the photon balls, keep the wave; clicks (if wanted) = Poisson
              samples of the flux j_z.  Error of the photon-ball picture vs the wave = shot noise
              1/√(clicks per cell): the wave is the n -> ∞ limit.
LEVEL 3 -> 4  rays: when the aperture is ≫ λ (Fresnel number N_F = a²/(λL) ≫ 1) the wave is not needed:
              a photon ball goes straight except in edge zones of width ~√(λL).  Error vs N_F measured
              (exact Fresnel wave and pilot balls); applied to gamma rays.
"""
import time, sys, os
import numpy as np
from multiprocessing import Pool
from scipy.special import fresnel

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from optics import K, W0, L, inc_amp, in_slit, screen_reference, coarse_corr, chi2, hist, XS, EDGES, BINS, COARSE  # noqa
from photon_rule import lattice_fan, edge_rule, guide, field_lattice, TH_FULL, SLITS, A  # noqa

T0 = time.time()
NPROC = 4
DS_SLIT, DELTA_SLIT = 0.04, 0.2           # REGOLA variant G (fast); final rule uses 0.02 / 0.1


def slit_fan():
    return lattice_fan([c - A / 2 for c in SLITS], [c + A / 2 for c in SLITS], DS_SLIT, inc_amp, TH=TH_FULL, dth=0.01)


def _ds_worker(args):
    x0, seed = args
    fan = slit_fan()
    rng = np.random.default_rng(seed)
    keep, _ = edge_rule(x0, fan, DELTA_SLIT, inc_amp(x0), "E2", rng)
    out = np.full(x0.size, np.nan)
    t = time.time()
    xl, _ = guide(x0[keep], fan, DELTA_SLIT, L, "E2", 2.0, 0.1, 0.01)
    out[keep] = xl
    return out, time.time() - t


def _slit1_worker(args):
    x0, a, Lz, seed = args
    ds = 0.05
    fan = lattice_fan([-a], [a], ds, lambda s: np.ones_like(s), TH=TH_FULL, dth=0.01)
    rng = np.random.default_rng(seed)
    keep, _ = edge_rule(x0, fan, 5 * ds, np.ones_like(x0), "E2", rng)
    out = np.full(x0.size, np.nan)
    xl, _ = guide(x0[keep], fan, 5 * ds, Lz, "E2", 2.0, 0.1, 0.01)
    out[keep] = xl
    return out


def fresnel_slit(x, a, Lz):
    """Exact paraxial (Fresnel) intensity behind a 2-D slit |x|<a lit by a unit plane wave."""
    s = np.sqrt(2 / Lz)                       # λ = 1
    S1, C1 = fresnel(s * (-a - x)); S2, C2 = fresnel(s * (a - x))
    return 0.5 * ((C2 - C1) ** 2 + (S2 - S1) ** 2)


def ray_error(NF, a=1.0):
    """Total-variation distance between the wave pattern and the geometric (ray) pattern, and fraction
    of the light outside the geometric beam, for a slit of half-width a at Fresnel number a²/(λL)."""
    Lz = a * a / NF
    half = a + 12 * np.sqrt(Lz) + 2 * Lz / a      # wide enough for the diffraction wings
    x = np.linspace(-half, half, 400001)
    I = fresnel_slit(x, a, Lz)
    R = (np.abs(x) <= a).astype(float)
    dx = x[1] - x[0]
    tv = 0.5 * np.abs(I - R).sum() * dx / (2 * a)
    out = I[np.abs(x) > a].sum() * dx / (I.sum() * dx)
    return tv, out


def main():
    res = []
    pool = Pool(NPROC)
    # ───────────────── level 1 -> 2: laser pulses on the double slit ─────────────────
    print("ZOOM 1->2  laser pulses (coherent state) on the double slit, one shared set of wave balls", flush=True)
    rng = np.random.default_rng(42)
    mu, npulse = 4.0, 4000
    n = rng.poisson(mu, npulse)
    ne = n.sum()
    xe = rng.normal(0, W0 / 2, ne)                        # source: each photon ball born with |U|²
    ins = in_slit(xe)
    t = time.time()
    parts = np.array_split(np.where(ins)[0], 4 * NPROC)
    outs = pool.map(_ds_worker, [(xe[p], 1000 + i) for i, p in enumerate(parts)])
    xf = np.full(ne, np.nan)
    for p, (o, _) in zip(parts, outs):
        xf[p] = o
    t_guide = time.time() - t
    clicks = xf[np.isfinite(xf)]
    ref = screen_reference("flux")
    corr = coarse_corr(clicks, ref)
    print(f"  {npulse} pulses, μ = {mu}: {ne} photon balls emitted, {ins.sum()} enter the slits, {clicks.size} clicks "
          f"[{t_guide:.0f}s on {NPROC} processes]", flush=True)
    ok = corr > 0.98
    res.append(ok)
    print(f"  pattern of all clicks vs single-photon flux: corr {corr:.4f}, χ²/dof {chi2(clicks, ref):.2f}  "
          f"(REGOLA variant G: 0.9938)  {'PASS' if ok else 'FAIL'}", flush=True)
    pid = np.repeat(np.arange(npulse), n)
    for name, sel in (("anywhere on the screen", np.isfinite(xf)), ("central fringe |x|<5", np.abs(np.nan_to_num(xf, nan=1e9)) < 5),
                      ("right half x>0", np.nan_to_num(xf, nan=-1) > 0)):
        c = np.bincount(pid, sel, npulse)
        F = c.var() / c.mean()
        # error of the Fano factor for Poisson counts: sqrt(2/npulse)
        sdF = np.sqrt(2 / npulse) * np.sqrt(1 + 1 / max(c.mean(), 1e-9))
        okF = abs(F - 1) < 3 * sdF
        res.append(okF)
        print(f"  clicks per pulse, {name:24s}: mean {c.mean():.3f}  Fano var/mean = {F:.3f} ± {sdF:.3f} (Poisson: 1)  "
              f"{'PASS' if okF else 'FAIL'}", flush=True)
    cl = np.bincount(pid, np.isfinite(xf) & (np.nan_to_num(xf) < 0), npulse)
    cr = np.bincount(pid, np.isfinite(xf) & (np.nan_to_num(xf) > 0), npulse)
    g = np.mean(cl * cr) / (cl.mean() * cr.mean())
    blocks = np.array_split(np.arange(npulse), 40)
    sg = np.std([np.mean(cl[b] * cr[b]) / (cl[b].mean() * cr[b].mean()) for b in blocks]) / np.sqrt(40)
    ok = abs(g - 1) < 3 * sg
    res.append(ok)
    print(f"  photons of the same pulse are independent: g2(left, right) = {g:.3f} ± {sg:.3f} (theory 1)  {'PASS' if ok else 'FAIL'}", flush=True)

    # cost: shared wave vs every photon carrying its own wave balls
    fan = slit_fan()
    nb = fan.n_balls
    steps = []
    z = DELTA_SLIT
    while z < L - 1e-9:
        dz = min(2.0, max(0.01, 0.1 * z), L - z); steps.append(dz); z += dz
    th = -fan.TH + (np.arange(fan.nu) + 0.5) * fan.dth
    pos = np.repeat(fan.s, fan.nu); u = np.tile(np.tan(th), fan.s.size)
    t = time.time()
    for dz in steps:
        pos = pos + u * dz                                 # every wave ball flies straight
    t_balls = time.time() - t
    t_read = sum(o[1] for o in outs) / max(clicks.size, 1)
    print(f"  cost: one set of wave balls = {nb} balls ({nb*48/1e6:.1f} MB), flying them to the screen {t_balls*1e3:.0f} ms; "
          f"one photon ball reading them {t_read*1e3:.1f} ms (CPU)", flush=True)
    print("     n photons   shared wave (1 set + n balls)   naive (n sets)          memory shared / naive", flush=True)
    for m in (1, 10, 100, 1000, 10 ** 6):
        sh = t_balls + m * t_read; nv = m * (t_balls + t_read)
        print(f"     {m:>9d}   {sh:12.3g} s                  {nv:12.3g} s        {(nb + m)*48/1e6:9.3g} MB / {m*(nb+1)*48/1e6:9.3g} MB", flush=True)
    # tabulated wave (multi_core): reading cost per photon is O(1), not O(origins)
    from multi_core import System, bs_arm, H
    arm = bs_arm("BS")
    sysm = System([arm], ["BS"], [(1.0, [("a", H, None)])])
    print(f"  with the tabulated wave (multi_core, BS arm, {arm.n_balls} wave balls; tables {sysm.build_time:.1f}s):", flush=True)
    tt = []
    for m in (1000, 4000, 16000):
        X0 = sysm.sample(m, np.random.default_rng(m))
        t = time.time(); sysm.run(X0); tt.append(time.time() - t)
        print(f"     {m:6d} photon balls through the BS: {tt[-1]:6.1f} s  ({tt[-1]/m*1e3:.2f} ms per ball)", flush=True)
    slope = np.polyfit(np.log([1000, 4000, 16000]), np.log(tt), 1)[0]
    print(f"     time ∝ n^{slope:.2f} (fixed table cost {sysm.build_time:.1f}s paid once)", flush=True)

    # ───────────────── level 2 -> 3: drop the photon balls ─────────────────
    print("ZOOM 2->3  bright beam: wave only (Maxwell limit)", flush=True)
    xs_f = np.linspace(-60, 60, BINS * 8 + 1); xs_f = 0.5 * (xs_f[1:] + xs_f[:-1])
    U, Ux, Uz = field_lattice(xs_f, L, fan)
    jz = np.imag(np.conj(U) * Uz).sum(-1)
    wave = jz.reshape(BINS, 8).sum(1); wave /= wave.sum()
    bright = ref.reshape(-1, COARSE).sum(1) > 0.01
    e_wave = np.abs(wave.reshape(-1, COARSE).sum(1)[bright] / ref.reshape(-1, COARSE).sum(1)[bright] - 1)
    ok = e_wave.max() < 0.02
    res.append(ok)
    print(f"  intensity from the wave balls alone (no photon balls) vs exact Rayleigh–Sommerfeld flux: "
          f"max error {100*e_wave.max():.2f}% on the bright coarse bins (discretization ds = {DS_SLIT})  {'PASS' if ok else 'FAIL'}", flush=True)
    np.save(os.path.join(HERE, "out", "multi_beams_clicks.npy"), clicks)
    rr = np.random.default_rng(5)
    ns, errs = (100, 300, 1000), []
    pc = ref.reshape(-1, COARSE).sum(1)
    rel = lambda h, m: np.sqrt(np.mean((h[bright] / (m * pc[bright]) - 1) ** 2))
    for m in ns:
        errs.append(np.mean([rel(hist(rr.choice(clicks, m, replace=False)).reshape(-1, COARSE).sum(1), m) for _ in range(60)]))
    shot = lambda m: np.sqrt(np.mean((1 - pc[bright]) / (m * pc[bright])))
    e_all = rel(hist(clicks).reshape(-1, COARSE).sum(1), clicks.size)
    floor = np.sqrt(max(e_all ** 2 - shot(clicks.size) ** 2, 0))
    print(f"  photon-ball clicks vs the wave: rms relative error per bright bin (8 λ wide, {bright.sum()} bins):", flush=True)
    for m, e in zip(ns + (clicks.size,), errs + [e_all]):
        print(f"     {m:5d} clicks: {100*e:5.1f}%   pure shot noise {100*shot(m):5.1f}%   shot ⊕ floor {100*np.hypot(shot(m), floor):5.1f}%", flush=True)
    sl = np.polyfit(np.log(ns + (clicks.size,)), np.log(errs + [e_all]), 1)
    ok = all(abs(e / shot(m) - 1) < 0.15 for m, e in zip(ns[:2], errs[:2]))
    res.append(ok)
    print(f"  few clicks (≤ 300): error = shot noise within 15%  {'PASS' if ok else 'FAIL'};  naive power law n^{sl[0]:.2f} "
          f"(not −0.5: the photon balls of this fast variant carry a systematic floor of {100*floor:.1f}% per bin vs the exact "
          f"wave; the wave alone has {100*e_wave.max():.2f}%)", flush=True)
    for metric in ("fringe windows",):
        from optics import fringe_errors
        e, r, mx = fringe_errors(clicks, ref)
        print(f"  same clicks with the REGOLA metric (bright-fringe windows): rms {100*r:.1f}%, max {100*mx:.1f}%", flush=True)
    need = 1 / (0.01 ** 2 * pc[bright].min())
    print(f"  shot noise alone reaches 1% in every bright bin at n ≈ {need:.2g} clicks; a bright beam (n̄ ≫ that) is "
          f"described by the wave alone with no loss", flush=True)
    # semiclassical clicks from the wave alone are statistically the same as the photon balls' clicks
    cdf = np.cumsum(wave); cdf /= cdf[-1]
    samp = XS[np.searchsorted(cdf, rr.random(clicks.size))] + (rr.random(clicks.size) - 0.5) * (EDGES[1] - EDGES[0])
    h1 = hist(clicks).reshape(-1, COARSE).sum(1); h2 = hist(samp).reshape(-1, COARSE).sum(1)
    m = (h1 + h2) > 10
    c2 = (((h1 - h2) ** 2) / (h1 + h2))[m].sum() / (m.sum() - 1)
    ok = c2 < 1 + 3 * np.sqrt(2 / (m.sum() - 1))
    res.append(ok)
    print(f"  clicks drawn from the wave alone (Poisson sampling of j_z) vs photon-ball clicks: two-sample χ²/dof = {c2:.2f} "
          f"(same distribution: 1 ± {np.sqrt(2/(m.sum()-1)):.2f})  {'PASS' if ok else 'FAIL'}", flush=True)
    print("  criterion: the photon balls are needed only for NON-classical light (single photons, pairs, entangled,"
          " sub-Poisson); for classical light the mean number of photons per mode n̄ decides only the noise:", flush=True)
    hnu = 6.626e-34 * 3e8 / 633e-9
    print(f"     He-Ne laser 1 mW, linewidth 1 MHz: n̄ = (P/hν)·τc = {1e-3/hnu*1e-6:.1e} photons per mode -> wave only", flush=True)
    x = 6.626e-34 * 3e8 / (500e-9 * 1.381e-23 * 5800)
    print(f"     sunlight at 500 nm (T = 5800 K): n̄ = 1/(e^(hν/kT)−1) = {1/np.expm1(x):.3f} per mode -> clicks are"
          f" Poisson draws of a fluctuating wave (thermal, level 2/3), never needs pair-level balls", flush=True)

    # ───────────────── level 3 -> 4: rays ─────────────────
    print("ZOOM 3->4  ray (geometric-optics) limit, slit of half-width a at distance L, N_F = a²/(λL)", flush=True)
    NFs = (0.25, 1, 4, 16, 64, 256, 1024, 4096, 16384)
    tvs = []
    for NF in NFs:
        tv, outf = ray_error(NF)
        tvs.append(tv)
        print(f"     N_F = {NF:7g}:  wave vs rays, total-variation distance {100*tv:7.3f}%   light outside the geometric beam {100*outf:7.3f}%", flush=True)
    outs_ = [ray_error(NF)[1] for NF in NFs]
    sl = np.polyfit(np.log(NFs[2:]), np.log(outs_[2:]), 1)
    ok = abs(sl[0] + 0.5) < 0.05
    res.append(ok)
    print(f"  light outside the geometric beam ∝ N_F^{sl[0]:.3f} (expected −1/2: edge zone √(λL) wide vs beam 2a)  {'PASS' if ok else 'FAIL'}", flush=True)
    sl2 = np.polyfit(np.log(NFs[2:]), np.log(tvs[2:]), 1)
    Xf = np.column_stack([np.ones(len(NFs) - 2), np.log(NFs[2:])]) / np.sqrt(np.array(NFs[2:]))[:, None]
    cfit = np.linalg.lstsq(Xf, np.array(tvs[2:]), rcond=None)[0]
    tvfit = lambda nf: (cfit[0] + cfit[1] * np.log(nf)) / np.sqrt(nf)
    worst = max(abs(tvfit(nf) / tv - 1) for nf, tv in zip(NFs[2:], tvs[2:]))
    ok = worst < 0.05
    res.append(ok)
    print(f"  total-variation error: naive power law N_F^{sl2[0]:.3f}; the Fresnel ripples inside the beam decay only as "
          f"1/distance from the edge, so TV = ({cfit[0]:.3f} + {cfit[1]:.3f}·ln N_F)/√N_F  (fit within {100*worst:.1f}%)  {'PASS' if ok else 'FAIL'}", flush=True)
    nf1 = 10 ** np.interp(-2, np.log10([tvfit(10 ** q) for q in np.linspace(5, 1, 400)]), np.linspace(5, 1, 400))
    print(f"  rule of thumb: rays are good to 1% (TV) when N_F > {nf1:.0f}; to 10% when N_F > "
          f"{10 ** np.interp(-1, np.log10([tvfit(10 ** q) for q in np.linspace(5, 0, 500)]), np.linspace(5, 0, 500)):.0f}", flush=True)
    # pilot balls through a single slit: they go straight except near the edges
    a = 4.0
    print(f"  pilot photon balls through a slit of half-width {a:g} λ (same single-photon rule):", flush=True)
    rows = []
    jobs = []
    for i, NF in enumerate((1, 4, 16)):
        Lz = a * a / NF
        x0 = np.random.default_rng(70 + i).uniform(-a, a, 1500)
        for j, part in enumerate(np.array_split(x0, NPROC)):
            jobs.append((part, a, Lz, 700 + 10 * i + j))
    t = time.time()
    outs1 = pool.map(_slit1_worker, jobs)
    for i, NF in enumerate((1, 4, 16)):
        xf1 = np.concatenate(outs1[i * NPROC:(i + 1) * NPROC])
        x0 = np.concatenate([jb[0] for jb in jobs[i * NPROC:(i + 1) * NPROC]])
        ok_ = np.isfinite(xf1)
        dev = np.abs(xf1[ok_] - x0[ok_])
        tv_ball = 0.5 * np.abs(np.histogram(xf1[ok_], np.linspace(-3 * a, 3 * a, 61))[0] / ok_.sum()
                               - np.histogram(np.random.default_rng(1).uniform(-a, a, 200000), np.linspace(-3 * a, 3 * a, 61))[0] / 200000).sum()
        tvw, _ = ray_error(NF)
        rows.append((NF, (dev > 0.1 * a).mean(), np.median(dev) / a, tv_ball, tvw))
        print(f"     N_F = {NF:3d} (L = {a*a/NF:5.2f} λ): balls deviating > 0.1a from the straight ray {100*(dev > 0.1*a).mean():5.1f}%, "
              f"median deviation {np.median(dev)/a:.3f} a;  click histogram vs rays TV {100*tv_ball:.1f}% (exact wave {100*tvw:.1f}%)", flush=True)
    ok = rows[-1][1] < rows[0][1] and rows[-1][3] < rows[0][3]
    res.append(ok)
    print(f"  deviations shrink as N_F grows {'PASS' if ok else 'FAIL'}  [{time.time()-t:.0f}s]", flush=True)
    pool.close()

    print("GAMMA RAYS: λ = 1e-12 m, aperture 1 mm (a = 0.5 mm), screen at L = 1 m", flush=True)
    lam, a_g, Lg = 1e-12, 0.5e-3, 1.0
    NF = a_g ** 2 / (lam * Lg)
    err, outg = ray_error(NF, a=1.0)
    orig2d = 2 * a_g / (0.1 * lam)
    orig3d = np.pi * a_g ** 2 / (0.1 * lam) ** 2
    print(f"  N_F = {NF:.3g};  ray error (exact Fresnel) TV {100*err:.3f}%, light outside the geometric beam {100*outg:.4f}% (edge zone √(λL) = {np.sqrt(lam*Lg)*1e6:.1f} µm vs a = 500 µm)", flush=True)
    print(f"  wave balls with the pilot rule (origins every 0.1 λ, ≥ 31 directions each):", flush=True)
    print(f"     2-D: {orig2d:.1e} origins -> {orig2d*31:.1e} balls;   3-D: {orig3d:.1e} origins -> {orig3d*31:.1e} balls "
          f"(× 48 bytes = {orig3d*31*48:.1e} B);  photon-ball reading cost {orig3d:.1e} wavelets per step", flush=True)
    print(f"     even only the edge zone in 3-D: {2*np.pi*a_g*3*np.sqrt(lam*Lg)/(0.1*lam)**2:.1e} origins", flush=True)
    print(f"  ray limit: 1 ray per photon, cost O(1) per photon per object -> mandatory for gamma rays", flush=True)
    print(f"RESULT: {sum(res)}/{len(res)} PASS   total {time.time()-T0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
