"""ZOOM test a: 1-D chain of N balls with springs (mass m, stiffness K, spacing a), periodic.

Exact answers: sound speed c = a sqrt(K/m); dispersion w(q) = 2 sqrt(K/m) |sin(q a / 2)|; a block of k^n
balls has mass k^n m; springs in series give long-wavelength stiffness K / k^n.  The best possible
coarse rule for the centres of mass (potential of mean force, PMF) is computed exactly from the Gaussian
covariance and printed next to the learned one.

Pipeline (zoom_core, NO code changes between levels):
  level 0  brute-force Langevin MD of the fine chain (thermal samples + fine forces)
  level 1  block(k=2) -> coarse forces = sum of member forces -> learn('force')
  level 2  Langevin MD of the LEARNED level-1 rule -> block(k=2) -> learn   (true repetition: the level-2
           rule never sees the fine chain)          ... level 3 likewise.
Checks: learned couplings vs exact PMF, K_eff vs K/k^n, sound speed (linearised and measured by NVE MD
of a standing wave), dispersion up to the coarse Nyquist q = pi/(k^n a), energy conservation (NVE),
spot-check of the bond-length variance in random regions, cost of MD per unit physical time.
Run: timeout 590 python lab/zoom_chain.py
"""
import os, sys, time, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zoom_core as zc

T0 = time.time()
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)
N0, M0, K0, A0, TEMP, KB = 1024, 1.0, 1.0, 1.0, 0.01, 2
LEVELS = 3
REPL, NSAMP = 8, 150
rng = np.random.default_rng(20261009)


def log(*a):
    print(f"[{time.time() - T0:6.1f}s]", *a, flush=True)


# ------------------------------------------------------------------ fine (brute force) system
BOX = N0 * A0


def fine_force(st):
    x = st.pos[:, 0]
    dr = zc.min_image(np.roll(x, -1) - x, BOX)
    dl = zc.min_image(x - np.roll(x, 1), BOX)
    return (K0 * (dr - A0) - K0 * (dl - A0))[:, None]


def fine_equilibrium(rng):
    """Exact thermal sample of the periodic harmonic chain (used only as the MD starting point)."""
    q = 2 * np.pi * np.fft.rfftfreq(N0)
    w2 = (2 - 2 * np.cos(q)) * K0 / M0
    z = rng.normal(size=q.size) + 1j * rng.normal(size=q.size)
    z[0] = 0
    z[1:] *= np.sqrt(TEMP / (M0 * w2[1:]) / 2 * N0)
    z[-1] = z[-1].real * np.sqrt(2)
    u = np.fft.irfft(z, n=N0)
    ref = (np.arange(N0) + 0.5) * A0
    v = rng.normal(size=N0) * np.sqrt(TEMP / M0)
    v -= v.mean()
    return zc.State(((ref + u) % BOX)[:, None], np.full(N0, M0), A0, box=[BOX], ref=ref[:, None],
                    vel=v[:, None], mass=np.full(N0, M0))


def omega_max(model, st):
    """Generic stability estimate from the learned pair stiffness: w_max^2 = 4 sum|kappa| / M."""
    u = np.arange(1, int(model.rcut) + 1)
    kap = -model.dphi(u.astype(float)) / st.h
    return np.sqrt(4 * np.abs(kap).sum() / st.mass.mean())


def sample(force_or_model, st0_list, dt, rng, nsteps, every):
    out = []
    for st0 in st0_list:
        zc.langevin(force_or_model, st0, dt, nsteps, TEMP, 0.5 / (dt * 10), rng, every=every, out=out)
    return out


# ------------------------------------------------------------------ exact PMF of block centres
def exact_pmf(kn):
    """Exact coarse couplings kappa_n (PMF Hessian of block centres) for blocks of kn fine balls."""
    D = 2 * np.eye(N0) - np.roll(np.eye(N0), 1, 0) - np.roll(np.eye(N0), -1, 0)
    C = np.linalg.pinv(D * K0)                    # covariance / T
    Bm = np.zeros((N0 // kn, N0))
    for I in range(N0 // kn):
        Bm[I, I * kn:(I + 1) * kn] = 1 / kn
    Hc = np.linalg.pinv(Bm @ C @ Bm.T)
    row = Hc[0]
    return np.array([-row[n] for n in range(1, N0 // kn // 2)])


def disp(kap, M, D, q):
    n = np.arange(1, kap.size + 1)
    return np.sqrt(np.maximum(2 / M * (kap[None, :] * (1 - np.cos(np.outer(q, n) * D))).sum(1), 0))


# ------------------------------------------------------------------ MD measurements
def standing_wave(model_force, st_lat, dt, tmax):
    """NVE Verlet from a small standing wave of the longest wavelength; returns (w, rel energy err, wall)."""
    st = st_lat.copy()
    D = st.h
    L = BOX
    q = 2 * np.pi / L
    X0 = st.ref[:, 0]
    amp = 0.01 * D
    st.pos = ((X0 + amp * np.sin(q * X0)) % L)[:, None]
    st.vel = np.zeros_like(st.pos)
    m = st.mass[:, None]
    proj = []
    t0 = time.time()
    F = model_force(st)
    nst = int(tmax / dt)
    for t in range(nst):
        st.vel += 0.5 * dt * F / m
        st.pos += dt * st.vel
        F = model_force(st)
        st.vel += 0.5 * dt * F / m
        u = zc.min_image(st.pos[:, 0] - X0, L)
        proj.append((u * np.sin(q * X0)).sum())
    wall = time.time() - t0
    p = np.array(proj)
    tt = dt * (np.arange(nst) + 1)
    s = np.sign(p)
    zc_idx = np.where(s[:-1] * s[1:] < 0)[0]
    tz = tt[zc_idx] - p[zc_idx] * dt / (p[zc_idx + 1] - p[zc_idx])
    halfp = np.diff(tz).mean()
    return np.pi / halfp, wall, nst


def nve_energy(model, st, dt, nsteps):
    st = st.copy()
    m = st.mass[:, None]
    E = []
    F = model.force(st)
    for t in range(nsteps):
        st.vel += 0.5 * dt * F / m
        st.pos += dt * st.vel
        F = model.force(st)
        st.vel += 0.5 * dt * F / m
        if t % 10 == 0:
            E.append(model.energy(st) + 0.5 * (m * st.vel ** 2).sum())
    E = np.array(E)
    kin = 0.5 * TEMP * st.N
    return float(np.abs(E - E[0]).max() / kin), float((E[-1] - E[0]) / kin)


def bond_var(states, idx):
    v = []
    for st in states:
        x = st.pos[:, 0]
        b = zc.min_image(np.roll(x, -1) - x, BOX)[idx]
        v.append(b)
    return np.var(np.concatenate(v))


# ------------------------------------------------------------------ main
def main():
    res = dict(levels=[])
    # level 0: brute force
    starts = [fine_equilibrium(rng) for _ in range(REPL)]
    log(f"fine chain N={N0}, m=K=a=1, T={TEMP}, {REPL} replicas x {NSAMP} samples (Langevin dt=0.1)")
    data = sample(fine_force, starts, 0.1, rng, NSAMP * 20, 20)
    states = [d[0] for d in data]; forces = [d[1] for d in data]
    log(f"fine samples: {len(states)}")
    # sanity: learner on the fine level itself (no blocking) must recover K = 1 nearest neighbour only
    m0 = zc.LocalModel()
    info = m0.learn_force(states[:40], forces[:40])
    kap0 = -m0.dphi(np.arange(1, 7.0)) / A0
    log(f"level 0 self-learn: kappa_1..6 = {np.round(kap0, 4)} (exact 1,0,0,...) rel.residual {info['rel_residual']:.2e}")
    res["level0_kappa"] = kap0.tolist()
    models = [m0]
    lat0 = zc.State(starts[0].ref.copy(), np.full(N0, M0), A0, box=[BOX], ref=starts[0].ref, vel=np.zeros((N0, 1)), mass=np.full(N0, M0))
    lattices = [lat0]
    blk_states, blk_forces = states, forces
    fine_blocked = states            # fine samples, blocked repeatedly for the direct comparison
    for lev in range(1, LEVELS + 1):
        kn = KB ** lev
        # block the samples of the level below (fine for level 1, LEARNED model for levels >= 2)
        cs, cf = [], []
        B = None
        for st, F in zip(blk_states, blk_forces):
            c, B = zc.block(st, KB, blocking=B)
            cs.append(c); cf.append(B.sum(F))
        model = zc.LocalModel()
        info = model.learn_force(cs, cf)
        D = cs[0].h; Mb = cs[0].mass.mean()
        u = np.arange(1, 7.0)
        kap = -model.dphi(u) / D
        kex = exact_pmf(kn)
        Keff = (u ** 2 * kap).sum()
        Kex_inf = (np.arange(1, kex.size + 1) ** 2 * kex).sum()
        c_lin = D * np.sqrt(Keff / Mb)
        # direct route: fine samples blocked lev times, learned at once (to see error accumulation)
        dsx, dfx = [], []
        Bd = None
        for st, F in zip(states, forces):
            c = st; FF = F
            for _ in range(lev):
                c, Bt = zc.block(c, KB); FF = Bt.sum(FF)
            dsx.append(c); dfx.append(FF)
        mdir = zc.LocalModel(); mdir.learn_force(dsx, dfx)
        kapd = -mdir.dphi(u) / D
        Keffd = (u ** 2 * kapd).sum()
        # dispersion
        qn = np.pi / D
        qs = qn * np.array([1 / 16, 1 / 8, 1 / 4, 1 / 2, 3 / 4, 1.0])
        w_ex = 2 * np.sqrt(K0 / M0) * np.abs(np.sin(qs * A0 / 2))
        w_l = disp(kap, Mb, D, qs)
        w_p = disp(kex, Mb, D, qs)
        qf = np.linspace(qn / 400, qn, 400)
        err_f = np.abs(disp(kap, Mb, D, qf) / (2 * np.abs(np.sin(qf / 2))) - 1)
        q5 = qf[np.argmax(err_f > 0.05)] / qn if (err_f > 0.05).any() else None
        q1 = qf[np.argmax(err_f > 0.01)] / qn if (err_f > 0.01).any() else None
        errp = np.abs(disp(kex, Mb, D, qf) / (2 * np.abs(np.sin(qf / 2))) - 1)
        q5p = qf[np.argmax(errp > 0.05)] / qn if (errp > 0.05).any() else None
        # MD: thermal samples of THIS learned level (used to learn the next level), NVE energy check
        wmax = omega_max(model, cs[0])
        dt = 0.2 / wmax
        lat = zc.State(cs[0].ref.copy(), cs[0].q, D, box=[BOX], n=cs[0].n, ref=cs[0].ref, vel=np.zeros_like(cs[0].pos), mass=cs[0].mass)
        e_max, e_drift = nve_energy(model, cs[0], dt, 2000)
        starts_c = cs[::NSAMP][:REPL]
        every = max(int(round(5.0 / dt)), 1)
        t_s = time.time()
        smp = sample(model, starts_c, dt, rng, every * NSAMP, every)
        t_s = time.time() - t_s
        # spot checks: bond-length variance in random regions, coarse model vs blocked fine samples
        fb = []
        for st in fine_blocked:
            fb.append(zc.block(st, KB)[0])
        fine_blocked = fb
        checks, nulls = [], []
        nI = cs[0].N
        half = len(fine_blocked) // 2
        for r in range(6):
            i0 = rng.integers(nI)
            region = (np.arange(32) + i0) % nI
            ck = zc.check(region, lambda R: bond_var(fine_blocked, R), lambda R: bond_var([s for s, _ in smp], R),
                          lambda a, b: abs(b / a - 1), 0.05)
            checks.append(ck)
            nulls.append(abs(bond_var(fine_blocked[:half], region) / bond_var(fine_blocked[half:], region) - 1))
        noise = float(np.mean(nulls))      # same check, fine data against fine data (two independent halves)
        L = dict(level=lev, k_total=kn, N=nI, mass=float(Mb), mass_exact=kn * M0, spacing=D,
                 kappa_learned=kap.tolist(), kappa_exact_pmf=kex[:6].tolist(), kappa_direct=kapd.tolist(),
                 Keff=float(Keff), Keff_direct=float(Keffd), Keff_exact=K0 / kn, Keff_pmf_full=float(Kex_inf),
                 c_lin=float(c_lin), q=list(qs / qn), w_exact=w_ex.tolist(), w_learned=w_l.tolist(), w_pmf=w_p.tolist(),
                 q_err5=q5, q_err1=q1, q_err5_pmf=q5p, dt=dt, wmax=float(wmax), nve_maxdev=e_max, nve_drift=e_drift,
                 eint_mean=float(np.mean([s.eint.mean() for s in cs])), eint_equip=(kn - 1) * TEMP / 2,
                 checks=[dict(err=c["err"], flagged=c["flagged"]) for c in checks], check_noise=float(noise),
                 rel_residual=info["rel_residual"], sample_time=t_s)
        res["levels"].append(L)
        log(f"LEVEL {lev}: {nI} super-balls of {kn} balls, mass {Mb:.3f} (exact {kn}), spacing {D}")
        log(f"   kappa_1..6 learned {np.round(kap, 4)}")
        log(f"   kappa_1..6 exactPMF {np.round(kex[:6], 4)}   direct-from-fine {np.round(kapd, 4)}")
        log(f"   K_eff = sum n^2 kappa_n: learned {Keff:.4f}  direct {Keffd:.4f}  exact K/k^n = {K0 / kn:.4f}  (PMF to n<=6: {(np.arange(1, 7) ** 2 * kex[:6]).sum():.4f})")
        log(f"   sound speed (linearised) {c_lin:.4f}  exact {A0 * np.sqrt(K0 / M0):.4f}  err {abs(c_lin - 1) * 100:.2f}%")
        log(f"   dispersion q/q_Nyq {np.round(qs / qn, 3)}")
        log(f"      exact fine  {np.round(w_ex, 4)}")
        log(f"      learned     {np.round(w_l, 4)}")
        log(f"      exact PMF   {np.round(w_p, 4)}")
        log(f"   error >1% from q/q_Nyq = {q1}, >5% from {q5} (exact-PMF rule: >5% from {q5p})")
        log(f"   NVE dt={dt:.3f} (w_max {wmax:.3f}): max|dE|/(N T/2) = {e_max:.2e}, drift {e_drift:.2e}")
        log(f"   internal energy per super-ball {L['eint_mean']:.5f} (equipartition (k^n-1)T/2 = {L['eint_equip']:.5f}) -> heat lost to the coarse rule")
        log(f"   spot checks (bond variance in 32 super-balls, tol 5%; fine-vs-fine noise {noise * 100:.1f}%): " +
            ", ".join(f"{c['err'] * 100:.1f}%{'!' if c['flagged'] else ''}" for c in checks))
        blk_states = [s for s, _ in smp]; blk_forces = [f for _, f in smp]
        models.append(model); lattices.append(lat)
    # --------------------------------------------------------- sound speed by MD + cost per level
    log("standing wave (longest wavelength), NVE, same physical duration 2.2 periods at every level")
    rows = []
    for lev, (model, lat) in enumerate(zip(models, lattices)):
        st = lat
        wmax = omega_max(model, st)
        dt = 0.2 / wmax
        w, wall, nst = standing_wave(model.force, st, dt, 2.2 * BOX / 1.0)
        qmin = 2 * np.pi / BOX
        c = w / qmin
        w_ex = 2 * np.sin(qmin / 2)
        rows.append(dict(level=lev, N=st.N, dt=dt, steps=nst, wall=wall, c=c, err=abs(w / w_ex - 1)))
        log(f"   level {lev}: N={st.N:5d} dt={dt:.3f} steps={nst:6d} wall {wall:6.2f}s  c={c:.5f} (exact sound speed at this q {w_ex / qmin:.5f}) err {abs(w / w_ex - 1) * 100:.3f}%")
    res["md"] = rows
    json.dump(res, open(os.path.join(OUT, "zoom_chain.json"), "w"), indent=1)
    # figure
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.4), dpi=90)
        qf = np.linspace(1e-3, np.pi, 400)
        ax[0].plot(qf, 2 * np.sin(qf / 2), "k", lw=2, label="exact fine chain")
        for L, col in zip(res["levels"], ["C0", "C1", "C2"]):
            D = L["spacing"]; q = np.linspace(1e-3, np.pi / D, 200)
            ax[0].plot(q, disp(np.array(L["kappa_learned"]), L["mass"], D, q), col, label=f"level {L['level']} (learned)")
            ax[0].axvline(np.pi / D, color=col, ls=":", lw=0.8)
        ax[0].set_xlabel("q a"); ax[0].set_ylabel("ω √(m/K)"); ax[0].legend(fontsize=7); ax[0].set_title("dispersion")
        for L, col in zip(res["levels"], ["C0", "C1", "C2"]):
            D = L["spacing"]; q = np.linspace(np.pi / D / 200, np.pi / D, 200)
            e = disp(np.array(L["kappa_learned"]), L["mass"], D, q) / (2 * np.sin(q / 2)) - 1
            ax[1].plot(q * D / np.pi, 100 * e, col, label=f"level {L['level']}")
        ax[1].set_xlabel("q / q_Nyquist(coarse)"); ax[1].set_ylabel("error ω (%)"); ax[1].axhline(-5, color="gray", lw=0.5); ax[1].legend(fontsize=7)
        ax[1].set_title("where the coarse rule breaks")
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "zoom_chain_dispersion.png")); plt.close(fig)
    except Exception as e:
        log("figure failed:", e)
    log("done")


if __name__ == "__main__":
    main()
