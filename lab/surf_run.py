"""surf_run.py - production runs of the study SURFISTI_ESATTI (exams in lab/out/surf_prereg.json).
usage: python surf_run.py f6 | cont | h1s | beat | a3d <1s|1s2p|free> | orb
Every run writes lab/out/surf_<run>.json; the 3-D runs are checkpointed (exit code 3 = run again).
"""
import json
import pickle
import sys
import time
import numpy as np
from surf_core import *

CK = os.path.join(OUT, "surf_ckpt")
os.makedirs(CK, exist_ok=True)
SMOKE = os.environ.get("SURF_SMOKE") == "1"   # code smoke test only (tiny runs, outputs *_smoke)
PRE = json.load(open(os.path.join(OUT, "surf_prereg.json"))) if os.path.exists(os.path.join(OUT, "surf_prereg.json")) else None


def save(name, res):
    if SMOKE:
        name = name + "_smoke"
    with open(os.path.join(OUT, f"surf_{name}.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    print("written", f"surf_{name}.json", flush=True)


def quantile_of(cell, xi, rho):
    """exact CDF of the piecewise-uniform measure at the surfer positions (1-D)"""
    C = np.concatenate([[0.0], np.cumsum(rho)])
    return (C[cell[:, 0]] + xi[:, 0] * rho[cell[:, 0]]) / C[-1]


def quantile_bins_H(q, nb=20):
    """coarse H = sum p ln(p / (1/nb)) of the surfers' quantiles in nb equal-probability bins"""
    c = np.bincount(np.minimum((q * nb).astype(int), nb - 1), minlength=nb) / len(q)
    m = c > 0
    return float(np.sum(c[m] * np.log(c[m] * nb)))


# ======================================================================================================
#  E1  the old F6 test (1-D free walk, mixed packet) with the exact flow
# ======================================================================================================
def run_f6():
    rng = np.random.default_rng(SEED)
    th = 0.1
    a = th
    n = int(round(400.0 / a))
    x = (np.arange(n) - n / 2 + 0.5) * a
    g = np.exp(-x**2 / (4 * 5.0**2)).astype(complex)
    psi = np.array([g, 0 * g])
    psi /= np.linalg.norm(psi)
    W = Walk1D(n, th)
    Ns = 2000
    rho = density(psi)
    cell, xi = born_sample(rng, rho, Ns)
    q0 = quantile_of(cell, xi, rho)
    # old reader (zz_free F6, verbatim logic): RK4 over two ticks, linear interpolation, same initial points
    pos = x[cell[:, 0]] + (xi[:, 0] - 0.5) * a
    def vel(ps, p):
        R = np.interp(p, x, ps[0].real) + 1j * np.interp(p, x, ps[0].imag)
        L = np.interp(p, x, ps[1].real) + 1j * np.interp(p, x, ps[1].imag)
        return (np.abs(R)**2 - np.abs(L)**2) / (np.abs(R)**2 + np.abs(L)**2)
    def qold(ps, p):        # the old quantile estimator (zz_free)
        r = density(ps); c = np.cumsum(r) - 0.5 * r
        return np.interp(p, x, c / r.sum())
    qo0 = qold(psi, pos)
    ticks = int(round(40.0 / a)) if not SMOKE else 10
    err_new, err_old, err_old_exactq = 0.0, 0.0, 0.0
    vmax_new = 0.0
    W2 = Walk1D(n, th)
    for t in range(0, ticks, 2):
        p_start = psi
        for _ in range(2):
            r0 = density(psi)
            psi, subs = W.tick(psi)
            r1 = density(psi)
            F = net_flux(0, subs[0][1], subs[0][2])
            c_prev = cell.copy(); x_prev = xi.copy()
            flow_tick(cell, xi, r0, r1, [F], periodic=True)
            disp = (cell[:, 0] + xi[:, 0]) - (c_prev[:, 0] + x_prev[:, 0])
            vmax_new = max(vmax_new, float(np.max(np.abs(disp))))
        p1 = W2.tick(p_start)[0]; p2 = W2.tick(p1)[0]
        dt = 2 * a
        k1 = vel(p_start, pos); k2 = vel(p1, pos + 0.5 * dt * k1); k3 = vel(p1, pos + 0.5 * dt * k2); k4 = vel(p2, pos + dt * k3)
        pos = pos + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        q = quantile_of(cell, xi, density(psi))
        err_new = max(err_new, float(np.max(np.abs(q - q0))))
        err_old = max(err_old, float(np.max(np.abs(qold(p2, pos) - qo0))))
    res = dict(theta=th, n=n, surfers=Ns, T=40.0, quantile_err_flow=err_new, quantile_err_old_reader=err_old,
               max_cells_per_tick_flow=vmax_new, note="flow = S4 (exact); old = zz_free F6 reader (linear interp + RK4)")
    print(res)
    save("f6", res)


# ======================================================================================================
#  E2  continuum limit in 1-D: lattice flow vs continuum de Broglie-Bohm-Dirac (spectral)
# ======================================================================================================
def cont_density(t, L=400.0, h=0.025, sigma=5.0):
    """continuum 1-D Dirac (H = sigma_z p + sigma_x, units hbar = c = m = 1), R-Gaussian at t = 0, exact."""
    M = int(round(L / h))
    xc = (np.arange(M) - M / 2 + 0.5) * h
    g = np.exp(-xc**2 / (4 * sigma**2)).astype(complex)
    k = 2 * np.pi * np.fft.fftfreq(M, d=h)
    Rk, Lk = np.fft.fft(g), np.zeros(M, complex)
    E = np.sqrt(k**2 + 1)
    c, s = np.cos(E * t), np.sin(E * t) / E
    R2 = c * Rk - 1j * s * (k * Rk + Lk)
    L2 = c * Lk - 1j * s * (Rk - k * Lk)
    R, Lw = np.fft.ifft(R2), np.fft.ifft(L2)
    return xc, R, Lw


def quantile_function_cont(xc, rho, u):
    h = xc[1] - xc[0]
    C = np.concatenate([[0], np.cumsum(rho)]) ; C = C / C[-1]
    edges = np.concatenate([xc - h / 2, [xc[-1] + h / 2]])
    return np.interp(u, C, edges)


def quantile_function_latt(a, rho, u):
    n = len(rho)
    edges = (np.arange(n + 1) - n / 2) * a
    C = np.concatenate([[0], np.cumsum(rho)]); C = C / C[-1]
    return np.interp(u, C, edges)


def run_cont():
    rng = np.random.default_rng(SEED + 1)
    u = np.linspace(0.01, 0.99, 99)
    out = {}
    times = np.arange(0, 41, 4.0) if not SMOKE else np.array([0.0, 1.0])
    for a in (0.2, 0.1, 0.05):
        n = int(round(400.0 / a))
        x = (np.arange(n) - n / 2 + 0.5) * a
        g = np.exp(-x**2 / (4 * 25.0)).astype(complex)
        psi = np.array([g, 0 * g]); psi /= np.linalg.norm(psi)
        W = Walk1D(n, a)
        Ns = 4000
        cell, xi = born_sample(rng, density(psi), Ns)
        cb = cell.copy()                         # Bell (S2b) surfers, same start cells
        D, spread, bits = [], [], 0
        tick = 0
        for T in times:
            while tick * a < T - 1e-9:
                r0 = density(psi)
                psi, subs = W.tick(psi)
                r1 = density(psi)
                F = net_flux(0, subs[0][1], subs[0][2])
                flow_tick(cell, xi, r0, r1, [F])
                bellflux_tick(cb, r0, r1, [F], rng)
                tick += 1
            xc, R, Lw = cont_density(T)
            Qc = quantile_function_cont(xc, np.abs(R)**2 + np.abs(Lw)**2, u)
            Ql = quantile_function_latt(a, density(psi), u)
            D.append(float(np.max(np.abs(Qc - Ql))))
            xs = (cell[:, 0] + xi[:, 0] - n / 2) * a
            xbell = (cb[:, 0] + 0.5 - n / 2) * a
            spread.append(float(np.sqrt(np.mean((xbell - xs)**2))))
        out[str(a)] = dict(D_max=max(D), D_t=D, bell_rms_vs_flow_T40=spread[-1], bell_rms_t=spread)
        print(a, out[str(a)]["D_max"], spread[-1], flush=True)
    Ds = [out[k]["D_max"] for k in ("0.2", "0.1", "0.05")]
    order = float(np.log2(Ds[1] / Ds[2]))
    # check of the identity "continuum 1-D Bohm-Dirac keeps the quantile": integrate 20 trajectories directly
    xc, R0, L0 = cont_density(0.0)
    rho0 = np.abs(R0)**2 + np.abs(L0)**2
    u20 = np.linspace(0.05, 0.95, 19)
    X = quantile_function_cont(xc, rho0, u20)
    dt = 0.02
    def vfield(t):
        xc_, R, Lw = cont_density(t)
        r = np.abs(R)**2 + np.abs(Lw)**2
        return xc_, (np.abs(R)**2 - np.abs(Lw)**2) / np.maximum(r, 1e-300), r
    t = 0.0
    for k in range(int(round((40.0 if not SMOKE else 0.1) / dt))):
        xa, va, _ = vfield(t); xb, vb, _ = vfield(t + dt / 2); xd, vd, _ = vfield(t + dt)
        k1 = np.interp(X, xa, va); k2 = np.interp(X + dt / 2 * k1, xb, vb); k3 = np.interp(X + dt / 2 * k2, xb, vb)
        k4 = np.interp(X + dt * k3, xd, vd)
        X = X + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
    _, _, r40 = vfield(40.0 if not SMOKE else round(0.1 / dt) * dt)
    h = xc[1] - xc[0]
    C = np.concatenate([[0], np.cumsum(r40)]); C /= C[-1]
    edges = np.concatenate([xc - h / 2, [xc[-1] + h / 2]])
    q40 = np.interp(X, edges, C)
    res = dict(per_a=out, order=order, D=Ds, bohm_quantile_identity_err=float(np.max(np.abs(q40 - u20))))
    print(res["order"], Ds, res["bohm_quantile_identity_err"])
    save("cont", res)


# ======================================================================================================
#  radial helpers (E3, E4)
# ======================================================================================================
def radial_tick_all(rad, chi, ens, rng, stats):
    """one tick of the radial rule; move every ensemble. ens: dict name -> dict(kind, cell, xi)"""
    r0 = density(chi)
    chi2, subs = rad.tick(chi)
    r1 = density(chi2)
    _, A, B = subs[0]
    F = net_flux(0, A, B, periodic=False)
    for name, e in ens.items():
        k = e["kind"]
        if k == "flow":
            flow_tick(e["cell"], e["xi"], r0, r1, [F], periodic=False)
        elif k == "bellflux":
            nj = bellflux_tick(e["cell"], r0, r1, [F], rng, periodic=False)
            stats[name]["jumps"] += int(nj.sum())
            e["ever"] |= nj > 0
        elif k in ("bell", "ballistic"):
            c, nj, H = bell_substep(e["cell"], 0, A, B, r0, 0 if k == "bell" else 1, rng, periodic=False)
            e["ever"] |= np.any(c != e["cell"], axis=1)
            e["cell"] = c
            stats[name]["jumps"] += nj
            stats[name]["bits"] += H
        elif k == "old":
            # old reader along r: v = (|chi+|^2 - |chi-|^2)/rho from linear interpolation, midpoint, tick boundary
            r = rad.rw.r / rad.rw.a                          # site centres in cells (i + 0.5)
            def v(p):
                cp = np.interp(p, r, chi2[0].real) + 1j * np.interp(p, r, chi2[0].imag)
                cm = np.interp(p, r, chi2[1].real) + 1j * np.interp(p, r, chi2[1].imag)
                return (np.abs(cp)**2 - np.abs(cm)**2) / np.maximum(np.abs(cp)**2 + np.abs(cm)**2, 1e-300)
            p = e["cell"][:, 0] + e["xi"][:, 0]
            v1 = v(p); v2 = v(p + 0.5 * v1)
            p = np.clip(p + v2, 0, len(r) - 1e-9)
            e["cell"][:, 0] = np.floor(p).astype(np.int64); e["xi"][:, 0] = p - np.floor(p)
    return chi2


def make_ens(rng, rho, kinds, N, start="born"):
    ens = {}
    w = rho if start == "born" else rho**2
    for k in kinds:
        c, x = born_sample(rng, w, N)
        ens[k] = dict(kind=k, cell=c, xi=x, ever=np.zeros(N, bool))
    return ens


def run_h1s():
    """E3: hydrogen (Z = 1, alpha) 1s of the radial rule: stationary real state -> surfers must stay still"""
    rng = np.random.default_rng(SEED + 3)
    a = 0.5
    om, chi, rw = radial_state(-1, 1, a, ALPHA, 2000.0)
    rad = Radial(rw)
    rho = density(chi)
    kinds = ("flow", "bellflux", "bell", "ballistic", "old")
    ens = make_ens(rng, rho, kinds, 20000)
    q0 = {k: quantile_of(e["cell"], e["xi"], rho) for k, e in ens.items()}
    stats = {k: dict(jumps=0, bits=0.0) for k in kinds}
    T = 10000 if not SMOKE else 20
    t0 = time.time()
    for t in range(T):
        chi = radial_tick_all(rad, chi, ens, rng, stats)
    rho = density(chi)
    res = dict(omega=om, dirac=float(1 / np.sqrt(1 + (ALPHA / (1 + np.sqrt(1 - ALPHA**2)))**2)), a=a, ticks=T,
               surfers=20000, time_s=time.time() - t0)
    for k, e in ens.items():
        q = quantile_of(e["cell"], e["xi"], rho)
        res[k] = dict(max_dq=float(np.max(np.abs(q - q0[k]))), rms_dq=float(np.sqrt(np.mean((q - q0[k])**2))),
                      jumps_per_surfer_per_tick=stats[k]["jumps"] / (20000 * T),
                      bits_per_surfer_per_tick=stats[k]["bits"] / (20000 * T),
                      fraction_ever_jumped=float(e["ever"].mean()), chi2_z=chi2_z(e["cell"], rho, 8)[0])
        print(k, res[k], flush=True)
    save("h1s", res)


def run_beat():
    """E4 + E8 (1-D part): radial Z = 60, (1s + 2s)/sqrt 2: the density breathes with the beat."""
    rng = np.random.default_rng(SEED + 4)
    a, Za = 0.05, 60 * ALPHA
    o1, c1, rw = radial_state(-1, 1, a, Za, 40.0)
    o2, c2, _ = radial_state(-1, 2, a, Za, 40.0)
    c1 /= np.linalg.norm(c1); c2 /= np.linalg.norm(c2)
    chi = (c1 + c2) / np.sqrt(2)
    rad = Radial(rw)
    rho = density(chi)
    per = 2 * np.pi / abs(o2 - o1)
    T = int(round(20 * per / a)) if not SMOKE else 40
    kinds = ("flow", "bellflux", "bell", "ballistic", "old")
    Nf, Ns = 100000, 20000
    ens = make_ens(rng, rho, ("flow",), Nf)
    ens.update(make_ens(rng, rho, kinds[1:], Ns))
    neq = make_ens(rng, rho, ("flow", "bellflux", "bell", "ballistic"), Ns, start="rho2")
    allens = dict(ens); allens.update({"neq_" + k: e for k, e in neq.items()})
    stats = {k: dict(jumps=0, bits=0.0) for k in allens}
    q0 = quantile_of(ens["flow"]["cell"], ens["flow"]["xi"], rho)
    qn0 = {k: quantile_of(e["cell"], e["xi"], rho) for k, e in neq.items()}
    H0 = {k: quantile_bins_H(v) for k, v in qn0.items()}
    rec = {k: [] for k in kinds}
    recH = {k: [] for k in neq}
    rec_qerr = []
    checks = np.unique(np.linspace(0, T, 41).astype(int))
    t0 = time.time()
    for t in range(T + 1):
        if t in checks:
            rho = density(chi)
            for k in kinds:
                rec[k].append(chi2_z(ens[k]["cell"], rho, 4)[0])
            rec_qerr.append(float(np.max(np.abs(quantile_of(ens["flow"]["cell"], ens["flow"]["xi"], rho) - q0))))
            for k, e in neq.items():
                recH[k].append(quantile_bins_H(quantile_of(e["cell"], e["xi"], rho)))
        if t < T:
            chi = radial_tick_all(rad, chi, allens, rng, stats)
    res = dict(omega_1s=o1, omega_2s=o2, beat_period=per, ticks=T, periods=20, a=a, Za=Za, time_s=time.time() - t0,
               chi2_z={k: v for k, v in rec.items()}, max_abs_z={k: float(np.max(np.abs(v))) for k, v in rec.items()},
               flow_quantile_err_max=max(rec_qerr),
               neq_H_ratio={k: float(v[-1] / H0[k]) for k, v in recH.items()}, neq_H={k: v for k, v in recH.items()},
               old_final_z=rec["old"][-1],
               jumps_per_surfer_per_tick={k: stats[k]["jumps"] / (len(allens[k]["cell"]) * T) for k in allens},
               bits_per_surfer_per_tick={k: stats[k]["bits"] / (len(allens[k]["cell"]) * T) for k in allens})
    print(json.dumps({k: res[k] for k in ("max_abs_z", "flow_quantile_err_max", "neq_H_ratio", "jumps_per_surfer_per_tick")}, indent=1, default=float))
    save("beat", res)


# ======================================================================================================
#  3-D runs (E5, E6, E7, E8-3D), checkpointed
# ======================================================================================================
N3, A3 = 48, 0.5
ZA = 60 * ALPHA


def grid3():
    x = (np.arange(N3) - N3 / 2 + 0.5) * A3
    return np.meshgrid(x, x, x, indexing="ij")


def phys(cell, xi):
    return (cell + xi - N3 / 2) * A3


def radius(cell, xi):
    p = phys(cell, xi)
    return np.sqrt(np.sum(p**2, 1))


def shells_H(r, edges):
    nb = len(edges) - 1
    c = np.histogram(r, bins=edges)[0] / len(r)
    m = c > 0
    return float(np.sum(c[m] * np.log(c[m] * nb)))


def radial_edges(rho, nb=20):
    """equal-probability radial shells of the piecewise-uniform measure (cell centres, fine enough here)"""
    X, Y, Z = grid3()
    R = np.sqrt(X**2 + Y**2 + Z**2).ravel()
    o = np.argsort(R)
    C = np.cumsum(rho.ravel()[o]); C /= C[-1]
    e = np.interp(np.linspace(0, 1, nb + 1)[1:-1], C, R[o])
    return np.concatenate([[0], e, [1e9]])


def setup_3d(state):
    X, Y, Z = grid3()
    if state == "free":
        W = Walk3D(N3, A3, V=None)
        k0 = np.array([0.6, 0.3, 0.0])
        sig = 3.0
        g = np.exp(-(X**2 + Y**2 + Z**2) / (4 * sig**2) + 1j * (k0[0] * X + k0[1] * Y + k0[2] * Z))
        psi = np.array([g, 0 * g, g, 0 * g]) / np.sqrt(2)    # L = R = (1, 0): spin up, rest-frame positive energy
        psi /= np.linalg.norm(psi)
        return W, psi
    R = np.sqrt(X**2 + Y**2 + Z**2)
    W = Walk3D(N3, A3, V=-ZA / R * A3)
    s1 = np.load(os.path.join(CK, "state_1s.npy"))
    if state == "1s":
        return W, s1 / np.linalg.norm(s1)
    s2 = np.load(os.path.join(CK, "state_2p.npy"))
    psi = s1 / np.linalg.norm(s1) + s2 / np.linalg.norm(s2)
    return W, psi / np.linalg.norm(psi)


KINDS3 = ("flow2", "flow", "bellflux", "bell", "ballistic", "cdfasse", "old", "old_rebirth")


def tick3_all(W, psi, ens, rng, stats):
    r0 = density(psi)
    subs_rho = []                        # densities before each sub-step (for cdfasse)
    psi2, subs = W.tick(psi)
    r1 = density(psi2)
    F = [None] * 3
    for (j, A, B) in subs:
        F[j] = net_flux(j, A, B)
    # densities between sub-steps: rho_after_j = rho_before_j - div_j F_j
    rr = [r0]
    for (j, A, B) in subs:
        rr.append(rr[-1] - (F[j] - np.roll(F[j], 1, axis=j)))
    first_of_pair = (W.tick_count - 1) % 2 == 0          # this tick had axis order x, y, z
    for name, e in ens.items():
        k = e["kind"]
        if k == "flow":
            nc, st = flow_tick(e["cell"], e["xi"], r0, r1, F)
            stats[name]["cross"] += int(nc.sum()); stats[name]["stuck"] += int(st)
        elif k == "flow2":
            # S4 over the two-tick period of the rule (x,y,z then z,y,x): flux of the pair, rho(2m) -> rho(2m+2)
            if first_of_pair:
                e["F2"], e["r_start"] = [f.copy() for f in F], r0
            else:
                F2 = [e["F2"][j] + F[j] for j in range(3)]
                nc, st = flow_tick(e["cell"], e["xi"], e["r_start"], r1, F2)
                stats[name]["cross"] += int(nc.sum()); stats[name]["stuck"] += int(st)
        elif k == "bellflux":
            nj = bellflux_tick(e["cell"], r0, r1, F, rng)
            stats[name]["jumps"] += int(nj.sum())
        elif k in ("bell", "ballistic"):
            for s, (j, A, B) in enumerate(subs):
                c, nj, H = bell_substep(e["cell"], j, A, B, rr[s], 0 if k == "bell" else 1, rng)
                e["cell"] = c
                stats[name]["jumps"] += nj; stats[name]["bits"] += H
        elif k == "cdfasse":
            for s, (j, A, B) in enumerate(subs):
                cdf_substep(e["cell"], e["xi"], j, rr[s], rr[s + 1])
        elif k in ("old", "old_rebirth"):
            Xc = e["cell"] + e["xi"]
            Xn = old_reader_move(psi2, Xc)
            e["cell"] = np.floor(Xn).astype(np.int64) % N3
            e["xi"] = Xn - np.floor(Xn)
            if k == "old_rebirth":            # Rust: 40 of 20000 reborn from |psi|^2 every tick (0.2 %)
                m = max(1, int(round(0.002 * len(Xc))))
                idx = rng.choice(len(Xc), m, replace=False)
                c, x = born_sample(rng, r1, m)
                e["cell"][idx] = c; e["xi"][idx] = x
    return psi2


def run_a3d(state, budget=540.0):
    t_start = time.time()
    T = dict(free=300, **{"1s": 1000, "1s2p": 1670})[state] if not SMOKE else 12
    ck = os.path.join(CK, f"a3d_{state}.pkl")
    W, psi = setup_3d(state)
    if os.path.exists(ck):
        S = pickle.load(open(ck, "rb"))
        psi, W.tick_count = S["psi"], S["t"]
        rng = S["rng"]
    else:
        rng = np.random.default_rng(SEED + 5)
        rho = density(psi)
        ens = {}
        for k in KINDS3:
            N = 200000 if k in ("flow", "flow2") else 50000
            c, x = born_sample(rng, rho, N)
            ens[k] = dict(kind=k, cell=c, xi=x)
        if state == "1s":
            for k in ("flow2", "flow", "bellflux", "bell", "ballistic", "cdfasse"):
                c, x = born_sample(rng, rho**2, 50000)
                ens["neq_" + k] = dict(kind=k, cell=c, xi=x)
        S = dict(t=0, psi=psi, rng=rng, ens=ens, stats={k: dict(jumps=0, bits=0.0, cross=0, stuck=0) for k in ens},
                 rec=dict(t=[], z={k: [] for k in KINDS3}, tv={k: [] for k in KINDS3}), r0={}, rT={}, H={}, vphi=[],
                 tvnull=None)
        S["edges"] = radial_edges(rho)
        S["r_start"] = {k: radius(e["cell"], e["xi"]) for k, e in ens.items()}
        S["H"] = {k: [shells_H(S["r_start"][k], S["edges"])] for k in ens if k.startswith("neq_")}
        S["tvnull"] = tv_null(rng, rho, 200000, 2, 10)
        S["tvnull_small"] = tv_null(rng, rho, 50000, 2, 10)
        S["drift_r"] = {k: [] for k in ("flow2", "flow", "bellflux", "bell", "ballistic", "cdfasse", "old")}
    ens = S["ens"]
    checks = set((np.unique((np.linspace(0, T, 41) / 2).astype(int)) * 2).tolist())     # even ticks only
    while S["t"] <= T and time.time() - t_start < budget:
        t = S["t"]
        if t in checks:
            rho = density(psi)
            S["rec"]["t"].append(t)
            for k in KINDS3:
                z, _, _, tv = chi2_z(ens[k]["cell"], rho, 2)
                S["rec"]["z"][k].append(z)
                S["rec"]["tv"][k].append(tv)
            for k in S["H"]:
                S["H"][k].append(shells_H(radius(ens[k]["cell"], ens[k]["xi"]), S["edges"]))
            if state == "1s":
                for k in S["drift_r"]:
                    dr = radius(ens[k]["cell"], ens[k]["xi"]) - S["r_start"][k]
                    S["drift_r"][k].append(float(np.sqrt(np.mean(dr**2))))
        if state == "1s" and t in (500, 510):
            # azimuthal / radial velocity over 10 ticks (surfers with 1 < r < 6, |cos theta| < 0.5)
            if t == 500:
                S["p500"] = {k: phys(ens[k]["cell"], ens[k]["xi"]) for k in ("flow2", "flow", "old", "bellflux", "cdfasse")}
            else:
                out = {}
                for k, p0 in S["p500"].items():
                    p1 = phys(ens[k]["cell"], ens[k]["xi"])
                    r = np.sqrt(np.sum(p0**2, 1)); rc = np.sqrt(p0[:, 0]**2 + p0[:, 1]**2)
                    sel = (r > 1) & (r < 6) & (np.abs(p0[:, 2]) < 0.5 * r)
                    ph0 = np.arctan2(p0[:, 1], p0[:, 0]); ph1 = np.arctan2(p1[:, 1], p1[:, 0])
                    dphi = (ph1 - ph0 + np.pi) % (2 * np.pi) - np.pi
                    vphi = rc * dphi / (10 * A3)
                    r1 = np.sqrt(np.sum(p1**2, 1))
                    out[k] = dict(vphi_over_Za_sin=float(np.mean(vphi[sel] / (rc[sel] / r[sel])) / ZA),
                                  vr_rms_over_Za=float(np.sqrt(np.mean(((r1 - r) / (10 * A3))[sel]**2)) / ZA))
                S["vphi"] = out
                print("v_phi:", out, flush=True)
        if t == T:
            S["t"] += 1
            break
        psi = tick3_all(W, psi, ens, rng, S["stats"])
        S["t"] += 1
    S["psi"] = psi
    S["rng"] = rng
    if S["t"] <= T:
        pickle.dump(S, open(ck, "wb"))
        print(f"a3d {state}: checkpoint at tick {S['t']}/{T}", flush=True)
        sys.exit(3)
    # Liouville + reversibility checks for the flow (state-specific, a few ticks from the final state)
    res = dict(state=state, ticks=T, n=N3, a=A3, rec=S["rec"], tv_null_200k=S["tvnull"], tv_null_50k=S["tvnull_small"],
               max_abs_z={k: float(np.max(np.abs(S["rec"]["z"][k]))) for k in KINDS3},
               tv_ratio_final={k: S["rec"]["tv"][k][-1] / (S["tvnull"][0] if k in ("flow", "flow2") else S["tvnull_small"][0]) for k in KINDS3},
               stats={k: {q: v / (len(ens[k]["cell"]) * T) for q, v in s.items()} for k, s in S["stats"].items()})
    if state == "1s":
        res["H_neq"] = S["H"]
        res["H_ratio"] = {k: v[-1] / v[0] for k, v in S["H"].items()}
        res["drift_r_rms"] = S["drift_r"]
        res["vphi"] = S["vphi"]
        rt = S["rec"]["t"]
        i100 = int(np.argmin(np.abs(np.array(rt) - 100)))
        res["drift_growth_1000_over_100"] = {k: v[-1] / max(v[i100], 1e-300) for k, v in S["drift_r"].items()}
    res.update(liouville_reversibility(W, psi, rng))
    save(f"a3d_{state}", res)
    if os.path.exists(ck):
        os.remove(ck)


def liouville_reversibility(W, psi, rng, ticks=20, M=4000, h=1e-6):
    st = []
    for _ in range(ticks):
        r0 = density(psi)
        psi, subs = W.tick(psi)
        r1 = density(psi)
        F = [None] * 3
        for (j, A, B) in subs:
            F[j] = net_flux(j, A, B)
        st.append((F, r0, r1))
    # Liouville over ONE tick
    c0, x0 = born_sample(rng, st[0][1], M)
    x0 = np.clip(x0, 0.01, 0.99)
    pts = [(c0.copy(), x0.copy())]
    for a_ in range(3):
        xx = x0.copy(); xx[:, a_] += h
        pts.append((c0.copy(), xx))
    for (c, x) in pts:
        flow_tick(c, x, st[0][1], st[0][2], st[0][0])
    P = [c + x for c, x in pts]
    Jm = np.stack([((P[k + 1] - P[0] + N3 / 2) % N3 - N3 / 2) / h for k in range(3)], axis=2)
    J = np.abs(np.linalg.det(Jm))
    e = np.abs(st[0][2][tuple(pts[0][0].T)] * J / st[0][1][tuple(c0.T)] - 1)
    # reversibility over `ticks` ticks
    c, x = born_sample(rng, st[0][1], M)
    cs, xs = c.copy(), x.copy()
    for F, r0, r1 in st:
        flow_tick(c, x, r0, r1, F)
    for F, r0, r1 in st[::-1]:
        flow_tick(c, x, r1, r0, [-f for f in F])
    d = np.abs((c + x) - (cs + xs)); d = np.minimum(d, N3 - d).max(1)
    return dict(liouville_median=float(np.median(e)), liouville_q99=float(np.quantile(e, 0.99)),
                reverse_median=float(np.median(d)), reverse_q99=float(np.quantile(d, 0.99)), reverse_max=float(d.max()))


# ======================================================================================================
#  E9: analytic orbitals (Rust orbitali.rs): continuous-time psi(t) on a grid, flux = Dirac current at faces +
#      the minimal (Poisson) correction that makes it exactly consistent with rho(t) -> rho(t + dt)
# ======================================================================================================
def poisson_corrected_flux(sp_mid, P0, P1, dt_cells):
    """sp_mid (4, n, n, n) chiral spinor at mid-step; P0, P1 cell probabilities (sum 1). Face flux in probability
    per step: F_j = dt * (j_c + j_{c+e_j})/2 / S, then F -> F - grad phi with L phi = div F - (P0 - P1)."""
    n = P0.shape[0]
    L, R = sp_mid[0:2], sp_mid[2:4]
    S = np.sum(np.abs(sp_mid)**2)
    F = []
    for j, M in enumerate((SXm, SYm, SZm)):
        jr = np.real(np.einsum("a...,ab,b...->...", R.conj(), M, R))
        jl = np.real(np.einsum("a...,ab,b...->...", L.conj(), M, L))
        jc = (jr - jl) / S
        F.append(dt_cells * 0.5 * (jc + np.roll(jc, -1, axis=j)))
    div = sum(F[j] - np.roll(F[j], 1, axis=j) for j in range(3))
    rhs = div - (P0 - P1)
    k = 2 * np.pi * np.fft.fftfreq(n)
    lam = sum((2 * np.cos(kk) - 2) for kk in np.meshgrid(k, k, k, indexing="ij"))
    lam[0, 0, 0] = 1.0
    ph = np.fft.ifftn(np.fft.fftn(rhs) / lam).real            # L phi = rhs, L = div grad (negative definite)
    F = [F[j] - (np.roll(ph, -1, axis=j) - ph) for j in range(3)]
    return F, float(np.sqrt(np.mean(rhs**2)))


def run_orb():
    """analytic 1s + 2p1/2 (radial rule levels, continuous time, as orbitali.rs) on a 48^3 grid at a = 0.5"""
    rng = np.random.default_rng(SEED + 9)
    X, Y, Z = grid3()
    o1, c1, rw1 = radial_state(-1, 1, A3, ZA, 40.0)
    o2, c2, rw2 = radial_state(1, 2, A3, ZA, 40.0)
    s1 = spinor3d_from_radial(rw1, c1, -1, X, Y, Z); s1 /= np.linalg.norm(s1)
    s2 = spinor3d_from_radial(rw2, c2, 1, X, Y, Z); s2 /= np.linalg.norm(s2)
    per = 2 * np.pi / abs(o2 - o1)
    dt = per / 200.0
    def psi_at(t):
        return (s1 * np.exp(-1j * (o1 - 1) * t) + s2 * np.exp(-1j * (o2 - 1) * t)) / np.sqrt(2)
    P = lambda sp: density(sp) / np.sum(density(sp))
    steps = 2000 if not SMOKE else 4                 # 10 beat periods
    rho = P(psi_at(0.0))
    c, x = born_sample(rng, rho, 200000)
    co, xo = born_sample(rng, rho, 50000)
    zr, zo, zrb, res_rms = [], [], [], []
    cob, xob = co.copy(), xo.copy()
    for s in range(steps + 1):
        t = s * dt
        if s % 50 == 0:
            r = P(psi_at(t))
            zr.append(chi2_z(c, r, 2)[0]); zo.append(chi2_z(co, r, 2)[0]); zrb.append(chi2_z(cob, r, 2)[0])
        if s == steps:
            break
        P0, P1 = P(psi_at(t)), P(psi_at(t + dt))
        F, rr = poisson_corrected_flux(psi_at(t + dt / 2), P0, P1, dt / A3)
        res_rms.append(rr)
        flow_tick(c, x, P0, P1, F)
        # old orbital reader: 4 midpoint sub-steps with trilinear spinor (orbitali.rs), with and without rebirth
        for (cc, xx, reb) in ((co, xo, False), (cob, xob, True)):
            Xc = cc + xx
            sp = psi_at(t)
            hh = dt / A3 / 4
            for _ in range(4):
                v1 = bohm_dirac_v(trilinear(sp, Xc)); q = Xc + 0.5 * hh * v1
                v2 = bohm_dirac_v(trilinear(sp, q)); Xc = (Xc + hh * v2) % N3
            cc[:] = np.floor(Xc).astype(np.int64) % N3; xx[:] = Xc - np.floor(Xc)
            if reb:                                  # orbitali.rs: 25 of 5000 reborn per step (0.5 %)
                m = int(0.005 * len(cc)); idx = rng.choice(len(cc), m, replace=False)
                c2_, x2_ = born_sample(rng, P(psi_at(t + dt)), m); cc[idx] = c2_; xx[idx] = x2_
    res = dict(omega=[o1, o2], period=per, steps=steps, z_flow=zr, z_old=zo, z_old_rebirth=zrb,
               max_abs_z_flow=float(np.max(np.abs(zr))), max_abs_z_old=float(np.max(np.abs(zo))),
               max_abs_z_old_rebirth=float(np.max(np.abs(zrb))), poisson_correction_rms=float(np.mean(res_rms)))
    print({k: v for k, v in res.items() if k.startswith("max") or k.startswith("poisson")})
    save("orb", res)


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "f6":
        run_f6()
    elif what == "cont":
        run_cont()
    elif what == "h1s":
        run_h1s()
    elif what == "beat":
        run_beat()
    elif what == "a3d":
        run_a3d(sys.argv[2])
    elif what == "orb":
        run_orb()
