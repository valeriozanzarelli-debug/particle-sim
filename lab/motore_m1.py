"""motore_m1.py - exams A1, A2 (error model), B1 (level choice from the error on the screen), B2 (the memory of the
error: camera approaching, naive policy vs 'passport'). Pre-registration: lab/out/motore_prereg.json.
Usage: python lab/motore_m1.py A | B1 | B2      -> lab/out/motore_m1_<part>.json"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import motore_core as mc
import zz_core as zz

PRE = json.load(open(os.path.join(mc.OUT, "motore_prereg.json")))
CH = 1 << 16


# ---------------- chunked Fourier-exact evolutions (memory-light) ----------------
def kgrid_n(n):
    return mc.kgrid2(n)


def evolve_lattice(F, ticks):
    """exact `ticks` (even) ticks of the light walk, chunked over k"""
    n = F.shape[0]
    Fk = mc.fft2F(F)
    KA = kgrid_n(n)
    out = np.empty_like(Fk)
    for a in range(0, len(KA), CH):
        W = zz.bloch_batch(KA[a:a + CH], 0.0, spin=1)
        U = mc.matpow(W, ticks // 2)
        out[a:a + CH] = np.einsum("nab,nb->na", U, Fk[a:a + CH])
    return mc.ifft2F(out, n)


def evolve_continuum(F, t):
    n = F.shape[0]
    Fk = mc.fft2F(F)
    KA = kgrid_n(n)
    S = np.array(mc.SPIN1)
    out = np.empty_like(Fk)
    for a in range(0, len(KA), CH):
        K = KA[a:a + CH]
        k = np.linalg.norm(K, axis=1)
        kh = K / np.where(k > 0, k, 1.0)[:, None]
        M = np.einsum("nj,jab->nab", kh, S)
        th = (k * t)[:, None, None]
        U = np.eye(3)[None] - 1j * np.sin(th) * M + (np.cos(th) - 1) * (M @ M)
        out[a:a + CH] = np.einsum("nab,nb->na", U, Fk[a:a + CH])
    return mc.ifft2F(out, n)


def resample(F, n2):
    """spectral (conservative) transfer of a periodic field sampled at cell centres (i+1/2)h on an n1 grid to an n2
    grid of the SAME box: Fourier coefficients kept for |m| < min(n1, n2)/2 (Parseval: energy exact for band-limited
    content), Nyquist bin dropped."""
    n1 = F.shape[0]
    if n1 == n2:
        return F.copy()
    m1 = np.fft.fftfreq(n1) * n1
    ph1 = np.exp(-1j * np.pi * m1 / n1)                    # shift to cell centres: e^{-i k h/2}
    G = np.fft.fft2(F, axes=(0, 1)) * ph1[:, None, None] * ph1[None, :, None]
    mm = min(n1, n2)
    keep1 = np.abs(m1) < mm / 2
    m2 = np.fft.fftfreq(n2) * n2
    keep2 = np.abs(m2) < mm / 2
    H = np.zeros((n2, n2, 3), complex)
    i1, i2 = np.where(keep1)[0], np.where(keep2)[0]
    o1 = np.argsort(m1[i1]); o2 = np.argsort(m2[i2])
    i1, i2 = i1[o1], i2[o2]
    H[np.ix_(i2, i2)] = G[np.ix_(i1, i1)]
    ph2 = np.exp(1j * np.pi * m2 / n2)
    H = H * ph2[:, None, None] * ph2[None, :, None]
    return np.fft.ifft2(H, axes=(0, 1)) * (n2 * n2) / (n1 * n1)


def run_pair(N, ang, T_lambda, sigma_l, box_sig=4.5):
    """lattice (N cells per lambda) and continuum from the same data, co-moving frame; returns measures"""
    lam = N
    n = int(round(2 * box_sig * sigma_l * N))
    n += n % 2
    F = mc.light_packet(n, lam, ang, sigma_l * lam, sigma_l * lam)
    t = int(round(T_lambda * N))
    t += t % 2
    a = np.radians(ang)
    d = np.array([np.cos(a), np.sin(a)]) * t
    Fl = mc.shift_k(evolve_lattice(F, t), -d)
    Fc = mc.shift_k(evolve_continuum(F, t), -d)
    dx = (mc.centroid(Fl) - mc.centroid(Fc)) / N
    o = mc.overlap(Fc, Fl)
    pdx, pph = mc.predicted_drift_spectral(F, t)
    kdx, kph = mc.predicted_drift(N, ang, t / N)
    return dict(N=N, ang=ang, n=n, t=t, T_lambda=t / N, meas_dx=dx.tolist(), meas_phase=float(np.angle(o)),
                fidelity=float(abs(o) ** 2), spec_dx=(pdx / N).tolist(), spec_phase=pph, k0_dx=kdx.tolist(),
                k0_phase=kph, E_ratio=mc.energy(Fl) / mc.energy(F))


def wrap(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def judge_A(r, which):
    p = np.array(r[which + "_dx"]); m = np.array(r["meas_dx"])
    ok_dx = np.linalg.norm(p - m) <= 0.05 * np.linalg.norm(m) + 0.01
    ok_ph = abs(wrap(r[which + "_phase"] - r["meas_phase"])) <= 0.05
    return bool(ok_dx and ok_ph)


def part_A():
    res = []
    for N in (8, 12, 16, 32):
        for ang in (0.0, 22.5, 30.0, 45.0):
            t0 = time.time()
            r = run_pair(N, ang, 50.0, 3.0)
            r["A1_pass"] = judge_A(r, "spec")
            r["A2_pass"] = judge_A(r, "k0")
            r["wall"] = time.time() - t0
            res.append(r)
            print("N=%2d ang=%4.1f dx=%s spec=%s ph=%.4f/%.4f fid=%.4f A1=%s A2=%s %.0fs" % (
                N, ang, np.round(r["meas_dx"], 4), np.round(r["spec_dx"], 4), r["meas_phase"], wrap(r["spec_phase"]),
                r["fidelity"], r["A1_pass"], r["A2_pass"], r["wall"]), flush=True)
    n1 = sum(r["A1_pass"] for r in res); n2 = sum(r["A2_pass"] for r in res)
    out = dict(cases=res, A1=dict(passed=n1, of=len(res), verdict="PASS" if n1 == len(res) else "FAIL"),
               A2=dict(passed=n2, of=len(res), verdict="PASS" if n2 >= 14 else "FAIL"))
    json.dump(out, open(os.path.join(mc.OUT, "motore_m1_A.json"), "w"), indent=1)
    print("A1", out["A1"], "A2", out["A2"])


def screen_meas(r, p):
    dx = np.array(r["meas_dx"])
    ph = r["k0_phase"] + wrap(r["meas_phase"] - r["k0_phase"])      # unwrap near the prediction
    e, env, fr, res_ = mc.screen_error_px(dx, ph, p)
    return e, env, fr, res_


def part_B1():
    cfg = PRE["configs"]
    sel = PRE["predictions"]["B1_selection"]
    out = {}
    for tier in ("intel", "ultra"):
        rows = []
        for row in sel[tier]["rows"]:
            p, ang = row["footprint_lambda_per_px"], row["angle"]
            N = row["N"] if row["N"] is not None else 64
            t0 = time.time()
            r = run_pair(N, ang, cfg["horizon_lambda"], 2.0 if N >= 48 else 3.0)
            e, env, fr, resv = screen_meas(r, p)
            r.update(footprint=p, flagged=row["N"] is None, meas_px=e, env_px=env, fringe_px=fr, fringes_resolved=resv,
                     predicted_px=row["predicted_px"], wall=time.time() - t0)
            rows.append(r)
            print(tier, "p=%5.2f ang=%4.1f N=%d flagged=%s meas=%.3f px (env %.3f, fringe %.3f) pred=%s %.0fs" % (
                p, ang, N, r["flagged"], e, env, fr, row["predicted_px"], r["wall"]), flush=True)
        ok = [r for r in rows if not r["flagged"]]
        fl = [r for r in rows if r["flagged"]]
        maxpx = max(r["meas_px"] for r in ok)
        red = len(ok) * 64 ** 3 / sum(r["N"] ** 3 for r in ok)
        flagged_confirmed = all(r["meas_px"] > 0.5 for r in fl)
        verdict = (maxpx <= 1.0) and flagged_confirmed and (red >= 4 if tier == "intel" else True)
        out[tier] = dict(rows=rows, max_px_unflagged=maxpx, cost_reduction_unflagged=red,
                         flagged_confirmed=flagged_confirmed, verdict="PASS" if verdict else "FAIL")
        print(tier, "max px", maxpx, "reduction", red, "flagged confirmed", flagged_confirmed, out[tier]["verdict"])
    json.dump(out, open(os.path.join(mc.OUT, "motore_m1_B1.json"), "w"), indent=1)


def part_B2():
    """camera approaching: footprint p(t) = 5 -> 0.1 lambda/px linearly over 100 lambda; packet at 30 deg,
    sigma = 2 lambda; checks every 10 lambda. Box fixed in wavelengths (18 lambda); transfers spectral."""
    ang, T, p0, p1, every, sig = 30.0, 100.0, 5.0, 0.1, 10.0, 2.0
    ladder = [12, 16, 24, 32, 48, 64, 96, 128]
    boxl = 18.0
    pt = lambda tl: p0 + (p1 - p0) * tl / T
    a = np.radians(ang)
    vdir = np.array([np.cos(a), np.sin(a)])

    def choose_naive(p_now):
        for N in ladder:
            dx, dph = mc.predicted_drift(N, ang, every)
            if mc.screen_error_px(dx, dph, p_now)[0] <= 0.5:
                return N
        return ladder[-1]

    def choose_passport(acc_dx, acc_ph, tl):
        # worst case: the camera reaches p_min = p1 (bounded camera speed); error budget over the remaining time
        for N in ladder:
            dx, dph = mc.predicted_drift(N, ang, T - tl)
            if mc.screen_error_px(acc_dx + dx, acc_ph + dph, p1)[0] <= 0.5:
                return N
        return None

    results = {}
    for policy in ("naive", "passport"):
        t0 = time.time()
        tl = 0.0
        N = None
        F = None
        acc_dx, acc_ph = np.zeros(2), 0.0
        hist = []
        cost = 0.0
        F0 = None
        while tl < T - 1e-9:
            Nn = choose_naive(pt(tl)) if policy == "naive" else choose_passport(acc_dx, acc_ph, tl)
            if Nn is None:
                Nn = ladder[-1]
            n = int(round(boxl * Nn))
            if F is None:
                F = mc.light_packet(n, Nn, ang, sig * Nn, sig * Nn)
                F0 = (Nn, F.copy())
            elif Nn != N:
                F = resample(F, n)
            N = Nn
            ticks = int(round(every * N))
            F = mc.shift_k(evolve_lattice(F, ticks), -vdir * ticks)
            dx, dph = mc.predicted_drift(N, ang, every)
            acc_dx, acc_ph = acc_dx + dx, acc_ph + dph
            cost += (boxl * N) ** 2 * ticks
            hist.append(dict(t_lambda=tl, p=pt(tl), N=N))
            tl += every
        # reference: continuum from the initial data, sampled on the final grid
        Ni, Fi = F0
        Fi_f = resample(Fi, F.shape[0])
        tt = T * N
        Fc = mc.shift_k(evolve_continuum(Fi_f, tt), -vdir * tt)
        dxm = (mc.centroid(F) - mc.centroid(Fc)) / N
        o = mc.overlap(Fc, F)
        ph = acc_ph + wrap(float(np.angle(o)) - acc_ph)
        e, env, fr, resv = mc.screen_error_px(dxm, ph, p1)
        results[policy] = dict(history=hist, final_N=N, meas_dx_lambda=dxm.tolist(), meas_phase=ph,
                               predicted_dx=acc_dx.tolist(), predicted_phase=acc_ph, fidelity=float(abs(o) ** 2),
                               final_px=e, env_px=env, fringe_px=fr, cost_cell_ticks=cost, wall=time.time() - t0)
        print(policy, [h["N"] for h in hist], "final px %.2f (env %.2f fringe %.2f) pred phase %.3f meas %.3f cost %.3g %.0fs"
              % (e, env, fr, acc_ph, ph, cost, time.time() - t0), flush=True)
    results["verdict"] = dict(naive_fails=results["naive"]["final_px"] > 1.0, passport_pass=results["passport"]["final_px"] <= 1.0)
    results["verdict"]["B2"] = "PASS" if (results["verdict"]["naive_fails"] and results["verdict"]["passport_pass"]) else "FAIL"
    results["note"] = "passport ladder extended to N = 96, 128 (no level of {8..64} keeps 0.5 px at 0.1 lambda/px for 100 lambda: see B1 flags)"
    json.dump(results, open(os.path.join(mc.OUT, "motore_m1_B2.json"), "w"), indent=1)
    print(results["verdict"])


if __name__ == "__main__":
    part = sys.argv[1]
    {"A": part_A, "B1": part_B1, "B2": part_B2}[part]()
