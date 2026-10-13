"""motore_m4.py - 'dinamica conservata': the same scenes at the Ultra, Intel and below-the-floor levels against the
reference (continuum Maxwell for light, exact fine zig-zag walk for the electron). Pre-registration:
lab/out/motore_prereg.json (exams EL, EE).
Usage: python lab/motore_m4.py EL3 | EL4 | EE | judge   -> lab/out/motore_m4_<part>.json"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import motore_core as mc
from motore_m1 import evolve_lattice, evolve_continuum, wrap

PRE = json.load(open(os.path.join(mc.OUT, "motore_prereg.json")))
TIERS = PRE["tiers"]
OUT = lambda p: os.path.join(mc.OUT, "motore_m4_%s.json" % p)


def part_EL3():
    """diffraction: beam sigma_perp = 2 lambda (amplitude), sigma_par = 6 lambda, after 60 lambda"""
    res = []
    for tier in ("ultra", "intel", "sotto_pavimento"):
        N = TIERS[tier]["light_N"]
        for ang in (0.0, 45.0):
            t0 = time.time()
            n = int(round(54 * N)); n += n % 2
            F = mc.light_packet(n, N, ang, 6 * N, 2 * N)
            t = 60 * N
            a = np.radians(ang)
            d = np.array([np.cos(a), np.sin(a)]) * t
            Fl = mc.shift_k(evolve_lattice(F, t), -d)
            Fc = mc.shift_k(evolve_continuum(F, t), -d)
            w0 = mc.width_perp(F, ang)
            wl, wc = mc.width_perp(Fl, ang), mc.width_perp(Fc, ang)
            ratio = float(np.sqrt((wl ** 2 - w0 ** 2) / (wc ** 2 - w0 ** 2)))
            pred = PRE["predictions"]["bloch"][str(N)][str(ang)]["diffraction_factor"]
            row = dict(tier=tier, N=N, ang=ang, w0_cells=w0, w_lat=wl, w_cont=wc, spread_ratio=ratio, predicted=pred,
                       pass_=abs(ratio - 1) <= TIERS[tier]["tol"], wall=time.time() - t0)
            res.append(row)
            print("EL3", tier, N, ang, "ratio %.4f pred %.4f" % (ratio, pred), row["pass_"], "%.0fs" % row["wall"], flush=True)
    json.dump(res, open(OUT("EL3"), "w"), indent=1)


def part_EL4():
    import reticolo_r2 as r2
    res = []
    for tier in ("ultra", "intel", "sotto_pavimento"):
        N = TIERS[tier]["light_N"]
        for ang in (0, 30):
            for dr in ("f2c", "c2f"):
                t0 = time.time()
                r = r2.run(N, ang, dr, 2)
                tol = TIERS[tier]["tol"]
                v = r["v_after"]
                ok = (r["T"] >= 1 - tol / 20) and (abs(r["lost"]) <= tol / 20) and (r["hel_minus"] + r["hel_static"] <= tol / 20) \
                    and (v is not None and abs(v - 1) <= tol)
                r.update(tier=tier, pass_=bool(ok))
                res.append(r)
                print("EL4", tier, N, ang, dr, "T=%.5f R=%.1e lost=%.1e wrong=%.1e v=%s" % (
                    r["T"], r["R"], r["lost"], r["hel_minus"] + r["hel_static"], v), ok, "%.0fs" % (time.time() - t0), flush=True)
    json.dump(res, open(OUT("EL4"), "w"), indent=1)


def part_EE():
    th, lv = mc.load_levels_1d()
    nf = 196608
    xf = (np.arange(nf) - nf / 2 + 0.5) * th
    A0 = mc.A0
    res = dict(tiers={})
    # broad positive-energy state for the mode energies (sigma = a0/8, k0 = 1/a0) and the atomic packet (sigma = a0/2)
    psiA = mc.positive_spinor_packet(xf, 1 / A0, A0 / 2, th, x0=-10 * A0)
    psiK = mc.positive_spinor_packet(xf, 1 / A0, A0 / 8, th)
    # Zitterbewegung probe: pure R at rest, sigma = a0/2 (both bands)
    g = np.exp(-xf ** 2 / (4 * (A0 / 2) ** 2))
    psiZ = np.stack([g, 0 * g]).astype(complex)
    ts = 2 * (A0 / 2) ** 2
    for tier in ("ultra", "intel", "sotto_pavimento"):
        B = TIERS[tier]["electron_B"]
        tol = TIERS[tier]["tol"]
        rule = lv[B]
        nc = nf // B
        xc = (np.arange(nc) - nc / 2 + 0.5) * th * B
        t0 = time.time()
        m = int(round(ts / (th * B)))
        t = m * B * th
        # EE1, EE2, EE4
        cA0 = mc.block1(psiA, B)
        fA = mc.block1(mc.fine_evolve(psiA, th, m * B), B)
        lA = rule.step(cA0[None], m)[0]
        a0_, af, al = mc.x_mean_width(cA0, xc), mc.x_mean_width(fA, xc), mc.x_mean_width(lA, xc)
        v_f, v_l = (af[0] - a0_[0]) / t, (al[0] - a0_[0]) / t
        sp_f, sp_l = af[1] ** 2 - a0_[1] ** 2, al[1] ** 2 - a0_[1] ** 2
        ee1 = v_l / v_f - 1
        ee2 = np.sqrt(sp_l / sp_f) - 1
        ee4 = al[2] / a0_[2] - 1
        # EE3 mode energies
        mK = max(1, int(round(2000 / (th * B))))
        tK = mK * B * th
        cK0 = mc.block1(psiK, B)
        fK = mc.block1(mc.fine_evolve(psiK, th, mK * B), B)
        lK = rule.step(cK0[None], mK)[0]
        F0, Ff, Fl = (np.fft.fft(z.sum(0) if False else (z[0] + z[1]) / np.sqrt(2)) for z in (cK0, fK, lK))
        Kc = 2 * np.pi * np.fft.fftfreq(nc)
        ee3 = {}
        for ka in (0.5, 1.0, 2.0):
            Kt = ka / A0 * th * B
            i = int(np.argmin(np.abs(Kc - Kt)))
            i0 = 0
            ph = lambda Fz, j: np.angle(Fz[j] / F0[j])
            Ef = -wrap(ph(Ff, i) - ph(Ff, i0)) / tK
            El = -wrap(ph(Fl, i) - ph(Fl, i0)) / tK
            kk = Kc[i] / (th * B)
            Ed = np.sqrt(1 + kk ** 2) - 1
            ee3[str(ka)] = dict(k=kk, E_fine=Ef, E_level=El, E_dirac=Ed, rel_err_vs_fine=El / Ef - 1, rel_err_vs_dirac=El / Ed - 1)
        # EE5 ZB: P_R(t) at the coarse times over 40 time units
        nz = max(8, int(round(40 / (th * B))))
        cZ = mc.block1(psiZ, B)
        PRl, PRf = [], []
        st_l = cZ.copy()
        st_f = psiZ.copy()
        for s in range(nz):
            PRl.append(float(np.sum(np.abs(st_l[0]) ** 2) / np.sum(np.abs(st_l) ** 2)))
            fb = mc.block1(st_f, B)
            PRf.append(float(np.sum(np.abs(fb[0]) ** 2) / np.sum(np.abs(fb) ** 2)))
            st_l = rule.step(st_l[None], 1)[0]
            st_f = mc.fine_evolve(st_f, th, B)
        tau = th * B

        def domf(series):
            s = np.array(series) - np.mean(series)
            if np.allclose(s, 0):
                return 0.0, 0.0
            sp = np.abs(np.fft.rfft(s * np.hanning(len(s))))
            fr = np.fft.rfftfreq(len(s), d=tau) * 2 * np.pi
            j = int(np.argmax(sp[1:]) + 1)
            return float(fr[j]), float(np.std(s))
        zl, zf = domf(PRl), domf(PRf)
        ok = abs(ee1) <= tol and abs(ee2) <= tol and all(abs(v["rel_err_vs_fine"]) <= tol for v in ee3.values()) and abs(ee4) <= 1e-10
        res["tiers"][tier] = dict(B=B, cell_lambdaC=th * B, cell_a0=th * B / A0, tau=tau, t_spread=t, EE1_v_rel=ee1,
                                  v_fine=v_f, v_level=v_l, EE2_spread_rel=ee2, EE3=ee3, EE4_norm=ee4,
                                  EE5=dict(nyquist_omega=np.pi / tau, true_omega=2.0, level_peak_omega=zl[0], level_amp=zl[1],
                                           fine_sampled_peak_omega=zf[0], fine_sampled_amp=zf[1]),
                                  verdict="PASS" if ok else "FAIL", wall=time.time() - t0)
        print("EE", tier, json.dumps(res["tiers"][tier], default=float)[:900], flush=True)
    json.dump(res, open(OUT("EE"), "w"), indent=1, default=float)


def judge():
    A = json.load(open(os.path.join(mc.OUT, "motore_m1_A.json")))["cases"]
    E3 = json.load(open(OUT("EL3")))
    E4 = json.load(open(OUT("EL4")))
    EE = json.load(open(OUT("EE")))
    out = {}
    for tier in ("ultra", "intel", "sotto_pavimento"):
        N = TIERS[tier]["light_N"]
        tol = TIERS[tier]["tol"]
        cs = [c for c in A if c["N"] == N]
        sp, fr, en = {}, {}, []
        for c in cs:
            a = np.radians(c["ang"])
            v = np.array([np.cos(a), np.sin(a)]) * c["T_lambda"] + np.array(c["meas_dx"])
            sp[str(c["ang"])] = float(np.linalg.norm(v) / c["T_lambda"] - 1)
            ph = c["k0_phase"] + wrap(c["meas_phase"] - c["k0_phase"])
            fr[str(c["ang"])] = float(-ph / (2 * np.pi * c["T_lambda"]))
            en.append(abs(c["E_ratio"] - 1))
        el1 = all(abs(x) <= tol for x in sp.values())
        el2 = all(abs(x) <= tol for x in fr.values())
        el3 = all(r["pass_"] for r in E3 if r["tier"] == tier)
        el4 = all(r["pass_"] for r in E4 if r["tier"] == tier)
        el5 = max(en) <= 1e-10
        aniso = max(sp.values()) - min(sp.values())
        el6 = aniso <= tol and max(abs(x) for x in sp.values()) < tol / 2
        out[tier] = dict(N=N, tol=tol, EL1_speed=sp, EL1=el1, EL2_freq=fr, EL2=el2, EL3=el3,
                         EL3_rows=[(r["ang"], r["spread_ratio"]) for r in E3 if r["tier"] == tier], EL4=el4,
                         EL4_rows=[(r["ang"], r["direction"], r["T"], r["lost"], r["v_after"]) for r in E4 if r["tier"] == tier],
                         EL5=el5, EL5_max=max(en), EL6=el6, EL6_anisotropy=aniso,
                         EL_verdict="PASS" if all([el1, el2, el3, el4, el5, el6]) else "FAIL",
                         EE_verdict=EE["tiers"][tier]["verdict"])
        print(tier, json.dumps(out[tier], default=float))
    json.dump(out, open(OUT("judge"), "w"), indent=1, default=float)


if __name__ == "__main__":
    {"EL3": part_EL3, "EL4": part_EL4, "EE": part_EE, "judge": judge}[sys.argv[1]]()
