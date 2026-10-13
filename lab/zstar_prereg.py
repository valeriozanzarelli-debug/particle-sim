"""zstar_prereg.py - pre-registration of the study ZSTAR (docs/ZSTAR.md): why the radial rule with a point
nucleus makes the 1s dive at Z* = 131.37 instead of 137. Written and committed BEFORE the production walk runs
(lab/zstar_run_*.py). Writes lab/out/zstar_prereg.json.

THE PREDICTION (formula, lab/zstar_core.py):
  near the nucleus the rule is scale-free (angles kappa/(i+s), zeta/(i+s) per tick, no a). Its zero-energy massless
  stationary state is a product of SU(1,1) transfer matrices
        T_i = [[e^{i phi_i}/cos th_i, -tan th_i], [-tan th_i, e^{-i phi_i}/cos th_i]],  th_i = kappa/(i+s),
        phi_i = zeta/(i+s),   v_0 = (refl0, 1)
  and far out (a << r << 1) it is A (r/a)^gamma + B (r/a)^-gamma: the rule = continuum Dirac with the boundary
  condition tau = beta a^{2 gamma}, beta = B/A (a self-adjoint extension chosen by the lattice).
  Z* is where the lattice's zero mode has NO r^{+gamma} part:  F(zeta*) = lim_n n^-gamma l_+ . T_{n-1}..T_0 v_0 = 0
  (computed as a Wronskian with the subdominant solution, D_wronskian). Then zeta_cr(a) - zeta* ~ a^{2 gamma*}.
  The full curve Z_cr(a) is predicted by matching: lattice (A0, B0) at the origin + continuum Dirac (m = 1, E = -1)
  outside (predict_Zcr). Analytic corollaries: (i) Z*(kappa, refl0) = Z*(-kappa, -refl0) exactly (T(-kappa) =
  sz T(kappa) sz); (ii) the order of the two half-steps (W = M S vs S M) cannot change any eigenvalue (similar
  matrices); (iii) the Trotter split of the coin is identical at m = 0 (only O(a) changes); (iv) refl1 (outer
  wall at r = 40) is irrelevant.

SEEN before this file (declared):
  * beta(Z = 20, 60, 92) = 0.06628, 0.10301, 0.22037 (= NUCLEO_1S.md on origin/master: 0.06629, 0.10301, 0.22037);
  * the formula zeta* for the default rule (kappa = -1, s = 1/2, refl0 = -1): 0.9585907 (Z* = 131.3614) and the
    giro-1 walk values Z_cr(a) (lab/out/supercrit_c1.json), which the matching formula reproduces to
    0.041, 0.012, 0.004, 0.001 at a = 1e-3, 3e-4, 1e-4, 3e-5 (retrospective, NOT a blind test);
  * the formula roots for all variants below (they are the predictions).
  * the matching radius: r0 = 1e-8 (pure leading order) vs r0 = 10 a changes Z_match by 0.2 / 0.02 at
    a = 1e-3 / 1e-4 (base): this is the size of the neglected O(a) terms; r0 = 1e-8 is kept.
NOT seen: every walk run of a variant (refl0 != -1, s != 1/2, kappa != -1, coin split, order, refl1).
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.optimize import brentq
from zstar_core import D_wronskian, predict_Zcr, beta_fit
from zz_core import ALPHA

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "zstar_prereg.json")
S3 = float(np.sqrt(3) / 2)

VARIANTS = {   # name: (kappa, s, refl0, description)
    "base":   (-1, 0.5, -1.0, "frozen rule: 1s1/2, sites (i+1/2)a, reflection chi- -> -chi+ at r = 0"),
    "R+1":    (-1, 0.5, 1.0, "reflection sign +1"),
    "R+i":    (-1, 0.5, 1j, "reflection phase +i"),
    "R-i":    (-1, 0.5, -1j, "reflection phase -i"),
    "S0.25":  (-1, 0.25, -1.0, "sites (i+1/4)a: wall at -a/4 (nucleus inside the first cell)"),
    "S0.75":  (-1, 0.75, -1.0, "sites (i+3/4)a"),
    "S0.866": (-1, S3, -1.0, "first site at a sqrt(3)/2 = distance of the 3-D rule's first shell (nucleus at cube centre)"),
    "S1":     (-1, 1.0, -1.0, "sites (i+1)a: nucleus on the mirror cell (3-D analogue: nucleus on a lattice site)"),
    "K+1":    (1, 0.5, -1.0, "p1/2 (kappa = +1)"),
    "K-2":    (-2, 0.5, -1.0, "p3/2 (kappa = -2)"),
    "K+2":    (2, 0.5, -1.0, "d3/2 (kappa = +2)"),
    "K+1R+1": (1, 0.5, 1.0, "p1/2 with refl0 = +1 (duality check: must equal base)"),
}
A_RUN = [3e-3, 1e-3, 3e-4, 1e-4]


def zeta_roots(kappa, s, refl0, n=160, L=100000):
    k = abs(kappa)
    zz = np.linspace(0.01, k - 1e-3, n)
    d = np.array([D_wronskian(z, kappa=kappa, s=s, refl0=refl0, L=L) for z in zz])
    j = np.argmax(np.abs(d))
    ph = d[j] / abs(d[j])
    dr = (d * np.conj(ph))
    imag = float(np.max(np.abs(dr.imag)) / np.max(np.abs(dr)))
    dr = dr.real
    out = []
    for i in range(n - 1):
        if np.sign(dr[i]) != np.sign(dr[i + 1]):
            f = lambda z: (D_wronskian(z, kappa=kappa, s=s, refl0=refl0, L=L) * np.conj(ph)).real
            out.append(float(brentq(f, zz[i], zz[i + 1], xtol=1e-12)))
    return out, imag


if __name__ == "__main__":
    t0 = time.time()
    pred = {}
    for name, (kappa, s, refl0, desc) in VARIANTS.items():
        roots, imag = zeta_roots(kappa, s, refl0)
        # the root relevant for the dive: the one with 2 gamma < 1 if any (approached as a^{2 gamma}), else the lowest
        entry = dict(kappa=kappa, s=s, refl0=[float(np.real(refl0)), float(np.imag(refl0))], desc=desc,
                     zeta_roots=roots, Z_roots=[z / ALPHA for z in roots], D_imag_rel=imag,
                     two_gamma_at_roots=[float(2 * np.sqrt(kappa**2 - z**2)) for z in roots])
        if kappa == -1 or kappa == 1:
            zp = {}
            for a in A_RUN + [1e-6, 1e-9]:
                # zlo: below zeta ~ 0.8 the matching at r0 = 1e-8 is ill-conditioned (regular part lost under the
                # r^-gamma part when 2 gamma ~ 2): spurious sign flips seen for S1 in development -> start higher
                zlo = min(0.85, roots[0] - 0.05) if roots else 0.85
                z, e = predict_Zcr(a, kappa, s, refl0, zlo=zlo, zhi=1.25, n=60)
                zp[f"{a:.6g}"] = None if z is None else z / ALPHA
            entry["Zcr_matching"] = zp
        pred[name] = entry
        print(name, [f"{z:.4f}" for z in entry["Z_roots"]], entry.get("Zcr_matching"), f"{time.time() - t0:.0f}s",
              flush=True)
    # s-scan and reflection-phase scan of zeta* (kappa = -1), for the figure and for s_c
    s_scan = {}
    for s in np.round(np.arange(0.2, 1.21, 0.05), 3):
        r, _ = zeta_roots(-1, float(s), -1.0, n=80, L=50000)
        s_scan[f"{s:.3f}"] = r
    eta_scan = {}
    for eta in np.linspace(0, 2, 25)[:-1]:
        r, _ = zeta_roots(-1, 0.5, np.exp(1j * np.pi * eta), n=80, L=50000)
        eta_scan[f"{eta:.4f}"] = r
    # critical offset s_c(kappa = -1, refl0 = -1): zeta* -> 1
    def zs_of_s(s):
        r, _ = zeta_roots(-1, s, -1.0, n=60, L=50000)
        return r
    lo, hi = 0.5, 0.75
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        if zs_of_s(mid):
            lo = mid
        else:
            hi = mid
    s_c = 0.5 * (lo + hi)
    # graphene-analogue: 2-D massless Dirac radial equation has the same form with kappa -> j = +-1/2
    g2d = {}
    for kap in (-0.5, 0.5):
        r, _ = zeta_roots(kap, 0.5, -1.0, n=100)
        g2d[str(kap)] = dict(zeta_roots=r, ratio_to_continuum=[z / abs(kap) for z in r])

    P = {
        "study": "zstar", "doc": "docs/ZSTAR.md",
        "rule": "lab/zz_core.py RadialWalk (frozen) with point nucleus; variants only at the origin (refl0, offset s of "
                "the sites r_i = (i + s) a), kappa, order of steps, Trotter split of the coin, refl1",
        "inputs": "m_e, alpha, Z (point nucleus). No nuclear radius (the remedy S1/S0.866 uses a geometric choice "
                  "of the first site = DEBT if adopted)",
        "free_parameters": 0,
        "formula": "zeta* = root of the Wronskian D(zeta) = det[s_reg, s_sub] of the scale-free (a-free, m = 0, "
                   "lambda = 1) lattice recursion; Z_cr(a) from matching lattice (A0, B0) with continuum Dirac at "
                   "E = -1 (lab/zstar_core.py)",
        "predictions": pred,
        "s_scan_kappa-1": s_scan, "eta_scan_kappa-1_s0.5": eta_scan, "s_c_kappa-1": s_c,
        "graphene_analogue_j": g2d,
        "configs": {"a": A_RUN, "rmax": 40.0, "tracking": "supercrit_run_c1.track_to_critical (DZMAX 0.5), "
                    "start Z0 = 120 (or Z*_pred - 20 when lower; K-2 from 180)",
                    "K+2": "contact state: for a in A_RUN[1:] scan Z in [40, 44], gap eigenstates near omega = 0, "
                           "the one with > 50 % weight in r < 50 a; Z_x(a) where its omega = 0"},
        "exams": {
            "F1": "formula vs walk, retrospective (seen): base |Z_walk(a) - Z_match(a)| <= 0.01 at a = 3e-5 (C1 data); "
                  "recorded, not blind",
            "F2": "blind: every variant with a root (R+1, R+i, R-i, S0.25, K+1, K+1R+1): |Z_walk(1e-4) - Z_match(1e-4)| "
                  "<= 0.05 and |Z_walk(1e-3) - Z_match(1e-3)| <= 0.2",
            "F3": "blind: walk-only extrapolation Z_walk(a) = Z_inf + c1 a^(2 gamma*) + c2 a on the 3 finest a: "
                  "|Z_inf - Z*_formula| <= 0.1 for R+1, S0.25, K+1",
            "F4": "duality Z*(kappa, refl0) = Z*(-kappa, -refl0) (scale-free part only; mass and the continuum "
                  "differ at finite a): walk-only Z_inf (as F3) of K+1R+1 equal to base within 0.1, and of K+1 equal "
                  "to R+1 within 0.1; at finite a each one must match its OWN Z_match (F2 thresholds)",
            "U1": "universality test: Z* is universal for the rule iff max-min over {base, R+1, R+i, R-i, S0.25} < 1. "
                  "PREDICTION: NOT universal (spread ~36): Z* is a boundary artifact",
            "R1": "remedy: S1 and S0.866 have no root for kappa = +-1: Z_walk(a) >= 137.036 at every a in A_RUN and "
                  "decreasing as a decreases; |Z_walk(1e-4) - Z_match(1e-4)| <= 0.2",
            "K2": "p3/2 (K-2): the lowest kappa = -2 level dives at Z_walk(a) -> Z*_formula = 206.02: "
                  "|Z_walk(1e-4) - 206.02| <= 0.3 (2 gamma* = 2.6 > 1: approach dominated by O(a))",
            "K3": "d3/2 (K+2): a contact state crosses omega = 0 at Z_x(a) -> 42.05: |Z_x(1e-4) - 42.05| <= 0.1",
            "O1": "order S M vs M S: identical quasi-energies (|delta omega| < 1e-10) at base, a = 1e-3, Z = 132",
            "O2": "coin split (e^-iVa e^-i kappa/r sx a e^-i m sz a): |Z_walk(1e-4) - Z_walk(base, 1e-4)| <= 0.05 "
                  "(differs only at O(a))",
            "O3": "refl1 = +1 instead of -1: |delta Z_cr(1e-3)| <= 1e-3",
            "G1": "graphene analogue (j = +-1/2): reported, not an exam: zeta*/j for the radial rule",
        },
        "time": time.strftime("%Y-%m-%d %H:%M"),
    }
    json.dump(P, open(OUT, "w"), indent=1)
    print("s_c =", s_c, " graphene analogue:", g2d, f"{time.time() - t0:.0f}s")
