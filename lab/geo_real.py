"""geo_real.py - Part B2: real-data pilot (pre-registered in lab/out/geo_prereg.json, partB_real).

Rule: a measured mass is used for TRAINING only if the shape family has fewer parameters than the
numbers it must reproduce; the left-over measured numbers are the EXAM. A fit is not a prediction.
Data (PDG): charged leptons; quarks (scheme mix: light quarks MS-bar at 2 GeV, c and b MS-bar at their
own mass, top from direct measurement - a declared debt for any quark mass relation).
Tests:
  KK2        hidden 2-torus: charge circle (all charged leptons n1 = 1) x 'generation' circle n2.
             m(n2) = sqrt(M^2 + n2^2 mu^2). 2 parameters, trained on e (n2 = 0), mu (n2 = 1).
             exam: tau (n2 = 2) and the number of charged leptons below 2 GeV.
  KK2-free   same, but the level of the tau is left free (n2 = j, any integer); then with a second
             hidden circle (j = p^2 + q^2): counts how many discrete choices 'fit' - the look-elsewhere.
  KOIDE      three points on a circle: sqrt m_k = M (1 + sqrt2 cos(delta + 2 pi k / 3)).
             2 parameters per sector; train on the two lighter masses, exam the heaviest.
Writes lab/out/geo_real.json and .log.
"""
import os, json
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
ME, MMU = 0.51099895000, 105.6583755
MTAU, SMTAU = 1776.86, 0.12          # pre-registered value (PDG 2022)
MTAU24, SMTAU24 = 1776.93, 0.09      # PDG 2024 average, as a cross-check
QUARKS = {"up": (2.16, 1273.0, 172570.0), "down": (4.70, 93.5, 4183.0)}


def koide_third(m1, m2):
    """the circle family with amplitude sqrt2 <=> Q = (sum m)/(sum sqrt m)^2 = 2/3; larger root"""
    x1, x2 = np.sqrt(m1), np.sqrt(m2)
    s, q = x1 + x2, m1 + m2
    x3 = 2 * s + np.sqrt(6 * s**2 - 3 * q)
    return x3**2


def koide_Q(m):
    m = np.asarray(m, float)
    return m.sum() / np.sqrt(m).sum()**2


def main():
    log, res = [], {}
    say = lambda s: (print(s), log.append(s))
    # ---- KK2
    mu2 = MMU**2 - ME**2
    tau_kk = np.sqrt(ME**2 + 4 * mu2)
    say(f"KK2: train e, mu -> M = {ME} MeV, 1/R2 = {np.sqrt(mu2):.3f} MeV")
    say(f"     exam tau: predicted {tau_kk:.1f} MeV, measured {MTAU} +- {SMTAU}: off by factor "
        f"{MTAU / tau_kk:.2f} ({(MTAU - tau_kk) / SMTAU:.0f} sigma) -> FAIL")
    levels = sorted({abs(n): np.sqrt(ME**2 + n * n * mu2) for n in range(-30, 31)}.items())
    below = sum(1 for n in range(-30, 31) if np.sqrt(ME**2 + n * n * mu2) < 2000)
    say(f"     exam count: charged leptons below 2 GeV predicted {below} (|n2| <= 18, each n2 != 0 twice), "
        f"measured 3 -> FAIL")
    res["KK2"] = dict(params=2, train=["e", "mu"], pred_tau=tau_kk, meas_tau=MTAU, n_below_2GeV=below,
                      verdict="FAIL")
    # ---- KK2 with free level: look-elsewhere
    ratio = (MTAU**2 - ME**2) / mu2
    j = int(round(np.sqrt(ratio)))
    pj = np.sqrt(ME**2 + j * j * mu2)
    say(f"KK2-free: tau level n2 = j free: need j = {np.sqrt(ratio):.3f}; nearest j = {j} -> {pj:.1f} MeV "
        f"({(pj - MTAU) / SMTAU:.0f} sigma) -> FAIL even with the extra discrete choice")
    # second hidden circle: mu at level j1 = p1^2+q1^2, tau at j2 = p2^2+q2^2 (levels up to LMAX)
    hits = {}
    for LMAX in (2000, 20000, 200000):
        P = int(np.sqrt(LMAX)) + 1
        sums = np.unique([p * p + q * q for p in range(P) for q in range(P) if 0 < p * p + q * q <= LMAX])
        sset = set(sums.tolist())
        h = []
        for j1 in sums[sums <= LMAX / 250]:
            lo = (((MTAU - SMTAU)**2 - ME**2) / mu2) * j1
            hi = (((MTAU + SMTAU)**2 - ME**2) / mu2) * j1
            h += [(int(j1), j2) for j2 in range(int(np.ceil(lo)), int(np.floor(hi)) + 1) if j2 in sset]
        hits[LMAX] = h
        say(f"KK2 x second circle, levels = p^2+q^2 <= {LMAX}: {len(h)} discrete choices (j_mu, j_tau) put the "
            f"tau within 1 sigma, e.g. {h[:3]}")
    say("     -> the more hidden levels we allow, the more 'fits': with free levels it is a FIT, not a prediction")
    res["KK2_free"] = dict(j_needed=float(np.sqrt(ratio)), j_nearest=j, pred=pj, hits_within_1sigma={k: len(v) for k, v in hits.items()},
                           examples={k: v[:3] for k, v in hits.items()})
    # ---- Koide family
    res["KOIDE"] = {}
    pt = koide_third(ME, MMU)
    v1, v2 = (pt - MTAU) / SMTAU, (pt - MTAU24) / SMTAU24
    say(f"KOIDE leptons: train e, mu -> predicted tau {pt:.3f} MeV; measured {MTAU} +- {SMTAU} ({v1:+.1f} sigma), "
        f"PDG 2024 {MTAU24} +- {SMTAU24} ({v2:+.1f} sigma) -> {'PASS' if abs(v1) < 3 else 'FAIL'} (pre-reg 3 sigma)")
    say(f"     Q(measured) = {koide_Q([ME, MMU, MTAU]):.6f} (2/3 = {2/3:.6f})")
    res["KOIDE"]["leptons"] = dict(params=2, pred=pt, meas=MTAU, sigma=v1, sigma_pdg24=v2,
                                   verdict="PASS" if abs(v1) < 3 else "FAIL")
    for sec, (a, b, c) in QUARKS.items():
        p3 = koide_third(a, b)
        err = p3 / c - 1
        say(f"KOIDE {sec}-type quarks: train {a}, {b} MeV -> predicted heaviest {p3/1000:.2f} GeV, measured "
            f"{c/1000:.2f} GeV ({100*err:+.0f} %) -> {'PASS' if abs(err) < 0.10 else 'FAIL'}; Q(measured) = "
            f"{koide_Q([a, b, c]):.3f}")
        res["KOIDE"][sec] = dict(params=2, pred=p3, meas=c, rel_err=err, Q=koide_Q([a, b, c]),
                                 verdict="PASS" if abs(err) < 0.10 else "FAIL")
    # ---- bookkeeping
    res["ledger"] = [
        dict(family="KK2 (2-torus, levels 0,1,2)", params=2, trained_on=2, exams=2, passed=0),
        dict(family="KK2 with free levels", params="2 + discrete", trained_on=3, exams=0, passed=0),
        dict(family="Koide circle (amplitude sqrt2 fixed post hoc)", params="2 per sector (+1 hidden)",
             trained_on=6, exams=3, passed=1),
    ]
    with open(os.path.join(OUT, "geo_real.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1, default=float)
    with open(os.path.join(OUT, "geo_real.log"), "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
