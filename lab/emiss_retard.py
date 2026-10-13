"""emiss_retard.py - the retardation term (test N2 of lab/out/emiss_prereg2.json).
The light-atom coupling of the runs is dipole (the atom reads the Gaussian-averaged field at the nucleus); in 3-D the
light cell (b = 4 lambda_C at 20 cells/lambda) is larger than the atom (<r> = 2.7 lambda_C), so the lattice cannot
carry the e^{-ik.r} variation across the atom.  Here: the rule's OWN radial states (radial walk, a = 0.025, as P_N1)
in the full golden rule with alpha.eps e^{-ik.r} (same code as emiss_exact.py), giving F_ret(rule) = A_EX/A_LW;
then (A_N1 / W2) x F_ret(rule) vs Dirac A_EX.  Inside check: velocity form at k -> 0 vs length form.
Run: python lab/emiss_retard.py -> lab/out/emiss_retard.json"""
import json, os
import numpy as np
from scipy.special import roots_legendre
from zz_core import ALPHA, T_UNIT
from emiss_core import radial_state, spin_angle, OUT
import emiss_exact as EXA

SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]


def golden(r, wr, G1, F1, G2, F2, k2, w, kphot, nth=40, nph=48):
    """sum over initial m (averaged), final m, 2 polarisations of |<1s| alpha.eps e^{-i kphot z} |2p>|^2 -> rate"""
    u, wu = roots_legendre(nth)
    ph = (np.arange(nph) + 0.5) * 2 * np.pi / nph
    TH, PH = np.meshgrid(np.arccos(u), ph, indexing="ij")
    WA = np.outer(wu, np.full(nph, 2 * np.pi / nph))
    X, Y, Z = np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)
    j2 = abs(k2) - 0.5
    phase = np.exp(-1j * kphot * np.outer(r, Z.ravel())).reshape((len(r),) + Z.shape)
    tot = 0.0
    for mi in np.arange(-j2, j2 + 1):
        Om2u, Om2l = spin_angle(k2, mi, X, Y, Z), spin_angle(-k2, mi, X, Y, Z)
        for mf in (0.5, -0.5):
            Om1u, Om1l = spin_angle(-1, mf, X, Y, Z), spin_angle(1, mf, X, Y, Z)
            for e in (0, 1):
                A1 = np.einsum("a...,ab,b...->...", Om1u.conj(), SIG[e], Om2l)
                A2 = np.einsum("a...,ab,b...->...", Om1l.conj(), SIG[e], Om2u)
                # radial x angular, integrated: sum_r wr [ (i G1 F2) A1 + (-i F1 G2) A2 ] e^{-ikz}
                I1 = np.einsum("r,rtp,tp->", wr * (1j * G1 * F2), phase, A1 * WA)
                I2 = np.einsum("r,rtp,tp->", wr * (-1j * F1 * G2), phase, A2 * WA)
                tot += abs(I1 + I2)**2
    return ALPHA * w / (2 * np.pi) * 4 * np.pi * tot / (2 * j2 + 1)


def rule_states(a=0.025, rmax=40.0):
    Za = 60 * ALPHA
    o1, r, g1, f1 = radial_state(-1, 1, a, Za, rmax)
    o2, r2, g2, f2 = radial_state(-2, 2, a, Za, rmax)
    assert np.allclose(r, r2)
    return r, np.full_like(r, a), g1, f1, g2, f2, o2 - o1


if __name__ == "__main__":
    res = {}
    # 0. the code on Dirac's own functions (must reproduce emiss_exact.json)
    Za = 60 * ALPHA
    rq, wq = EXA.quad_radial(Za)
    E1, G1, F1 = EXA.radial(1, -1, Za, rq)
    E2, G2, F2 = EXA.radial(2, -2, Za, rq)
    n1 = np.sqrt(np.sum(wq * (G1**2 + F1**2))); n2 = np.sqrt(np.sum(wq * (G2**2 + F2**2)))
    G1, F1, G2, F2 = G1 / n1, F1 / n1, G2 / n2, F2 / n2
    wD = E2 - E1
    AEX_D = golden(rq, wq, G1, F1, G2, F2, -2, wD, wD)
    RD = np.sum(wq * (G1 * G2 + F1 * F2) * rq)
    ALW_D = 4 / 3 * ALPHA * wD**3 * RD**2 / 3
    print(f"Dirac with this code: A_LW {ALW_D:.5e}, A_EX {AEX_D:.5e}, F_ret {AEX_D / ALW_D:.5f}")
    res["dirac_check"] = dict(A_LW=ALW_D, A_EX=AEX_D, F=AEX_D / ALW_D)
    # 1. the rule's states
    r, wr, g1, f1, g2, f2, w = rule_states()
    R = np.sum(wr * (g1 * g2 + f1 * f2) * r)
    ALW = 4 / 3 * ALPHA * w**3 * R**2 / 3
    AEX = golden(r, wr, g1, f1, g2, f2, -2, w, w)
    AV0 = golden(r, wr, g1, f1, g2, f2, -2, w, 1e-7)          # velocity form, k -> 0
    F = AEX / ALW
    print(f"rule (a = 0.025): w {w:.6f}, R {R:.5f}, A_LW {ALW:.5e}, A_EX {AEX:.5e}, F_ret {F:.5f}; velocity k->0 / length {AV0 / ALW:.5f}")
    res["rule"] = dict(omega=w, R=R, A_LW=ALW, A_EX=AEX, F_ret=F, vel0_over_len=AV0 / ALW)
    # 2. test N2
    Rv = json.load(open(os.path.join(OUT, "emiss_results.json")))
    A_N1 = Rv["P_N1_numbers"]["A_over_W2"]
    EX = json.load(open(os.path.join(OUT, "emiss_exact.json")))
    AEXD = EX["Z60"]["2p3/2"]["A_EX"]
    val = A_N1 * F
    res["N2"] = dict(A_N1_over_W2=A_N1, value=val, dirac_A_EX=AEXD, ratio=val / AEXD, pass_=bool(abs(val / AEXD - 1) < 0.03 and abs(AV0 / ALW - 1) < 0.02))
    res["N2x"] = dict(F_rule=F, F_dirac=EX["Z60"]["2p3/2"]["ratio_EX_LW"])
    print(f"N2: (A_N1/W2) x F_ret = {val:.5e} ({val / T_UNIT:.4e} s^-1) vs Dirac A_EX {AEXD:.5e} ({100 * (val / AEXD - 1):+.2f} %) "
          f"-> {'PASS' if res['N2']['pass_'] else 'FAIL'};  N2x: F_rule {F:.5f} vs F_Dirac {EX['Z60']['2p3/2']['ratio_EX_LW']:.5f}")
    json.dump(res, open(os.path.join(OUT, "emiss_retard.json"), "w"), indent=1, default=float)
