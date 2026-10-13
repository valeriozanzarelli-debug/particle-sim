"""agg2s_analyze.py - exams of "2s all'agglomeratore" against lab/out/agg2s_prereg.json (commit 03dab0d).
Reads lab/out/agg2s_atom.json, agg2s_contact.json, agg2s_post.json (optional), and read-only the aggfull_ results
(first-order dV, learned Darwin coefficient). Writes lab/out/agg2s_verdict.json."""
import os, json
import numpy as np
import agg2s_core as c
from zz_core import dirac_formula

OUT = c.OUT
ld = lambda f: json.load(open(os.path.join(OUT, f))) if os.path.exists(os.path.join(OUT, f)) else {}
P, A, CT, PO = ld("agg2s_prereg.json"), ld("agg2s_atom.json"), ld("agg2s_contact.json"), ld("agg2s_post.json")
AK, LE = ld("aggfull_atomK.json"), ld("aggfull_learn.json")
MC2, RY, AL = c.MC2, c.RY, c.ALPHA
for _j, _v in PO.items():                     # observables recomputed on the twisted grid (agg2s_post.py) replace run-time ones
    if _j in A:
        A[_j]["k4_c_a0"], A[_j]["kin_remainder_eV"] = _v["k4_c_a0"], _v["kin_remainder_eV"]
BC = c.born_constants()
GK = {l["B"]: l["physics"]["gk"][0] for l in LE["levels"] if l.get("chosen")}
dV = dict(s1=float(np.mean(AK["qK8_1s"]["dV_eigen_eV"])), s2=float(np.mean(AK["qK8_2s"]["dV_eigen_eV"])),
          p_lo=float(np.mean(AK["qK8_2p"]["dV_eigen_eV"][:2])), p_hi=float(np.mean(AK["qK8_2p"]["dV_eigen_eV"][2:])),
          s1_a16=float(np.mean(AK["qK16_1s"]["dV_eigen_eV"])))
dirac = lambda n, k: (dirac_formula(n, k) - 1) * MC2
A_ex = lambda n: dirac(n, -1) + RY / n ** 2                     # alpha^4 (and higher) shift of ns: Dirac - Schroedinger
out = dict(born_constants=BC, dV=dV, exams={}, numbers={})
ex = out["exams"]
num = out["numbers"]


def check(name, ok, detail):
    ex[name] = dict(passed=bool(ok), detail=detail)
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")


def have(*jobs):
    return all(A.get(j, {}).get("converged") for j in jobs)


def avg(pfx, key="E_eV"):
    return 0.5 * (A[pfx + "_000"][key] + A[pfx + "_111"][key])


rel = lambda E, Eref: E / Eref - 1                                # > 0: more bound than the reference


def dEB(rho_c, h_a0, m=1.0):
    """Born cell energy of the nucleus cell: -2 m (4 pi Z alpha)^2 rho_c J (Hartree, a.u.) -> eV"""
    return -2 * m * (4 * np.pi) ** 2 * rho_c * BC["J_coef"] * h_a0 ** 3 * 2 * RY


def bornrho(rho_c, h_a0, m=1.0):
    return rho_c / (1 - 8 * np.pi * BC["I_coef"] * m * h_a0) ** 2


def bornk4(k4_c, rho_B, h_a0, m=1.0):
    return k4_c + (8 * np.pi * m) ** 2 * rho_B * BC["I_coef"] * h_a0


# ------------------------------------------------------------------ box
if have("K2s_000", "K2s_111"):
    ED = dirac(2, -1)
    sh = rel(A["K2s_000"]["E_eV"], ED) - rel(avg("K2s"), ED)
    pr = P["predictions_from_development"]["box_shift_2s_L32_rel"]
    num["box_2s_K"] = sh
    check("BOX", abs(sh / pr - 1) <= 0.15, f"2s box shift (K, L = 32 a0) {sh:+.4e} vs predicted {pr:+.4e} ({100 * (sh / pr - 1):+.1f} %)")
if have("K2p_000", "K2p_111"):
    ED = dirac(2, 1)
    sh = rel(A["K2p_000"]["E_eV"], ED) - rel(avg("K2p"), ED)
    pr = P["predictions_from_development"]["box_shift_2p_L32_rel"]
    num["box_2p_K"] = sh
    check("BOX2p", abs(sh / pr - 1) <= 0.20, f"2p box shift (K) {sh:+.4e} vs predicted {pr:+.4e} ({100 * (sh / pr - 1):+.1f} %)")
if have("K2s_000", "K2s_100", "K2s_110", "K2s_111"):
    E4 = [A[f"K2s_{t}"]["E_eV"] for t in ("000", "100", "110", "111")]
    num["K2s_band_eV"] = E4
    num["K2s_full_twist_avg_eV"] = (E4[0] + 3 * E4[1] + 3 * E4[2] + E4[3]) / 8
    num["K2s_pair_minus_full_rel"] = rel(0.5 * (E4[0] + E4[3]), dirac(2, -1)) - rel(num["K2s_full_twist_avg_eV"], dirac(2, -1))
    num["K2s_band_steps_eV"] = list(np.diff(E4))


# ------------------------------------------------------------------ E2X / E2XT / E2Xs
def levels(X):
    L = {}
    if have(f"{X}1s_000", f"{X}1s_111"):
        L["1s1/2"] = (avg(f"{X}1s") + dV["s1"], dirac(1, -1))
    if have(f"{X}2s_000", f"{X}2s_111"):
        L["2s1/2"] = (avg(f"{X}2s") + dV["s2"], dirac(2, -1))
    if have(f"{X}2p_000", f"{X}2p_111"):
        L["2p1/2"] = (avg(f"{X}2p") + dV["p_lo"], dirac(2, 1))
        L["2p3/2"] = (avg(f"{X}2p") + dV["p_hi"], dirac(2, -2))
    return L


for X, nm in (("K", "E2X"), ("T", "E2XT")):
    L = levels(X)
    num[f"levels_{X}"] = {k: dict(E_eV=v[0], dirac_eV=v[1], rel=rel(*v)) for k, v in L.items()}
    if len(L) == 4:
        worst = max(abs(rel(*v)) for v in L.values())
        check(nm, worst <= 1e-3, "; ".join(f"{k} {v[0]:.6f} eV ({100 * rel(*v):+.4f} %)" for k, v in L.items()))
    if X == "K" and "2s1/2" in L:
        r2 = rel(*L["2s1/2"])
        check("E2Xs", abs(r2) <= 5e-4, f"2s1/2 (K, twist avg + dV) {L['2s1/2'][0]:.7f} eV vs Dirac {L['2s1/2'][1]:.7f} ({100 * r2:+.4f} %)")
    # periodic-only, for comparison (= the aggfull numbers with this code)
    num[f"levels_{X}_periodic"] = {k: rel(A[f"{X}{s}_000"]["E_eV"] + d, D) for k, s, d, D in
                                  (("1s1/2", "1s", dV["s1"], dirac(1, -1)), ("2s1/2", "2s", dV["s2"], dirac(2, -1)),
                                   ("2p1/2", "2p", dV["p_lo"], dirac(2, 1)), ("2p3/2", "2p", dV["p_hi"], dirac(2, -2)))
                                  if A.get(f"{X}{s}_000", {}).get("converged")}

# ------------------------------------------------------------------ budget of the 2s (K)
if have("K2s_000", "K2s_111", "S2s_000", "S2s_111"):
    ED = dirac(2, -1)
    k = A["K2s_000"]
    mass = k["m_kin"] - 1
    cell = rel(avg("S2s"), -RY / 4)
    krem = avg("K2s", "kin_remainder_eV")
    a4 = (krem + dV["s2"] - A_ex(2)) / ED
    box = rel(k["E_eV"], ED) - rel(avg("K2s"), ED)
    obs0 = rel(k["E_eV"] + dV["s2"], ED)
    obsT = rel(avg("K2s") + dV["s2"], ED)
    pred0, predT = mass + cell + a4 + box, mass + cell + a4
    bud = dict(mass=mass, cell=cell, alpha4=a4, box=box, observed_periodic=obs0, predicted_periodic=pred0,
               observed_twisted=obsT, predicted_twisted=predT)
    if PO.get("K2s_000") and PO.get("K2s_111"):
        be = 0.5 * (PO["K2s_000"]["band_err_eV"] + PO["K2s_111"]["band_err_eV"])
        kd = 0.5 * (PO["K2s_000"]["krem_D_eV"] + PO["K2s_111"]["krem_D_eV"])
        bud["alpha4_band_shape"] = be / ED
        bud["alpha4_contact"] = (kd + dV["s2"] - A_ex(2)) / ED
    num["budget_2s_K"] = bud
    check("BUD", abs(pred0 - obs0) <= 3e-5 and abs(predT - obsT) <= 3e-5,
          f"periodic: observed {obs0:+.5e} vs sum {pred0:+.5e} (mass {mass:+.3e}, cell {cell:+.3e}, alpha4 {a4:+.3e}, box {box:+.3e}); "
          f"twisted: {obsT:+.5e} vs {predT:+.5e}")

# ------------------------------------------------------------------ E2B (Born cell energy)
rows, okB = [], True
for X in ("K", "T"):
    for s, n, d in (("1s", 1, dV["s1"]), ("2s", 2, dV["s2"])):
        if not have(f"{X}{s}_000", f"{X}{s}_111"):
            continue
        a = A[f"{X}{s}_000"]
        e = avg(f"{X}{s}") + d + dEB(avg(f"{X}{s}", "rho_c_a0"), a["h_a0"], a["m_kin"])
        r = e / (a["m_kin"] * dirac(n, -1)) - 1
        num[f"E2B_{X}{s}"] = dict(E_eV=e, dE_B_eV=dEB(avg(f"{X}{s}", "rho_c_a0"), a["h_a0"], a["m_kin"]), rel_to_mkin_dirac=r,
                                  rel_without_born=(avg(f"{X}{s}") + d) / (a["m_kin"] * dirac(n, -1)) - 1)
        okB &= abs(r) <= 5e-5
        rows.append(f"{X}{s} {100 * r:+.5f} % (without Born {100 * num[f'E2B_{X}{s}']['rel_without_born']:+.5f} %)")
if len(rows) == 4:
    check("E2B", okB, "; ".join(rows) + " (vs m_kin x Dirac)")

# ------------------------------------------------------------------ BORN_S (exact Schroedinger at a0/8)
rows, okS = [], True
for s, n in (("1s", 1), ("2s", 2)):
    if not have(f"S{s}_000", f"S{s}_111"):
        continue
    h_a0 = A[f"S{s}_000"]["h_a0"]
    rc, k4 = avg(f"S{s}", "rho_c_a0"), avg(f"S{s}", "k4_c_a0")
    rB = bornrho(rc, h_a0)
    kB = bornk4(k4, rB, h_a0)
    rex, kex = 1 / (np.pi * n ** 3), 8 / n ** 3 - 3 / n ** 4
    E = avg(f"S{s}") + dEB(rc, h_a0)
    re = rel(E, -RY / n ** 2)
    num[f"BORN_S_{s}"] = dict(rho_c_ratio=rc / rex, rho_B_ratio=rB / rex, k4_c_ratio=k4 / kex, k4_B_ratio=kB / kex,
                              E_rel=rel(avg(f"S{s}"), -RY / n ** 2), E_born_rel=re)
    okS &= abs(rB / rex - 1) <= 5e-3 and abs(kB / kex - 1) <= 5e-3 and abs(re) <= 2e-5
    rows.append(f"{s}: rho {rc / rex:.5f} -> {rB / rex:.5f}, <k^4> {k4 / kex:.5f} -> {kB / kex:.5f}, E {100 * rel(avg(f'S{s}'), -RY / n ** 2):+.5f} % -> {100 * re:+.5f} %")
if len(rows) == 2:
    check("BORN_S", okS, "; ".join(rows))

# ------------------------------------------------------------------ BORN_3s (blind, contact runs a0/4, L = 48)
if all(CT.get(f"S{s}_{t}", {}).get("converged") for s in ("1s", "2s", "3s") for t in ("000", "111")):
    g = lambda s, k: 0.5 * (CT[f"S{s}_000"][k] + CT[f"S{s}_111"][k])
    lost = {}
    for s in ("1s", "2s", "3s"):
        n = int(s[0])
        h_a0 = CT[f"S{s}_000"]["h_a0"]
        rc, k4 = g(s, "rho_c_a0"), g(s, "k4_c_a0")
        rB = bornrho(rc, h_a0)
        kB = bornk4(k4, rB, h_a0)
        rex, kex = 1 / (np.pi * n ** 3), 8 / n ** 3 - 3 / n ** 4
        E = g(s, "E_eV")
        num[f"BORN_a4_{s}"] = dict(rho_c_ratio=rc / rex, rho_B_ratio=rB / rex, k4_c_ratio=k4 / kex, k4_B_ratio=kB / kex,
                                   E_rel=rel(E, -RY / n ** 2), E_born_rel=rel(E + dEB(rc, h_a0), -RY / n ** 2),
                                   E_rel_000=CT[f"S{s}_000"]["E_rel"], E_rel_111=CT[f"S{s}_111"]["E_rel"])
        lost[s] = (rex - rc) / rc
    num["contact_lost_per_unit_rho_c"] = lost
    b3 = num["BORN_a4_3s"]
    uni = max(lost.values()) / min(lost.values()) - 1
    check("BORN_3s", abs(b3["rho_B_ratio"] - 1) <= 0.02 and abs(b3["k4_B_ratio"] - 1) <= 0.02 and abs(b3["E_born_rel"]) <= 2e-4 and uni <= 0.05,
          f"3s: rho {b3['rho_c_ratio']:.4f} -> {b3['rho_B_ratio']:.4f}, <k^4> {b3['k4_c_ratio']:.4f} -> {b3['k4_B_ratio']:.4f}, "
          f"E {100 * b3['E_rel']:+.4f} % -> {100 * b3['E_born_rel']:+.4f} %; lost contact per unit rho_c 1s/2s/3s "
          + " / ".join(f"{v:.4f}" for v in lost.values()) + f" (spread {100 * uni:.1f} %)")

# ------------------------------------------------------------------ F3S (contact rule of the nucleus cell)
def f3s(X):
    rows = {}
    cases = [("1s a0/16", f"{X}1s_000_a16", None, 1, 256), ("1s a0/8", f"{X}1s_000", f"{X}1s_111", 1, 512),
             ("2s a0/8", f"{X}2s_000", f"{X}2s_111", 2, 512)]
    for nm, j0, j1, n, B in cases:
        if not have(j0) or (j1 and not have(j1)):
            continue
        g = (lambda k: 0.5 * (A[j0][k] + A[j1][k])) if j1 else (lambda k: A[j0][k])
        a = A[j0]
        m, c4, h_a0 = a["m_kin"], a["c4"], a["h_a0"]
        d = -GK[B]
        rc, k4 = g("rho_c_a0"), g("k4_c_a0")
        rB = bornrho(rc, h_a0, m)
        kB = bornk4(k4, rB, h_a0, m)
        MV, D = c4 * AL ** 4 * kB * MC2, 4 * np.pi * d * AL ** 4 * rB * MC2
        MVc, Dc = c4 * AL ** 4 * k4 * MC2, 4 * np.pi * d * AL ** 4 * rc * MC2
        dv = dV["s1_a16"] if B == 256 else (dV["s1"] if n == 1 else dV["s2"])
        Aagg = g("kin_remainder_eV") + dv
        Ae = A_ex(n)
        rows[nm] = dict(A_rule_eV=MV + D, A_plain_eV=MVc + Dc, A_aggfull_style_eV=Aagg, A_exact_eV=Ae,
                        err_rule=(MV + D) / Ae - 1, err_plain=(MVc + Dc) / Ae - 1, err_aggfull_style=Aagg / Ae - 1,
                        MV_rule_eV=MV, D_rule_eV=D, MV_exact_eV=-(8 / n ** 3 - 3 / n ** 4) / 8 * AL ** 4 * MC2,
                        D_exact_eV=AL ** 4 * MC2 / (2 * n ** 3), c4=c4, d=d)
    return rows


for X in ("T", "K"):
    R = f3s(X)
    num[f"F3S_{X}"] = R
    if len(R) == 3:
        worst = max(abs(v["err_rule"]) for v in R.values())
        det = "; ".join(f"{k}: {100 * v['err_rule']:+.2f} % (plain {100 * v['err_plain']:+.1f} %, aggfull style "
                        f"{100 * v['err_aggfull_style']:+.1f} %)" for k, v in R.items())
        check(f"F3S_{X}", worst <= 0.05, det)
        check(f"F3Sstrong_{X}", worst <= 0.015, f"worst {100 * worst:.2f} %")

json.dump(out, open(os.path.join(OUT, "agg2s_verdict.json"), "w"), indent=1, default=float)
print(f"\n{sum(v['passed'] for v in ex.values())} PASS / {sum(not v['passed'] for v in ex.values())} FAIL")
