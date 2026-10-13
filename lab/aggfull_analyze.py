"""aggfull_analyze.py - verdicts of AGGLOMERATORE COMPLETO against lab/out/aggfull_prereg.json.
Reads aggfull_learn.json, aggfull_atom.json, aggfull_fineref.json, aggfull_proton.pkl; writes aggfull_verdict.json."""
import json, os, pickle
import numpy as np
import aggfull_core as ac
import aggfull_atom as aa
from zz_core import dirac_formula

OUT = ac.OUT
L = json.load(open(os.path.join(OUT, "aggfull_learn.json")))
R = json.load(open(os.path.join(OUT, "aggfull_atom.json")))
FR = json.load(open(os.path.join(OUT, "aggfull_fineref.json")))
MC2, ALPHA, HEV = ac.MC2, ac.zz.ALPHA, 4.135667696e-15
V, rep = [], {}


def check(name, ok, msg):
    V.append((name, bool(ok), msg))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def dirac_eV(n, kappa):
    return -(1 - dirac_formula(n, kappa)) * MC2


# ---- L1-L3
for lv in L["levels"]:
    if lv.get("chosen") is None:
        print(f"  B={lv['B']}: no admissible coupling (flag): {[(c['stencil'], round(c['tame_max'], 3)) for c in lv['candidates']]}")
        continue
    B, ph = lv["B"], lv["physics"]
    check(f"L1 B={B} coupling vs fine walk", lv["vs_fine_rel"] <= 0.02 and lv["vs_fine_abs"] <= 1e-4,
          f"rel(k-dep) {lv['vs_fine_rel']:.2e}, abs {lv['vs_fine_abs']:.2e} ({lv['chosen']}, {lv['n_sites']} sites)")
    gk, gx, gso = ph["gk"][0], ph["gx"][0], ph["gso"][0]
    check(f"L2 B={B} Darwin", abs(gk / -0.125 - 1) <= 0.02 and abs(gx / 0.25 - 1) <= 0.02,
          f"gk {gk:+.5f} ({100 * (gk / -0.125 - 1):+.2f} %), gx {gx:+.5f} ({100 * (gx / 0.25 - 1):+.2f} %)")
    check(f"L3 B={B} spin-orbit", abs(gso / 0.25 - 1) <= 0.02, f"gso {gso:+.5f} ({100 * (gso / 0.25 - 1):+.2f} %)")

# ---- levels
def level(job, qjob, kind):
    E0 = R[job]["E_eV"]
    q = R.get(qjob)
    if q is None:
        return None
    ev = np.array(q["dV_eigen_eV"])
    if kind == "2p":
        return dict(E0=E0, half=E0 + ev[:2].mean(), three=E0 + ev[2:].mean(), R=q["kinetic_remainder_eV"], ev=ev)
    return dict(E0=E0, half=E0 + ev.mean(), R=q["kinetic_remainder_eV"], ev=ev)


lev = {}
for job, qj, kind, nm in (("e16_1s", "q16_1s", "1s", "1s_a16"), ("e8_1s", "q8_1s", "1s", "1s_a8"),
                          ("e8_2s", "q8_2s", "2s", "2s_a8"), ("e8_2p", "q8_2p", "2p", "2p_a8")):
    if job in R and qj in R:
        lev[nm] = level(job, qj, kind)
rep["levels"] = {k: {kk: (vv.tolist() if isinstance(vv, np.ndarray) else vv) for kk, vv in v.items()} for k, v in lev.items()}
if "1s_a16" in lev:
    e = lev["1s_a16"]["half"]
    d = dirac_eV(1, -1)
    fr = -FR["levels"]["1s1/2"]["E_bind_eV"]
    check("E1 1s a0/16 vs Dirac", abs(e / d - 1) <= 1e-3, f"{e:.6f} eV vs {d:.6f} ({100 * (e / d - 1):+.4f} %); vs fine radial walk "
          f"{fr:.6f} ({100 * (e / fr - 1):+.4f} %)")
tab = []
if all(k in lev for k in ("1s_a8", "2s_a8", "2p_a8")):
    items = [("1s1/2", lev["1s_a8"]["half"], dirac_eV(1, -1)), ("2s1/2", lev["2s_a8"]["half"], dirac_eV(2, -1)),
             ("2p1/2", lev["2p_a8"]["half"], dirac_eV(2, 1)), ("2p3/2", lev["2p_a8"]["three"], dirac_eV(2, -2))]
    worst = max(abs(e / d - 1) for _, e, d in items)
    check("E2 four levels a0/8 vs Dirac", worst <= 1e-3, "; ".join(f"{n} {e:.6f} ({100 * (e / d - 1):+.4f} %)" for n, e, d in items))
    e, d = items[1][1], items[1][2]
    check("S1 2s a0/8 vs Dirac (1 %)", abs(e / d - 1) <= 1e-2, f"{e:.6f} eV vs {d:.6f} ({100 * (e / d - 1):+.4f} %)")
    q = R["q8_2p"]
    fs = q["splitting_eV"] / HEV
    fsd = (dirac_eV(2, -2) - dirac_eV(2, 1)) / HEV
    check("F1 fine structure 2p vs Dirac", abs(fs / fsd - 1) <= 0.01, f"{fs / 1e9:.5f} GHz vs {fsd / 1e9:.5f} GHz ({100 * (fs / fsd - 1):+.3f} %); "
          f"fine walk {FR['fs_2p_eV'] / HEV / 1e9:.5f} GHz")
    spl = q["splitting_eV"]
    check("F2 2+4 pattern and Lande interval", q["lower_pair_spread_eV"] <= 0.01 * spl and q["upper_quartet_spread_eV"] <= 0.01 * spl
          and abs(q["interval_ratio"] / 0.5 - 1) <= 0.01,
          f"pair spread {q['lower_pair_spread_eV'] / spl:.1e}, quartet spread {q['upper_quartet_spread_eV'] / spl:.1e} of the "
          f"splitting; interval ratio {q['interval_ratio']:.5f} (0.5)")
    # alpha^4 shifts
    sh = []
    for nm, key, part, n, kap in (("1s (a0/16)", "1s_a16", "half", 1, -1), ("1s", "1s_a8", "half", 1, -1), ("2s", "2s_a8", "half", 2, -1),
                                  ("2p1/2", "2p_a8", "half", 2, 1), ("2p3/2", "2p_a8", "three", 2, -2)):
        if key not in lev:
            continue
        x = lev[key]
        ours = x["R"] + (x[part] - x["E0"])
        dirac = dirac_eV(n, kap) - (-ac.RY / n ** 2)
        sh.append((nm, ours, dirac))
    worst = max(abs(o / d - 1) for _, o, d in sh)
    check("F3 alpha^4 shifts vs Dirac - Schroedinger", worst <= 0.05,
          "; ".join(f"{n}: {o:.4e} vs {d:.4e} eV ({100 * (o / d - 1):+.1f} %)" for n, o, d in sh))
    rep["alpha4"] = sh

# ---- recoil
st = pickle.load(open(os.path.join(OUT, "aggfull_proton.pkl"), "rb"))
mks = []
for Bp, h in ((2 ** 19, 256 * ac.TH_A), (2 ** 20, 512 * ac.TH_A)):
    r = [l for l in st["levels"] if l["B"] == Bp][0]
    mk, alld = aa.kinetic_mass(r, h, m_ref=ac.MP)
    mks.append((Bp, mk, np.array(alld)))
check("R1 proton kinetic mass", all(abs(m - 1) <= 1e-3 for _, m, _ in mks),
      "; ".join(f"cell {'a0/16' if Bp == 2 ** 19 else 'a0/8'}: m_kin/m_p {m:.6f} (dirs {a[:, 0].round(6).tolist()})" for Bp, m, a in mks))
want = -1 / (1 + ac.MP)
rr = []
for e, p in (("e16_1s", "p16_1s"), ("e8_1s", "p8_1s"), ("e8_2s", "p8_2s"), ("e8_2p", "p8_2p")):
    if e in R and p in R:
        rr.append((e[1:], (R[p]["E_eV"] / R[e]["E_eV"] - 1) / want))
if rr:
    check("R2 recoil shift vs reduced mass", all(abs(x - 1) <= 0.01 for _, x in rr),
          "; ".join(f"{n}: {x:.5f}" for n, x in rr) + f"  (expected shift {want:.4e})")
if "m8_1s" in R and "M_over_me_plus_mp" in R["m8_1s"]:
    m = R["m8_1s"]["M_over_me_plus_mp"]
    check("R3 moving atom mass", abs(m - 1) <= 1e-3, f"M/(m_e + m_p) = {m:.6f}")

# ---- amendment exams (aggfull_amend.json / aggfull_amend2.json): the same exams with the re-learned kinetic
JK = os.path.join(OUT, "aggfull_atomK.json")
if os.path.exists(JK):
    RK = json.load(open(JK))
    g = lambda j: RK.get(j[0] + "K" + j[1:])
    levK = {}
    for job, qj, kind, nm in (("e16_1s", "q16_1s", "1s", "1s_a16"), ("e8_1s", "q8_1s", "1s", "1s_a8"),
                              ("e8_2s", "q8_2s", "2s", "2s_a8"), ("e8_2p", "q8_2p", "2p", "2p_a8")):
        e, q = g(job), g(qj)
        if e and q and e.get("converged"):
            ev = np.array(q["dV_eigen_eV"])
            levK[nm] = dict(E0=e["E_eV"], half=e["E_eV"] + (ev[:2].mean() if kind == "2p" else ev.mean()),
                            three=e["E_eV"] + ev[2:].mean() if kind == "2p" else None, R=q["kinetic_remainder_eV"], ev=ev)
    rep["levels_K"] = {k: {kk: (vv.tolist() if isinstance(vv, np.ndarray) else vv) for kk, vv in v.items()} for k, v in levK.items()}
    if "1s_a16" in levK:
        e, d = levK["1s_a16"]["half"], dirac_eV(1, -1)
        check("E1K 1s a0/16 vs Dirac (re-learned kinetic)", abs(e / d - 1) <= 1e-3, f"{e:.6f} eV vs {d:.6f} ({100 * (e / d - 1):+.4f} %)")
    itemsK = [(n, levK[k][p], dirac_eV(nn, kap)) for n, k, p, nn, kap in (("1s1/2", "1s_a8", "half", 1, -1), ("2s1/2", "2s_a8", "half", 2, -1),
              ("2p1/2", "2p_a8", "half", 2, 1), ("2p3/2", "2p_a8", "three", 2, -2)) if k in levK]
    if len(itemsK) == 4:
        worst = max(abs(e / d - 1) for _, e, d in itemsK)
        check("E2K four levels a0/8 vs Dirac (re-learned kinetic)", worst <= 1e-3,
              "; ".join(f"{n} {e:.6f} ({100 * (e / d - 1):+.4f} %)" for n, e, d in itemsK))
    shK = []
    for nm, key, part, n, kap in (("1s (a0/16)", "1s_a16", "half", 1, -1), ("1s", "1s_a8", "half", 1, -1), ("2s", "2s_a8", "half", 2, -1),
                                  ("2p1/2", "2p_a8", "half", 2, 1), ("2p3/2", "2p_a8", "three", 2, -2)):
        if key in levK:
            x = levK[key]
            shK.append((nm, x["R"] + (x[part] - x["E0"]), dirac_eV(n, kap) - (-ac.RY / n ** 2)))
    if shK:
        worst = max(abs(o / d - 1) for _, o, d in shK)
        check(f"F3K alpha^4 shifts, re-learned kinetic ({len(shK)} of 5 levels)", worst <= 0.05,
              "; ".join(f"{n}: {o:.4e} vs {d:.4e} eV ({100 * (o / d - 1):+.1f} %)" for n, o, d in shK))
        rep["alpha4_K"] = shK
    rrK = []
    for e, p in (("e16_1s", "p16_1s"), ("e8_1s", "p8_1s"), ("e8_2s", "p8_2s"), ("e8_2p", "p8_2p")):
        if g(e) and g(p) and g(e).get("converged") and g(p).get("converged"):
            rrK.append((e[1:], (g(p)["E_eV"] / g(e)["E_eV"] - 1) / want))
    if rrK:
        check(f"R2K recoil, re-learned kinetic ({len(rrK)} of 4)", all(abs(x - 1) <= 0.01 for _, x in rrK),
              "; ".join(f"{n}: {x:.5f}" for n, x in rrK))
# ---- box and old-restriction diagnostics (not exams)
diag = {}
for j in ("e8_1s_N128", "e8_1s_N192", "e8_2s_N192", "ca8_1s_N192", "ca8_2s_N192", "ca8_1s", "ca8_2s", "e8_1s", "e8_2s"):
    if j in R:
        diag[j] = dict(E_eV=R[j]["E_eV"], rel=R[j]["rel_to_schroedinger"], box_a0=R[j]["box_a0"], converged=R[j]["converged"])
rep["diagnostics"] = diag
for j, d in diag.items():
    print(f"  diag {j}: {d['E_eV']:.6f} eV ({100 * d['rel']:+.4f} % vs Schroedinger, box {d['box_a0']:.0f} a0) conv {d['converged']}")
rep["verdict"] = V
json.dump(rep, open(os.path.join(OUT, "aggfull_verdict.json"), "w"), indent=1, default=float)
n = sum(v[1] for v in V)
print(f"=== {n} PASS / {len(V) - n} FAIL")
