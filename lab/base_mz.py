"""Mach–Zehnder with the frozen base rule, real 2-D geometry (beam splitters and mirrors at 45°).

  python lab/base_mz.py phase      single-photon interference vs phase φ (P_D = cos²(φ/2))
  python lab/base_mz.py bomb       Elitzur–Vaidman bomb tester (P_explode 1/2, P_C 1/4, P_D 1/4)
  python lab/base_mz.py wheeler    delayed choice: BS2 inserted / removed after the photon passed BS1

Geometry (x, z), λ = 1:  source at z = −20 → BS1 at (0,0) → arm A (+z) with phase plate at z = 12.5 and
mirror M_A at (0,25) → arm B (+x) with compensating plate at x = 12.5 and mirror M_B at (25,0) → BS2 at
(25,25) → detector C (output travelling +x, screen at x = 40) and D (output travelling +z, screen at z = 40).
BS: t = 1/√2, r = i/√2 (Jones ×1, both faces); mirror: r = diag(−1, +1) (s = E_y, p = H_y); plate: e^{iφ}.
QM: amplitude to D = t·r·e^{iφ} + r·t → P_D = cos²(φ/2), P_C = sin²(φ/2).
"""
import sys
import time
import numpy as np
from base_rule import *
from base_exp2d import report

TOL = dict(p_abs=0.02)            # every probability: |P − P_QM| ≤ max(3σ, 0.02)

LA, LD, ZS = 25.0, 15.0, 20.0
D45 = np.array([1.0, 1.0]) / np.sqrt(2)
BS_T = np.eye(2) / np.sqrt(2)
BS_R = 1j * np.eye(2) / np.sqrt(2)
MIRROR = np.diag([-1.0, 1.0])


def mz_scene(phi, bomb=False, bs2=True):
    sc = Scene2D()
    ids = {}
    ids["source"] = sc.add(Obj2("source", "source", (0, -ZS), (1, 0), 20.0, src=dict(w0=3.0, zw=ZS + LA, jones=(1, 0))))
    ids["BS1"] = sc.add(Obj2("BS1", "bs", (0, 0), D45, 18.0, T=BS_T, R=BS_R))
    ids["plateA"] = sc.add(Obj2("plate A", "phase", (0, LA / 2), (1, 0), 10.0, T=np.exp(1j * phi)))
    ids["plateB"] = sc.add(Obj2("plate B", "phase", (LA / 2, 0), (0, 1), 10.0, T=1.0))
    if bomb:
        ids["bomb"] = sc.add(Obj2("bomb", "stop", (0.7 * LA, 0), (0, 1), 10.0))
    ids["MA"] = sc.add(Obj2("M_A", "mirror", (0, LA), D45, 18.0, R=MIRROR))
    ids["MB"] = sc.add(Obj2("M_B", "mirror", (LA, 0), D45, 18.0, R=MIRROR))
    if bs2:
        ids["BS2"] = sc.add(Obj2("BS2", "bs", (LA, LA), D45, 18.0, T=BS_T, R=BS_R))
    ids["armA"] = sc.add(Obj2("arm A marker", "plane", (0, 0.8 * LA), (1, 0), 9.0))       # passive: records only
    ids["armB"] = sc.add(Obj2("arm B marker", "plane", (0.8 * LA, 0), (0, 1), 9.0))
    ids["C"] = sc.add(Obj2("detector C", "screen", (LA + LD, LA), (0, 1), 12.0))
    ids["D"] = sc.add(Obj2("detector D", "screen", (LA, LA + LD), (1, 0), 12.0))
    return sc, ids


BOX = [-30, 60, -ZS - 2, 60]


def outcome(res, ids):
    st, wh = res["status"], res["where"]
    n = res["n"]
    c = int(((st == CLICK) & (wh == ids["C"])).sum())
    d = int(((st == CLICK) & (wh == ids["D"])).sum())
    other = {k: int(v) for k, v in res["counts"].items() if v}
    mir = int(((st == ABSORBED) & np.isin(wh, [ids["MA"], ids["MB"]])).sum())
    return n, c, d, mir, other


def check(rows, name, k, n, p, extra=""):
    P = k / n
    sig = np.sqrt(max(p * (1 - p), 1 / n) / n)
    report(rows, name, f"{P:.4f}", f"{p:.4f}", f"max(3σ={3 * sig:.4f}, {TOL['p_abs']})",
           abs(P - p) <= max(3 * sig, TOL["p_abs"]), f"{k}/{n} {extra}")


def exp_phase(n=1500, degs=(0, 90, 180, 270)):
    rows, metas = [], []
    for phi in np.radians(degs):
        sc, ids = mz_scene(phi)
        if not metas:
            log("MZ:", balls_summary(sc))
        res = simulate(sc, ids["source"], n, BOX, seed=400 + int(100 * phi), batch=1500, verbose=False)
        nb, c, d, mir, other = outcome(res, ids)
        m = nb - res["src_coin"]
        log(f"  φ={phi:.3f}: C {c}  D {d}  of {m} (lost at mirrors {mir}; {other}) [{res['time']:.0f}s]")
        check(rows, f"Mach–Zehnder φ={np.degrees(phi):.0f}°: P(D)", d, m, np.cos(phi / 2) ** 2, f"(C: {c})")
        metas.append(dict(phi=phi, **{k: v for k, v in dict(photons=nb, clicks=c + d, counts=res["counts"],
                                                              runtime=res["time"] + sc.build_time).items()},
                          balls=balls_summary(sc)))
    save_results("mz_phase_" + "_".join(str(int(d)) for d in degs), rows, metas)
    return rows, metas


def exp_bomb(n=3000):
    rows = []
    sc, ids = mz_scene(0.0, bomb=True)
    log("bomb:", balls_summary(sc))
    res = simulate(sc, ids["source"], n, BOX, seed=501, batch=1500)
    nb, c, d, mir, other = outcome(res, ids)
    m = nb - res["src_coin"]
    boom = int(((res["status"] == ABSORBED) & (res["where"] == ids["bomb"])).sum())
    log(f"  explode {boom}, C {c}, D {d} of {m}; lost at mirrors {mir}; {other}")
    check(rows, "bomb tester: P(explosion)", boom, m, 0.5)
    check(rows, "bomb tester: P(C) = bomb detected without interaction (dark port)", c, m, 0.25)
    check(rows, "bomb tester: P(D)", d, m, 0.25)
    # without the bomb (same φ = 0): the dark port C never clicks
    sc0, ids0 = mz_scene(0.0)
    res0 = simulate(sc0, ids0["source"], n // 3, BOX, seed=502, batch=1500, verbose=False)
    nb0, c0, d0, mir0, _ = outcome(res0, ids0)
    check(rows, "same interferometer without bomb: P(C) (dark port)", c0, nb0 - res0["src_coin"], 0.0)
    meta = dict(photons=nb, counts=res["counts"], balls=balls_summary(sc), runtime=res["time"] + sc.build_time + res0["time"])
    save_results("bomb", rows, meta)
    return rows, meta


def exp_wheeler(n=1500):
    """Each photon: the choice (insert / remove BS2) is made at t_c, when the photon is in the middle of the
    arms (it crossed BS1 at t ≈ 20, it reaches BS2 at t ≈ 70).  'insert': BS2 absent before t_c, present
    after; 'remove': the opposite.  The fans of BS2 exist only for wave balls that reached BS2 while it was
    there; rays crossing BS2's place are blocked only while it is there (base_kernel: presence gating)."""
    rows, metas = [], []
    tc = ZS + LA / 2
    rng = np.random.default_rng(601)
    n_ins = rng.binomial(n, 0.5)
    for mode, nn, init in (("inserted after the photon passed BS1", n_ins, False),
                           ("removed after the photon passed BS1", n - n_ins, True)):
        sc, ids = mz_scene(0.0)
        sc.set_toggle(ids["BS2"], init)
        res = simulate(sc, ids["source"], nn, BOX, seed=602 + int(init), batch=1500,
                       tog_fn=lambda r_, m: np.full(m, tc), verbose=False)
        nb, c, d, mir, other = outcome(res, ids)
        m = nb - res["src_coin"]
        # which arm did each photon take? (x of the photon when it was at the height of the arms, recorded by
        # the time it crossed z = LA/2 or x = LA/2: reconstructed from its final detector when BS2 is absent)
        log(f"  BS2 {mode}: C {c}, D {d} of {m}; lost at mirrors {mir}; {other} [{res['time']:.0f}s]")
        if init:     # removed → which-path: 1/2, 1/2, and the detector tells the arm
            check(rows, f"Wheeler, BS2 {mode}: P(D)", d, m, 0.5, f"(C: {c}) — QM with BS2 absent: 1/2")
            arm = {}
            for i_, pl, _ in res["rec"]:
                arm[i_] = "A" if pl == ids["armA"] else "B"
            clk = np.nonzero(res["status"] == CLICK)[0]
            okc = sum(1 for i_ in clk if arm.get(i_) == ("A" if res["where"][i_] == ids["C"] else "B"))
            # NOT a QM prediction (QM says nothing about an unobserved path): de Broglie–Bohm trajectories
            # cannot cross, so where the two beams cross (BS2's empty place) the photon from arm A turns into the
            # output of arm B ("surrealistic trajectories", Englert–Scully–Süssmann–Walther 1992; Hiley 2000).
            check(rows, f"Wheeler, BS2 {mode}: fraction of photons whose detector is the straight continuation of their arm "
                        f"(C⇔A, D⇔B)", okc, len(clk), 0.0,
                  "— dBB prediction 0 (no crossing of trajectories); QM: no prediction")
        else:        # inserted → interference: φ = 0 → P(D) = 1
            check(rows, f"Wheeler, BS2 {mode}: P(D)", d, m, 1.0, f"(C: {c}) — QM with BS2 present: cos²(0) = 1")
        metas.append(dict(mode=mode, photons=nb, counts=res["counts"], balls=balls_summary(sc), runtime=res["time"] + sc.build_time))
    save_results("wheeler", rows, metas)
    return rows, metas


if __name__ == "__main__":
    what = sys.argv[1:] or ["phase"]
    if what[0] == "phase" and len(what) > 1:
        t = time.time()
        rows, meta = exp_phase(degs=tuple(float(a) for a in what[1:]))
        log(f"== phase {what[1:]}: {'PASS' if all(r['ok'] for r in rows) else 'FAIL'}  [{time.time() - t:.0f}s]")
        what = []
    for w_ in what:
        t = time.time()
        rows, meta = {"phase": exp_phase, "bomb": exp_bomb, "wheeler": exp_wheeler}[w_]()
        log(f"== {w_}: {'PASS' if all(r['ok'] for r in rows) else 'FAIL'}  [{time.time() - t:.0f}s]")
