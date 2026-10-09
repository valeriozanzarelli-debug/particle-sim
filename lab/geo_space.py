"""geo_space.py - Part A exams: read space (and the hidden circle) off the balls alone.

usage: python geo_space.py vac|quench SHAPE      SHAPE in chain48 ring48 torus16 cyl36x18 stadium30x16
Rules frozen in lab/out/geo_prereg.json (geo_prereg.py). The reconstruction sees ONLY the mutual
information matrix (vacuum) or the arrival-time matrix (quench). Grading uses the drawn shape and the
causal graph of the rule (who can reach whom in one tick) - never passed to the reconstruction.
Writes lab/out/geo_space_<mode>_<shape>.json and .npz (for the figure).
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import geo_core as g

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
PRE = json.load(open(os.path.join(OUT, "geo_prereg.json"), encoding="utf-8"))
THETA = PRE["partA"]["theta"]

SHAPES = {
    "chain48": lambda: g.chain(48),
    "ring48": lambda: g.ring(48),
    "torus16": lambda: g.torus(16, 16),
    "cyl36x18": lambda: g.cylinder(36, 18),
    "stadium30x16": lambda: g.stadium_with_hole(30, 16),
}


def main(mode, key):
    t0 = time.time()
    S = SHAPES[key]()
    ex = PRE["partA"]["exams"][key]
    U = g.rule_unitary(S, THETA)
    D_drawn = S.true_dist()
    D_causal = g.causal_distance(U, S.Ns)
    maxdim = 2 if S.axes == (0, 1) else 1
    if mode == "vac":
        C, gap, nocc = g.vacuum(U)
        I = g.mutual_information(C, S.Ns)
        npieces, rs = g.reconstruct_from_I(I, maxdim=maxdim)
        extra = {"gap_min_abs_sin_omega": gap, "filled": nocc, "modes": 4 * S.Ns}
        raw = I
    else:
        A = g.arrival_times(U, S.Ns, ex["T"], seeds=(0, 1))
        npieces, rs = g.reconstruct_from_arrival(A, maxdim=maxdim)
        extra = {"T": ex["T"], "arrived_fraction": float(np.mean(np.isfinite(A)))}
        raw = A
    pieces = []
    for r in rs:
        idx = r["sites"]
        fd = g.fidelity(r["G"], D_drawn[np.ix_(idx, idx)])
        fc = g.fidelity(r["G"], D_causal[np.ix_(idx, idx)])
        H1 = r["bars"][1]
        pers = H1[:, 1] - H1[:, 0]
        o = np.argsort(-pers)[:4]
        pieces.append(dict(n=len(idx), betti=r["betti"], D_small=r["dim"], D_large=r["dim_large"],
                           h=r["h"], diam=r["diam"], R_est=r["R"],
                           R_over_diam=(r["R"] / r["diam"]) if np.isfinite(r["R"]) else None,
                           big_hole_death=r["big_death"],
                           pearson_drawn=fd[0], spearman_drawn=fd[1],
                           pearson_causal=fc[0], spearman_causal=fc[1],
                           top_H1_bars=H1[o].tolist(), mds_eigs=r["ev"][:6].tolist()))
    # ---- grading against the pre-registered references
    dv = ex["doubled_vacuum" if mode == "vac" else "doubled_quench"]
    it = ex["intended"]
    p0 = pieces[0]
    grade = {
        "pieces_vs_doubled": npieces == dv["pieces"],
        "pieces_vs_drawn": npieces == it["pieces"],
        "b01_vs_doubled": all(p["betti"][:2] == dv["b"][:2] for p in pieces),
        "b01_vs_drawn": all(p["betti"][:2] == it["b"][:2] for p in pieces),
        "fidelity_causal>=0.9": all(p["spearman_causal"] >= 0.9 for p in pieces),
        "dimension": all(round(p["D_small"]) == (it.get("D") or it.get("D_small")) for p in pieces),
    }
    if "D_large" in it:
        grade["D_large=1"] = all(round(p["D_large"]) == 1 for p in pieces)
    if "R_over_extent" in it:
        errs = [p["R_over_diam"] / it["R_over_extent"] - 1 if p["R_over_diam"] else None for p in pieces]
        grade["radius_err"] = errs
        grade["radius<15%"] = all(e is not None and abs(e) < 0.15 for e in errs)
    res = dict(mode=mode, shape=key, name=S.name, Ns=S.Ns, pieces_found=npieces, pieces=pieces,
               grade=grade, seconds=time.time() - t0, **extra)
    tag = f"geo_space_{mode}_{key}"
    with open(os.path.join(OUT, tag + ".json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1, default=float)
    np.savez_compressed(os.path.join(OUT, tag + ".npz"), xy=S.xy, X0=rs[0]["X"], sites0=rs[0]["sites"],
                        G0=rs[0]["G"], Dc=D_causal, Dd=D_drawn, raw=raw)
    print(f"[{mode}] {S.name}: Ns={S.Ns} pieces={npieces} ({extra})")
    for p in pieces[:4]:
        print(f"   n={p['n']} betti={p['betti']} D_small={p['D_small']:.2f} D_large={p['D_large']:.2f} "
              f"spearman drawn={p['spearman_drawn']:.3f} causal={p['spearman_causal']:.3f} "
              f"R/diam={p['R_over_diam']} H1 top={np.round(p['top_H1_bars'][:2], 2).tolist()}")
    print("   grade:", grade, f"({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
