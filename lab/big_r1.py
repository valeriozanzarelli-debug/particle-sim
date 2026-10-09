"""big_r1.py - R1: the repaired shape reader (H1 persistence / diameter > 0.4, full filtration) on FRESH nulls
(4 classes x 25, 576 sites) and on 3 known shapes. Pre-registered in lab/out/big_prereg_R.json.
usage: python big_r1.py [boot]   ('boot' = the pre-declared alternative, only if the first reader fails)
Writes lab/out/big_r1[_boot].json (~2-3 min).
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import geo_core as G
import big_core as B
from pre_core import LatticeWalk
from pre_p2 import mi_one_mode, random_regular
from pre_p4 import vacuum_bloch

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
NPER = 25


def strip(r):
    return {k: v for k, v in r.items() if k not in ("X", "ev", "sites", "d")}


if __name__ == "__main__":
    boot = len(sys.argv) > 1 and sys.argv[1] == "boot"
    reader = (lambda I: B.read_shape_boot(I)) if boot else (lambda I: B.read_shape(I, with_mds=False))
    key = "holes_boot" if boot else "holes"
    t0 = time.time()
    res = {"reader": "bootstrap" if boot else "plain"}
    # ---------- known shapes
    shapes = {}
    for name, S in (("ring", G.ring(48)), ("line_x_hidden_circle", G.cylinder(36, 18))):
        C, gap, _ = G.vacuum(G.rule_unitary(S, 0.5))
        shapes[name] = G.mutual_information(C, S.Ns)
    shapes["torus_wilson"] = G.mutual_information(vacuum_bloch(16, 0.5, 1.0), 256)
    tp = {}
    for name, I in shapes.items():
        r = reader(I)
        tp[name] = strip(r)
        print(name, tp[name], f"{time.time()-t0:.0f}s", flush=True)
    # hidden-circle bar of the cylinder (prediction, not graded): all H1 bars, geo_core significance 1.5 h
    rr = B.read_shape(shapes["line_x_hidden_circle"], with_mds=False)
    res["cylinder_all_rel"] = rr["rel"]
    nc, rs = G.reconstruct_from_I(shapes["line_x_hidden_circle"], maxdim=1)
    big = max(rs, key=lambda q: len(q["sites"]))
    D1 = big["bars"][1]; pers = D1[:, 1] - D1[:, 0]; sig = pers > 1.5 * big["h"]
    first = np.argmin(np.where(sig, D1[:, 1], np.inf))
    res["cylinder_hidden_bar"] = dict(death=float(D1[first, 1]), rel=float(pers[first] / big["diam"]),
                                      R_est=float(big["R"]), R_true=18 / (2 * np.pi),
                                      geo_core_betti=[int(b) for b in big["betti"]])
    print("hidden bar", res["cylinder_hidden_bar"], flush=True)
    # ---------- base for the nulls: 24x24 box, theta = 0.125, w = 0.3
    n = 24
    wk = LatticeWalk(n, 0.125, periodic=False, w=[np.full((n, n), 0.3)] * 2)
    Ib, gap = B.vacuum_mi(wk)
    res["base_reader_not_graded"] = strip(reader(Ib))
    print("base box", res["base_reader_not_graded"], f"{time.time()-t0:.0f}s", flush=True)
    Ns = n * n
    iu = np.triu_indices(Ns, 1)
    off = Ib[iu]
    pos = off > G.I_FLOOR
    mu, sd, frac = np.log(off[pos]).mean(), np.log(off[pos]).std(), pos.mean()
    res["base"] = dict(frac_above_floor=float(frac), mean_logI=float(mu), std_logI=float(sd), gap=float(gap))
    rng = np.random.default_rng(20261009 + 7777)
    nulls = {k: [] for k in ("shuffle", "noise", "random_graph", "random_state")}
    for k in range(NPER):
        v = off.copy(); rng.shuffle(v)
        I = np.zeros_like(Ib); I[iu] = v; I = I + I.T; np.fill_diagonal(I, np.diag(Ib))
        nulls["shuffle"].append(strip(reader(I)))
        v = np.where(rng.random(len(off)) < frac, np.exp(rng.normal(mu, sd, len(off))), 0.0)
        I = np.zeros_like(Ib); I[iu] = v; I = I + I.T; np.fill_diagonal(I, np.diag(Ib))
        nulls["noise"].append(strip(reader(I)))
        A = random_regular(Ns, 4, rng)
        w, V = np.linalg.eigh(-A)
        Vo = V[:, :Ns // 2]
        nulls["random_graph"].append(strip(reader(mi_one_mode(Vo @ Vo.conj().T))))
        Z = rng.normal(size=(Ns, Ns)) + 1j * rng.normal(size=(Ns, Ns))
        Q, _ = np.linalg.qr(Z)
        Vo = Q[:, :Ns // 2]
        nulls["random_state"].append(strip(reader(mi_one_mode(Vo @ Vo.conj().T))))
        print(k, {c: (L[-1][key], round(L[-1]["rel"][0], 3) if L[-1]["rel"] else 0) for c, L in nulls.items()},
              f"{time.time()-t0:.0f}s", flush=True)
    rates = {}
    for c, L in nulls.items():
        rates[c] = dict(n=len(L), fp=float(np.mean([r[key] >= 1 for r in L])),
                        max_rel=float(max(r["rel"][0] if r["rel"] else 0 for r in L)),
                        median_rel=float(np.median([r["rel"][0] if r["rel"] else 0 for r in L])),
                        mean_pieces=float(np.mean([r["pieces"] for r in L])),
                        mean_size=float(np.mean([r["size"] for r in L])))
    res["rates"] = rates
    res["known_shapes"] = tp
    res["nulls"] = nulls
    ntp = sum(tp[s][key] >= 1 for s in tp)
    res["verdict"] = {"fp_le_5pct_all_classes": bool(all(r["fp"] <= 0.05 for r in rates.values())),
                      "true_positives": f"{ntp}/3",
                      "R1_pass": bool(all(r["fp"] <= 0.05 for r in rates.values()) and ntp == 3)}
    res["time_s"] = time.time() - t0
    print(json.dumps(rates, indent=1), res["verdict"], res["time_s"])
    with open(os.path.join(OUT, "big_r1_boot.json" if boot else "big_r1.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
