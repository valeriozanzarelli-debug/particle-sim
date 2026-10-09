"""pre_p2.py - P2: null tests (false-positive rate) of the geo_core blind MI reconstruction.
Inputs with NO shape: shuffled MI, random-regular-graph vacuum, log-normal noise with matched statistics,
Haar-random Gaussian state. Positive control: vacuum of the rule on a 12x12 torus.
Pre-registered in lab/out/pre_prereg.json. Writes lab/out/pre_p2.json (~1-2 min).
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import geo_core as G

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SEED = 20261009


def h(nu):
    nu = np.clip(nu.real, 1e-300, 1); mu = np.clip(1 - nu, 1e-300, 1)
    return -(nu * np.log(nu) + mu * np.log(mu))


def mi_one_mode(C):
    """MI matrix of a Gaussian fermion state with ONE mode per site"""
    n = len(C)
    d = np.real(np.diag(C))
    S1 = h(d)
    a = d[:, None]; b = d[None, :]; c2 = np.abs(C)**2
    # eigenvalues of [[a, c],[c*, b]]
    tr, det = a + b, a * b - c2
    disc = np.sqrt(np.maximum(tr**2 / 4 - det, 0))
    S2 = h(tr / 2 + disc) + h(tr / 2 - disc)
    I = S1[:, None] + S1[None, :] - S2
    np.fill_diagonal(I, 2 * S1)
    return np.maximum(I, 0)


def random_regular(n, k, rng):
    while True:
        stubs = np.repeat(np.arange(n), k)
        rng.shuffle(stubs)
        e = stubs.reshape(-1, 2)
        if np.any(e[:, 0] == e[:, 1]):
            continue
        A = np.zeros((n, n))
        ok = True
        for u, v in e:
            if A[u, v]:
                ok = False; break
            A[u, v] = A[v, u] = 1
        if ok:
            return A


def b2_capped(Gm, h, frac=0.5, min_pers=1.5):
    """DEVIATION (declared): geo_core.homology's H2 runs the Rips filtration to infinity and does not finish on
    shapeless metrics (> 8 min for 72 points). Same landmarks and persistence rule, but the filtration is cut at
    frac * diameter; a cavity still alive at the cut dies there."""
    from ripser import ripser
    L = G.landmarks(Gm)
    eps = float(np.max(np.min(np.where(np.isfinite(Gm[:, L]), Gm[:, L], 1e9), axis=1)))
    sub = Gm[np.ix_(L, L)]
    cut = frac * float(np.max(sub[np.isfinite(sub)]))
    D2 = ripser(sub, distance_matrix=True, maxdim=2, thresh=cut)["dgms"][2]
    if len(D2) == 0:
        return 0
    death = np.where(np.isfinite(D2[:, 1]), D2[:, 1], cut)
    return int(((death - D2[:, 0]) > min_pers * max(h, 2 * eps)).sum())


def summary(I):
    nc, rs = G.reconstruct_from_I(I, maxdim=1)
    big = max(rs, key=lambda r: len(r["sites"]))
    b2 = b2_capped(big["G"], big["scale"])
    return dict(pieces=int(nc), size=int(len(big["sites"])), betti=[int(b) for b in big["betti"]] + [b2],
                D=float(big["dim"]), D_large=float(big["dim_large"]),
                # POST-HOC diagnostic: largest H1 persistence relative to the geodesic diameter
                p1_rel=float(np.max(np.r_[0, big["bars"][1][:, 1] - big["bars"][1][:, 0]]) / big["diam"]),
                p1_rel_2nd=float(np.sort(np.r_[0, 0, big["bars"][1][:, 1] - big["bars"][1][:, 0]])[-2] / big["diam"]))


if __name__ == "__main__":
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    res = {}
    shp = G.torus(12, 12)
    U = G.rule_unitary(shp, 0.5)
    C, gap, nocc = G.vacuum(U)
    Iv = G.mutual_information(C, shp.Ns)
    res["positive_control_torus12"] = summary(Iv)
    nc, rs = G.reconstruct_from_I(Iv)
    res["positive_control_torus12"]["all_pieces_betti_full_geo_core"] = [[int(b) for b in r["betti"]] for r in rs]
    print("control", res["positive_control_torus12"], flush=True)
    # one sheet only (72 sites) as the matched base for the nulls
    A_idx = rs[0]["sites"]
    I0 = Iv[np.ix_(A_idx, A_idx)]
    n0 = len(A_idx)
    iu = np.triu_indices(n0, 1)
    off = I0[iu]
    pos = off > G.I_FLOOR
    mu, sd, frac = np.log(off[pos]).mean(), np.log(off[pos]).std(), pos.mean()
    res["base"] = dict(sites=int(n0), frac_above_floor=float(frac), mean_logI=float(mu), std_logI=float(sd))
    nulls = {"shuffle": [], "noise": [], "random_graph": [], "random_state": []}
    for k in range(20):
        v = off.copy(); rng.shuffle(v)
        I = np.zeros_like(I0); I[iu] = v; I = I + I.T; np.fill_diagonal(I, np.diag(I0))
        nulls["shuffle"].append(summary(I))
        v = np.where(rng.random(len(off)) < frac, np.exp(rng.normal(mu, sd, len(off))), 0.0)
        I = np.zeros_like(I0); I[iu] = v; I = I + I.T; np.fill_diagonal(I, np.diag(I0))
        nulls["noise"].append(summary(I))
        A = random_regular(144, 4, rng)
        w, V = np.linalg.eigh(-A)
        Vo = V[:, :72]
        nulls["random_graph"].append(summary(mi_one_mode(Vo @ Vo.conj().T)))
        if k < 10:
            Z = rng.normal(size=(144, 144)) + 1j * rng.normal(size=(144, 144))
            Q, _ = np.linalg.qr(Z)
            Vo = Q[:, :72]
            nulls["random_state"].append(summary(mi_one_mode(Vo @ Vo.conj().T)))
        print(k, f"{time.time()-t0:.0f}s", flush=True)
    rates = {}
    for key, L in nulls.items():
        holes = np.mean([(r["betti"][1] + r["betti"][2]) >= 1 for r in L])
        lowd = np.mean([0.8 <= r["D"] <= 3.5 for r in L])
        multi = np.mean([r["pieces"] > 1 for r in L])
        rates[key] = dict(n=len(L), fp_holes=float(holes), fp_lowdim=float(lowd), multi_piece=float(multi),
                          mean_b1=float(np.mean([r["betti"][1] for r in L])),
                          mean_b2=float(np.mean([r["betti"][2] for r in L])),
                          mean_D=float(np.nanmean([r["D"] for r in L])),
                          max_p1_rel=float(np.max([r["p1_rel"] for r in L])),
                          max_p1_rel_2nd=float(np.max([r["p1_rel_2nd"] for r in L])))
        print(key, rates[key], flush=True)
    res["nulls"] = nulls
    res["rates"] = rates
    res["verdict"] = {"P2_pass(fp_holes<=0.1 all)": bool(all(r["fp_holes"] <= 0.1 for r in rates.values()))}
    res["time_s"] = time.time() - t0
    print(res["verdict"], res["time_s"])
    with open(os.path.join(OUT, "pre_p2.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
