"""zoomknot_learn.py - learn the agglomerate couplings k2_l, k4_l of the Z2 hierarchy (pre-registered recipe) and run
the held-out check. For level l >= 1 (spacing h, shell = box R_l minus box R_{l-1}):
  fine truth   : bare lattice energy (e2, e4 separately) of the hedgehog f(r/lam) at spacing h/4 on the shell;
  agglomerate  : the same configuration blocked 4x4x4 -> 1 (two vertex full-weighting steps, 1D weights
                 [1,2,3,4,3,2,1]/16 at spacing h/4), renormalised, bare lattice energy at spacing h on the shell;
  learn        : k2 = sum E2f E2c / sum E2c^2 over lam in {0.7, 0.85, 1.2, 1.4}; same for k4.
  check (lam=1): per-level bare / learned energy vs continuum integral of the exact T00 on the shell (trapezoid h/8).
Octant grids (the hedgehog energy density is even in x, y, z). Writes lab/out/zoomknot_learn.json."""
import sys, os, json, math, time, numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot_core as kc
import zoomknot_core as zc

torch.set_num_threads(1)
DT = torch.float64
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
PRE = json.load(open(os.path.join(OUT, "zoomknot_prereg.json")))
LV = PRE["configs"]["Z2"]
prof = kc.radial_bvp(60.0)
O = np.load(os.path.join(OUT, "press_profiles.npz"))
LAMS = [0.7, 0.85, 1.2, 1.4]


def oct_axis(h, Rout):
    n = int(round(Rout / h)) + 2
    return torch.arange(-2, n + 1, dtype=DT) * h


def oct_weight(x, Rin, Rout):
    O1 = torch.where(x > 1e-9, torch.ones_like(x), torch.where(x > -1e-9, 0.5 * torch.ones_like(x), 0 * x))
    To = zc.trap1(x, Rout) * O1
    Ti = zc.trap1(x, Rin) * O1 if Rin is not None else 0 * x
    W = To[:, None, None] * To[None, :, None] * To[None, None, :]
    if Rin is not None:
        W = W - Ti[:, None, None] * Ti[None, :, None] * Ti[None, None, :]
    return 8 * W


def hh(x, lam, shift=(0.0, 0.0, 0.0)):
    X, Y, Z = torch.meshgrid(x + shift[0], x + shift[1], x + shift[2], indexing="ij")
    return kc.hedgehog(X, Y, Z, {"r": prof["r"] * lam, "f": prof["f"]})


def lat_E(phi, h, W):
    e2, e4, _ = zc.densities(phi, h)
    return float((W * e2).sum()) * h**3, float((W * e4).sum()) * h**3


def blocked(xc, hf, lam):
    w1 = np.array([1, 2, 3, 4, 3, 2, 1]) / 16.0
    acc = None
    for a in range(7):
        for b in range(7):
            for c in range(7):
                w = w1[a] * w1[b] * w1[c]
                t = hh(xc, lam, ((a - 3) * hf, (b - 3) * hf, (c - 3) * hf)) * w
                acc = t if acc is None else acc + t
    return kc.normalize(acc)


@torch.no_grad()
def main():
    t0 = time.time()
    res = {"levels": LV, "lams_train": LAMS, "per_level": []}
    k2s, k4s = [1.0], [1.0]
    for l in range(1, len(LV)):
        h = LV[l]["h"]; Rin = LV[l - 1]["R"]
        Rout = LV[l]["R"] if LV[l]["R"] is not None else LV[l]["N"] * h
        hf = h / 4
        xf = oct_axis(hf, Rout); Wf = oct_weight(xf, Rin, Rout)
        xc = oct_axis(h, Rout); Wc = oct_weight(xc, Rin, Rout)
        rows = []
        for lam in LAMS + [1.0]:
            E2f, E4f = lat_E(hh(xf, lam), hf, Wf)
            E2c, E4c = lat_E(blocked(xc, hf, lam), h, Wc)
            E2s, E4s = lat_E(hh(xc, lam), h, Wc)            # sampled (not blocked), for information
            rows.append(dict(lam=lam, E2f=E2f, E4f=E4f, E2c=E2c, E4c=E4c, E2s=E2s, E4s=E4s))
        tr = [r for r in rows if r["lam"] != 1.0]
        k2 = sum(r["E2f"] * r["E2c"] for r in tr) / sum(r["E2c"] ** 2 for r in tr)
        k4 = sum(r["E4f"] * r["E4c"] for r in tr) / sum(r["E4c"] ** 2 for r in tr)
        # continuum check at lam = 1: trapezoid of the exact T00 at h/8
        x8 = oct_axis(h / 8, Rout); W8 = oct_weight(x8, Rin, Rout)
        X, Y, Z = torch.meshgrid(x8, x8, x8, indexing="ij")
        rr = torch.sqrt(X**2 + Y**2 + Z**2).numpy()
        Econt = float((W8.numpy() * np.interp(rr, O["r"], O["T00"])).sum()) * (h / 8) ** 3
        te = rows[-1]
        Eb = te["E2c"] + te["E4c"]; El = k2 * te["E2c"] + k4 * te["E4c"]; Ef = te["E2f"] + te["E4f"]
        Es = te["E2s"] + te["E4s"]
        d = dict(level=l, h=h, block_in_finest_balls=h / LV[0]["h"], shell=[Rin, Rout], k2=k2, k4=k4,
                 E_cont=Econt, E_bare_blocked=Eb, E_learned_blocked=El, E_fine_h4=Ef, E_bare_sampled=Es,
                 err_bare=Eb / Econt - 1, err_learned=El / Econt - 1, err_fine_h4=Ef / Econt - 1,
                 err_bare_sampled=Es / Econt - 1, train_rows=rows, E4_share=te["E4f"] / Ef)
        res["per_level"].append(d)
        k2s.append(k2); k4s.append(k4)
        print(f"level {l} h {h:5.2f} shell {Rin}-{Rout}: k2 {k2:.5f} k4 {k4:.5f} | Econt {Econt:.5f} "
              f"err bare {d['err_bare']:+.4%} learned {d['err_learned']:+.4%} fine(h/4) {d['err_fine_h4']:+.4%} "
              f"sampled {d['err_bare_sampled']:+.4%}  t {time.time()-t0:.0f}s", flush=True)
    res["k2"], res["k4"] = k2s, k4s
    Eex = 1.231446 * 12 * math.pi**2
    pl = res["per_level"]
    res["sum_static_err_learned_over_E"] = sum(d["E_learned_blocked"] - d["E_cont"] for d in pl) / Eex
    res["sum_static_err_bare_over_E"] = sum(d["E_bare_blocked"] - d["E_cont"] for d in pl) / Eex
    res["n_levels_learned_better"] = int(sum(abs(d["err_learned"]) < abs(d["err_bare"]) for d in pl))
    res["Z2d_i_pass"] = bool(res["n_levels_learned_better"] >= 6)
    res["Z2e_pass"] = bool(abs(res["sum_static_err_learned_over_E"]) < 3e-3)
    res["n_learned_numbers"] = 2 * (len(LV) - 1)
    res["seconds"] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, "zoomknot_learn.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ("per_level", "levels")}, indent=1))


main()
