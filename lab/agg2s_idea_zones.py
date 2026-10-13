"""agg2s_idea_zones.py - IDEA TEST (post hoc, declared): how much of the nucleus cusp tail does a coarse level lose with
Brillouin zones of different SHAPE but the same number of cells (same zone volume)? The lost contact and <p^4> scale as
I = (2pi)^-3 int dOmega / R(Omega) and the cell energy as J = (2pi)^-3 int dOmega / (3 R^3), R(Omega) = distance from the
centre to the zone surface. Zones: cube (simple cubic lattice), truncated octahedron (body-centred cubic lattice),
rhombic dodecahedron (face-centred cubic lattice), sphere (ideal). Writes lab/out/agg2s_idea_zones.json"""
import json, os, itertools
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
rng = np.random.default_rng(7)
n = rng.normal(size=(4_000_000, 3))
n /= np.linalg.norm(n, axis=1)[:, None]


def polytope(G):
    """zone {k : k.g <= 1 for g in G}; R(n) = 1 / max_g(n.g)"""
    G = np.asarray(G, float)
    return 1.0 / np.max(n @ G.T, axis=1)


def volume_from_R(R):
    return 4 * np.pi * np.mean(R ** 3) / 3


cube = polytope([s * e for e in np.eye(3) for s in (1, -1)])                        # |k_i| <= 1
to_g = [s * e for e in np.eye(3) for s in (1, -1)] + [np.array(v) / 1.5 for v in itertools.product((1, -1), repeat=3)]
trunc_oct = polytope(to_g)                                                          # |k_i| <= 1, |kx|+|ky|+|kz| <= 1.5
rd_g = [np.array(v, float) for v in set(itertools.permutations((1, 1, 0))) | set(itertools.permutations((1, -1, 0))) |
        set(itertools.permutations((-1, -1, 0)))]
rhomb = polytope(rd_g)                                                              # |k_i| + |k_j| <= 1
out = {}
Vc = volume_from_R(cube)
for name, R in (("cube (simple cubic)", cube), ("truncated octahedron (bcc)", trunc_oct), ("rhombic dodecahedron (fcc)", rhomb)):
    s = (Vc / volume_from_R(R)) ** (1 / 3)                                          # rescale to the cube's volume
    Rs = R * s
    out[name] = dict(I_rel=float(np.mean(1 / Rs) / np.mean(1 / cube)), J_rel=float(np.mean(1 / Rs ** 3) / np.mean(1 / cube ** 3)),
                     R_min_over_R_max=float(Rs.min() / Rs.max()))
Rsph = (3 * Vc / (4 * np.pi)) ** (1 / 3)
out["sphere"] = dict(I_rel=float((1 / Rsph) / np.mean(1 / cube)), J_rel=float((1 / Rsph ** 3) / np.mean(1 / cube ** 3)), R_min_over_R_max=1.0)
json.dump(out, open(os.path.join(OUT, "agg2s_idea_zones.json"), "w"), indent=1)
for k, v in out.items():
    print(f"{k:30s} lost contact/<p^4> x {v['I_rel']:.4f}   cell energy x {v['J_rel']:.4f}")
