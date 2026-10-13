"""locale_cands.py - the candidates of the pre-registration (lab/out/locale_prereg.json), built with locale_core."""
import itertools
import numpy as np
from locale_core import (make_walk, composed_walk, zigzag_ticks, CUBE_AXES, LATTICES, S4_SEQ)

H_WILSON = 0.1
AX = CUBE_AXES
HEXD = [np.array([1, 0, 0.]), np.array([0, 1, 0.]), np.array([-1, 1, 0.])]
FCC = [np.array(v, float) for v in [(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]]
BCC = [np.array(v, float) for v in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]


def wil(dim, h=H_WILSON, axis="Y"):
    return dict(h=h, vecs=[2 * e for e in AX[:dim]], axis=axis)


def perms(ds):
    return [[ds[i] for i in p] for p in itertools.permutations(range(len(ds)))]


def fermion(name, dim=3, theta=0.1, h=H_WILSON):
    if name in ("F0", "F1"):
        return make_walk("dirac", zigzag_ticks(dim), theta=theta, dim=dim, name=name)
    if name == "F2":
        return make_walk("dirac", zigzag_ticks(dim), theta=theta, dim=dim, wilson=wil(dim, h, "Y"), name=name)
    if name == "F2b":
        return make_walk("dirac", zigzag_ticks(dim), theta=theta, dim=dim, wilson=wil(dim, h, "X"), name=name)
    if name == "F3":
        return composed_walk("dirac", S4_SEQ, dim=dim, theta=theta, name=name)
    raise KeyError(name)


# light: (walk, sites per unit volume (lattice length 1), dim)
def light(name):
    if name == "L0_3D":
        return make_walk("light", zigzag_ticks(3), name=name), 1.0, 3
    if name == "L0_2D":
        return make_walk("light", zigzag_ticks(2), dim=2, name=name), 1.0, 2
    if name == "L1_2D":
        return make_walk("light", perms(HEXD), basis=LATTICES["hex"][0], dim=2, name=name), 2 / np.sqrt(3), 2
    if name == "L2_3D":
        return composed_walk("light", S4_SEQ, dim=3, name=name), 1.0, 3
    if name == "L2_2D":
        return composed_walk("light", S4_SEQ, dim=2, name=name), 1.0, 2
    if name == "L3_3D":
        return make_walk("light", perms(AX), name=name), 1.0, 3
    if name == "L4a_3D":
        return make_walk("light", [FCC, FCC[::-1]], name=name), 0.5, 3
    if name == "L4b_3D":
        return make_walk("light", [BCC, BCC[::-1]], name=name), 0.25, 3
    if name == "L5_3D":
        return make_walk("light", [AX, AX[::-1]] * 6 + [BCC, BCC[::-1]], name=name), 1.0, 3
    raise KeyError(name)


LIGHTS = ["L0_3D", "L0_2D", "L1_2D", "L2_3D", "L2_2D", "L3_3D", "L4a_3D", "L4b_3D", "L5_3D"]
FERMIONS = ["F0", "F1", "F2", "F2b", "F3"]
