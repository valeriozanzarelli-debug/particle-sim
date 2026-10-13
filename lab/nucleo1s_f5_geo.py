"""nucleo1s_f5_geo.py - phase 5: GEOMETRIC predictions of the nucleus factors (s0 on the 2 nearest sites of a colour,
s1 on the 6 next ones), with no learning:
 P_A  colour-cell average: each site gets the average of -Z alpha/r over its own Voronoi cell in its colour (a bcc lattice,
      cell = truncated octahedron of volume 4 a^3). s_i = r_i <1/r>_cell. Independent of Z and a.
 P_B1 exact scattering of the first cell (radial form of the same rule, uniform cells a, first site at a/2): the factor
      s on the Coulomb phase of the first site that makes the radial 1s equal to the 1s of the same rule with the
      inner zone resolved (zoom, the rule's own limit). Depends on Z and a.
 P_B2 as P_B1 but the first radial cell is [0, sqrt3 a] (its site at r0 = sqrt3/2 a, the radius of the 3-D nearest
      sites), cells a beyond.
Units hbar = c = m = 1."""
import os, sys, itertools
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA  # noqa: E402
from nucleo1s_core import Mesh, GradedWalk, rot_chi, mesh_zoom, dirac_binding  # noqa: E402


def colour_cell_avg(c, N=240, half=2.0):
    g = (np.arange(N) + 0.5) / N * 2 * half - half
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    P = np.stack([X, Y, Z], -1).reshape(-1, 3) + c
    nb = [np.array(v) for v in itertools.product(range(-2, 3), repeat=3)
          if (all(x % 2 == 0 for x in v) or all(x % 2 for x in v)) and any(v)]
    d0 = np.sum((P - c)**2, 1)
    keep = np.ones(len(P), bool)
    for v in nb:
        keep &= d0 <= np.sum((P - c - v)**2, 1)
    r = np.linalg.norm(P[keep], axis=1)
    return float(np.linalg.norm(c) * np.mean(1 / r)), float(keep.sum() * (2 * half / N)**3)


def P_A():
    s0, v0 = colour_cell_avg(np.array([0.5, 0.5, 0.5]))
    s1, v1 = colour_cell_avg(np.array([1.5, 0.5, 0.5]))
    return dict(s0=s0, s1=s1, vol0=v0, vol1=v1)


def _binding(edges, Za, s, guess):
    m = Mesh(edges)
    M = rot_chi(m.h, -1, m.r, Za)
    M[0] = rot_chi(m.h[:1], -1, m.r[:1], Za * s)[0]
    gw = GradedWalk(m, -1, Za, M=M)
    w, x, info = gw.eigen(1.0 - guess)
    return 1.0 - w


def target_binding(a, Za, rmax):
    m = mesh_zoom(a, rmax, 12, 16)
    return _binding(m.e, Za, 1.0, dirac_binding(1, -1, Za))


def P_B(Z, a, variant=1):
    Za = Z * ALPHA
    rmax = max(40.0, 12 * 137.036 / Z)
    tgt = target_binding(a, Za, rmax)
    if variant == 1:
        edges = np.arange(int(np.ceil(rmax / a)) + 1) * a
    else:
        R0 = np.sqrt(3) * a
        edges = np.concatenate([[0.0], R0 + np.arange(int(np.ceil((rmax - R0) / a)) + 1) * a])
    f = lambda s: _binding(edges, Za, s, tgt) - tgt
    grid = np.linspace(0.05, 2.0, 40)
    vals = [f(s) for s in grid]
    roots = [brentq(f, grid[i], grid[i + 1], xtol=1e-10) for i in range(len(grid) - 1) if vals[i] * vals[i + 1] < 0]
    return dict(Z=Z, a=a, target=tgt, roots=roots, s=(min(roots, key=lambda x: abs(x - 1)) if roots else None),
                bare=_binding(edges, Za, 1.0, tgt))


if __name__ == "__main__":
    print(P_A())
    print(P_B(60, 0.5, 1))
