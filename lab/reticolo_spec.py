"""reticolo_spec.py - R5: export the lattice format of the sparse hierarchical light lattice (v2 engine,
lab/reticolo_core.py) as assets/atomo/reticolo_spec.json, for the Rust lab to implement later."""
import os, json, glob
import numpy as np
import reticolo_core as rc

HERE = os.path.dirname(os.path.abspath(__file__))
DST = os.path.join(HERE, "..", "assets", "atomo", "reticolo_spec.json")


def r2_summary():
    rows = []
    for p in sorted(glob.glob(os.path.join(rc.OUT, "reticolo_r2_cases", "*.json"))):
        r = json.load(open(p))
        rows.append({k: r[k] for k in ("d", "N", "ang", "direction", "T", "R", "lost", "fidelity", "phase",
                                        "hel_minus", "hel_static")})
    return rows


spec = {
    "name": "reticolo_universo",
    "version": 2,
    "engine_reference": "lab/reticolo_core.py (Python + numba), doc docs/RETICOLO_UNIVERSO.md",
    "units": "h0 = 1 (finest cell), c = 1; field F = E + i c B (Riemann-Silberstein), energy = sum |F|^2 h^d",
    "levels": {
        "index": "k = 0, 1, 2, ... (u8 is enough for the known universe: see docs)",
        "cell_size": "h_k = 2^k h0",
        "tick": "dt_k = h_k / c",
        "rule": "identical on every level for light (massless spin-1 zig-zag walker, src/zigzag/walk.rs Light3D::tick)",
        "rule_massive_hook": "Dirac walker with theta_k = m c h_k / hbar is valid while theta_k << 1; on coarser levels a "
                             "massive object must be an agglomerate with a learned level rule (lab/agg_zz_*): NOT part of v2",
    },
    "blocks": {
        "size": {"2d": 16, "3d": 8},
        "cell_order": "local index s = x + B*y (+ B*B*z), x fastest",
        "payload": "B^d cells x 3 complex (f64 in the prototype; f32 would halve memory)",
        "key": "integer block coordinates per level (prototype: 21 bits per axis = +-2^20 blocks; Rust: i64 or i128 "
               "per axis, see docs for the address range)",
        "missing_block": "exactly zero (one shared zero sentinel block)",
        "children": "block b of level k+1 covers blocks 2b + {0,1}^d of level k (octree / quadtree)",
    },
    "tick_of_one_level": {
        "axes_order": "even own tick: x, y, (z); odd own tick: reversed (palindromic pair)",
        "one_axis_j": "j1 = (j+1)%3, j2 = (j+2)%3; c+ = (F_j1 - i F_j2)/sqrt2 moves +1 cell along j; "
                      "c- = (F_j1 + i F_j2)/sqrt2 moves -1 cell; F_j stays; then F_j1 = (c+ + c-)/sqrt2, "
                      "F_j2 = i (c+ - c-)/sqrt2. Ghost values from the neighbour block (zero if missing).",
    },
    "time_schedule": {
        "kind": "Berger-Oliger, ratio 2 in space and time",
        "pseudocode": [
            "for each finest tick t = 0,1,2,...:",
            "  for i = top..0: if t % 2^i == 0: tick(level i)        # coarse levels first, they wait ahead",
            "  t += 1",
            "  for i = 1..top: if t % 2^i != 0: break; if (t >> i) % SYNC_EVERY != 0: break;",
            "                  sync(level i-1 -> level i); maybe regrid(level i-1)",
            "  if t % 2^top == 0: maybe regrid(top)"],
        "SYNC_EVERY": rc.SYNC_EVERY,
        "why_SYNC_EVERY_2": "a single tick (x then y) is not palindromic: syncing after every coarse tick couples the "
                            "split-step error O(kx ky h^2) of the two levels and creates/destroys ~1% of energy at 30 deg "
                            "(v1); syncing after every PAIR removes it (v2)",
    },
    "interface": {
        "band_width_fine_cells": rc.W_BAND,
        "band": "fine cells within W of any missing neighbour block (all 3^d directions)",
        "prolongation_into_band": {
            "kind": "tensor-product cubic Lagrange of coarse cell centres",
            "fine_index_even_2m": {"coarse_nodes": "m-2, m-1, m, m+1", "weights": rc.WP_EVEN.tolist()},
            "fine_index_odd_2m+1": {"coarse_nodes": "m-1, m, m+1, m+2", "weights": rc.WP_ODD.tolist()},
        },
        "restriction": {
            "kind": "4th-order midpoint interpolation, tensor product",
            "coarse_gc_from_fine": "2gc-1, 2gc, 2gc+1, 2gc+2",
            "weights": rc.WR.tolist(),
            "where": "only coarse cells whose whole 4^d stencil is fine INTERIOR (not band, not missing)",
        },
        "sync_order": "restriction first, then prolongation into the band",
        "energy_once": "count fine interior cells; count a coarse cell unless all its 2^d children are fine interior",
        "known_limit": "for content below ~4 cells per wavelength on the coarse level the interpolating restriction "
                       "ALIASES instead of averaging (a ghost of the gamma survives up to ~8 levels above it at a few % "
                       "amplitude); a conservative low-pass restriction for unresolved content is debt",
    },
    "regrid": {
        "every_own_ticks": 8,
        "flag": "level k < top: max over block of |F_k - P(F_k+1)| > tol_detail (evaluated on every 2nd cell per axis); "
                "top level: max |F| > tol_top",
        "tol_detail": 1e-5, "tol_top": 1e-4,
        "halo_blocks": 1,
        "parent_margin_coarse_blocks": 2,
        "nesting_check": "a fine block is kept only if its parent block and the 3^d neighbours of the parent exist",
        "new_blocks": "filled by prolongation from the coarser level; removed blocks: the coarser level already holds "
                      "the restricted field",
    },
    "parameters": {
        "physics": 0,
        "engineering": ["B", "W_band", "SYNC_EVERY", "interpolation orders (cubic/cubic)", "tol_detail", "tol_top",
                        "halo", "parent_margin", "regrid interval"],
        "learned": 0,
    },
    "validation_v2": r2_summary(),
}
os.makedirs(os.path.dirname(DST), exist_ok=True)
json.dump(spec, open(DST, "w"), indent=1)
print("written", os.path.normpath(DST), len(spec["validation_v2"]), "R2 cases")
