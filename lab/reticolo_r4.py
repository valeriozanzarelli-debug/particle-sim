"""reticolo_r4.py - R4: the cost of holding things in the sparse hierarchical lattice.
(a) the same pulse (lambda = 24 cells of its level, sigma = 36) at home level k = 0..40: active cells must be
    identical (the absolute scale costs nothing);
(b) a 'gamma' at level 0 + a 'radio' at level m (m = 2, 5, 10 built with the real engine; 20 from R3): the
    levels in between are a chimney whose cost must grow linearly with m; a block-only model of the chimney
    (nesting rule alone) is checked against the engine and extended to m = 40, 64, 128;
(c) honest estimates for a Milky-Way-size box and what does NOT fit.
Writes lab/out/reticolo_r4.json."""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reticolo_r3 as r3
from reticolo_core import packet, single_level_world, OUT

res = {}
t0 = time.time()
# ---------------------------------------------------------------- (a) one pulse, any scale
F = packet((512, 512), 2, 24, 60.0, [256, 256], 36.0, amp=1.0)
rows = []
for k in (0, 1, 5, 10, 20, 30, 40):
    w = single_level_world(k, 2, 16, np.array([-256, -256]), F, tol=1e-4)
    L = w.levels[0]
    rows.append(dict(k=k, h_over_h0=2.0 ** k, blocks=int(L.nb), cells=int(L.nb * L.Bd),
                     region_side_cells=int((np.ptp(L.coords[:, 0]) + 1) * 16), region_side_h0=float((np.ptp(L.coords[:, 0]) + 1) * 16 * L.h), bytes=int(L.F.nbytes)))
    print(rows[-1], flush=True)
res["a_one_pulse"] = rows
res["a_identical_cells"] = len(set(r["cells"] for r in rows)) == 1

# ---------------------------------------------------------------- (b) chimney
def chimney_model(level0_coords, m, margin=r3.PAR_MARGIN):
    c = level0_coords
    out = [len(c)]
    for k in range(1, m):
        c = r3.dil(np.unique(np.floor_divide(c, 2), axis=0), margin)
        out.append(len(c))
    return out


real = {}
lvl0 = None
for m in (2, 5, 10):
    t1 = time.time()
    w, _, _ = r3.build(m, (-40, -30))
    tab = r3.level_table(w)
    real[m] = dict(blocks=[r["blocks"] for r in tab], cells=w.cells(), bytes=sum(r["bytes"] for r in tab),
                   build_s=time.time() - t1)
    if lvl0 is None:
        lvl0 = w.levels[0].coords.copy()
        radio_blocks = w.levels[-1].nb
    print(m, real[m], flush=True)
res["b_engine"] = real
try:
    lad = json.load(open(os.path.join(OUT, "reticolo_r3_ladder.json")))
    res["b_engine_m20_from_R3"] = dict(blocks=[r["blocks"] for r in lad["levels_t0"]], cells=lad["cells_t0"])
except Exception:
    pass
model = {}
for m in (2, 5, 10, 20, 40, 64, 128):
    ch = chimney_model(lvl0, m)
    model[m] = dict(blocks_levels_0_to_m_minus_1=ch, chimney_blocks=int(sum(ch[1:])),
                    total_blocks_without_radio=int(sum(ch)))
res["b_model"] = model
# per-level chimney cost (fixed point of the nesting rule) and linear fit on the engine numbers
ms = np.array(sorted(real))
tot = np.array([real[m]["cells"] for m in ms])
slope = np.polyfit(ms, tot, 1)
res["b_linear_fit_cells_vs_m"] = dict(slope_cells_per_level=float(slope[0]), intercept=float(slope[1]))
res["b_chimney_fixed_point_blocks_per_level"] = int(chimney_model(lvl0, 40)[-1])

# ---------------------------------------------------------------- (c) estimates (3-D, blocks 8^3, F = 3 complex)
LY = 9.4607e15
MW = 1.0e5 * LY                     # Milky Way disc diameter ~ 1e5 ly
UNIV = 8.8e26                       # observable universe diameter (m)
BYTES_CELL = 48                     # 3 x complex128 (24 with complex64)
side = rows[0]["region_side_cells"]  # 2-D measured side of the active region for the pulse, in its cells
cells3d_pulse = side ** 3
lam_radio = 1.0                     # m
lam_gamma = 1.23984e-12             # m (1 MeV)
h_radio = lam_radio / 24
h_gamma = lam_gamma / 24
lv_gr = np.log2(lam_radio / lam_gamma)
blocks_chim = res["b_chimney_fixed_point_blocks_per_level"]
est = {
    "pulse_cells_3d": cells3d_pulse,
    "pulse_bytes_3d": cells3d_pulse * BYTES_CELL,
    "note_box": "a single pulse needs NO level above its own: the box size L does not enter the cost at all",
    "levels_gamma_to_radio": float(lv_gr),
    "chimney_cells_3d_per_level": float((np.sqrt(blocks_chim) * 16) ** 3),
    "gamma_plus_radio_cells_3d": float(2 * cells3d_pulse + lv_gr * (np.sqrt(blocks_chim) * 16) ** 3),
    "address_bits_level0_gamma_across_universe": float(np.log2(UNIV / h_gamma)),
    "address_bits_level0_gamma_across_milky_way": float(np.log2(MW / h_gamma)),
    "radio_filling_milky_way_cells": float((MW / h_radio) ** 3),
    "radio_filling_milky_way_bytes": float((MW / h_radio) ** 3 * BYTES_CELL),
    "gamma_filling_milky_way_cells": float((MW / h_gamma) ** 3),
    "ticks_radio_cross_milky_way": float(MW / h_radio),
    "ticks_gamma_cross_milky_way": float(MW / h_gamma),
    "gamma_ticks_per_radio_tick": float(2 ** round(lv_gr)),
    "radio_min_beam_width_cross_mw_m": float(np.sqrt(MW * lam_radio / np.pi)),
    "gamma_min_beam_width_cross_mw_m": float(np.sqrt(MW * lam_gamma / np.pi)),
}
est["radio_min_beam_cells_across"] = est["radio_min_beam_width_cross_mw_m"] / h_radio
est["gamma_min_beam_cells_across"] = est["gamma_min_beam_width_cross_mw_m"] / h_gamma
res["c_estimates"] = est
res["wall_s"] = time.time() - t0
json.dump(res, open(os.path.join(OUT, "reticolo_r4.json"), "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in res.items() if k in ("a_identical_cells", "b_linear_fit_cells_vs_m",
                                                        "b_chimney_fixed_point_blocks_per_level", "c_estimates")},
                 indent=1, default=float))
