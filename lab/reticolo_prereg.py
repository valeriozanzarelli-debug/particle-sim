"""reticolo_prereg.py - pre-registration of the SPARSE HIERARCHICAL LATTICE study (R1-R5), written BEFORE the
validation runs.  Engine: lab/reticolo_core.py (committed together with this file; only a kernel smoke test
- one level, periodic box, 10 ticks against the Bloch matrix of zz_core: max error 3e-14 in 2-D and 3-D - was run).
Writes lab/out/reticolo_prereg.json."""
import json, os, subprocess, datetime
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)
head = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()


def vg_err(N, ang):          # plane-wave group-velocity error of the rule in the x-y plane, cos w = cos kx cos ky
    a = np.radians(ang)
    om = lambda k: np.arccos(np.cos(k * np.cos(a)) * np.cos(k * np.sin(a)))
    k, dk = 2 * np.pi / N, 1e-6
    return (om(k + dk) - om(k - dk)) / (2 * dk) - 1


P = {
 "meta": {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": head,
          "hardware": "i5 4 cores shared with 3 other agents; numba 1 thread; each script <= 10 min"},
 "question": "Valerio: 'our lattice must potentially hold the whole Milky Way or the known universe; to see a radio "
             "wave whole you zoom far out, but we want to see it'. Build a SPARSE HIERARCHICAL lattice (levels "
             "h_k = 2^k h0, only blocks where something is non-zero), the same zig-zag light rule on every level, "
             "levels joined by a restriction/prolongation band (as ZOOM_NODO).",
 "design_frozen": {
   "rule": "massless spin-1 zig-zag walker (F = E + icB), identical on every level (it has no length scale)",
   "blocks": "2-D: 16x16 cells; 3-D: 8x8x8 cells; zero sentinel outside",
   "time": "Berger-Oliger: one level-k tick = two level-(k-1) ticks; each level keeps its own parity",
   "restriction": "coarse cell <- 4th-order midpoint interpolation (-1,9,9,-1)/16 per axis, only from fine INTERIOR",
   "prolongation": "fine band (outer W=4 fine cells) <- tensor-cubic Lagrange of the coarse level at every sync; "
                   "band evolves with the fine rule between syncs (no time interpolation)",
   "energy": "composite: fine interior + coarse cells not made only of fine-interior children (each point once)",
   "regrid": "level k block kept if max|F_k - P(F_k+1)| > tol (coarsest: max|F| > tol), dilated by 1 block, plus "
             "parents of the finer level's blocks dilated by 1 block; new blocks by prolongation",
   "free_parameters": "B, W=4, interpolation orders (cubic/cubic), tol, halo = engineering choices; 0 physics parameters "
                      "(c and the rule are those of docs/ZIGZAG.md)"},
 "R2a_speed": {
   "test": "packet lambda = 24 cells, sigma = 4 lambda, on one sparse level with regrid, levels k = 0,1,5,10,20; "
           "centroid velocity over 480 own ticks, in physical units (cells*h_k)/(ticks*dt_k)",
   "criterion": "|v/c - 1| <= 0.5 % on every level",
   "prediction": {"axis": "PASS on every level (axis shift is exact); arrays bitwise identical across levels",
                  "diagonal_45deg_24cells": "FAIL: plane-wave group velocity of the rule %.2f %% (needs >= 32 cells/lambda: %.2f %%)"
                  % (100 * vg_err(24, 45), 100 * vg_err(32, 45)),
                  "30deg_24cells": "FAIL: %.2f %%" % (100 * vg_err(24, 30))},
   "also": "centroid velocity measured on each side of a level interface in R2b (axis): |v-1| <= 0.5 %"},
 "R2b_interface": {
   "test": "2-D strip periodic in y; fine level 0 next to coarse level 1; Gaussian packet sigma = 3 lambda (axis) of "
           "helicity +1 crossing fine->coarse and coarse->fine; N = coarse cells per lambda in {2,3,4,6,8,12,16,24,32}; "
           "normal incidence and 30 deg; T = energy beyond the interface (+16 fine cells), R = energy behind it, "
           "lost = 1 - T - R (composite energy)",
   "criterion": "for N >= 16: T >= 99 %, R <= 1 %, both directions, normal incidence (headline); 30 deg reported with the same thresholds",
   "prediction": "PASS at N >= 16 with R < 1e-4 (cubic interpolation error ~ (k h_c)^4); short waves N <= 4 MUST fail "
                 "(T < 90 %): a coarse lattice cannot carry them, physics of every multigrid; reported as a curve"},
 "R2c_phase_helicity": {
   "test": "transmitted packet compared on the coarse grid with the same packet evolved on the uniform lattice of its "
           "origin level (no interface); helicity fractions from the lattice eigen-projectors of the two-tick walk",
   "criterion": "normal incidence, N >= 16: |phase error| <= 0.01 rad, fidelity |<ref|out>|/(|ref||out|) >= 0.999, "
                "wrong helicity + static mode <= 1e-3 of the transmitted energy; 30 deg: fidelity and helicity only "
                "(the phase there also contains the different dispersion of the two levels)",
   "prediction": "PASS"},
 "R2_3D": {"test": "same as R2b in 3-D (blocks 8^3, strip periodic in y and z), N = 16 and 8, normal incidence",
           "criterion": "N = 16: T >= 99 %, R <= 1 %", "prediction": "PASS", "status": "if budget allows"},
 "R3_demo": {
   "world": "ONE world, levels 0..20: gamma pulse lambda = 24 h0 at level 0 (sigma = 18 h0, moving +x), radio pulse "
            "lambda = 24 h20 = 24*2^20 h0 at level 20 (sigma = 36 h20, moving at 60 deg), the gamma inside the radio "
            "envelope; levels 1..19 exist only as a nested 'chimney' around the gamma",
   "criteria": {"cells": "total active cells <= 1e6 and field memory <= 100 MB at all times",
                "radio": "after n_r >= 1 level-20 ticks (target 2): radio on level 20 vs the same radio evolved alone: "
                         "relative L2 difference <= 1e-3",
                "gamma": "after 2^20 n_r level-0 ticks: centroid moved 2^20 n_r h0 +- 0.5 h0, fidelity vs the exactly "
                         "shifted initial packet >= 0.999, energy conserved within 1e-3"},
   "prediction": "PASS; wall time in Python+numba ~10-20 min per radio tick (2^20 gamma ticks): TIME, not memory, is the wall",
   "figure": "lab/out/reticolo_zoom_*.png, zoom-out ladder from level 0 to level 20"},
 "R4_scaling": {"claim": "active cells to hold a pulse of wavelength lambda in a box L = (cells at its home level, "
                         "independent of the absolute scale) + (number of levels) x (chimney blocks)",
                "test": "same packet at home level k = 0..40 (cells must be identical), plus worlds with two pulses "
                        "2^m apart in wavelength, m = 0..40 (overhead must grow linearly with m)",
                "prediction": "PASS (exactly identical cells; linear overhead)"},
 "R5": "assets/atomo/reticolo_spec.json + docs/RETICOLO_UNIVERSO.md",
}
json.dump(P, open(os.path.join(OUT, "reticolo_prereg.json"), "w"), indent=1)
print(json.dumps(P["R2a_speed"]["prediction"], indent=1))
