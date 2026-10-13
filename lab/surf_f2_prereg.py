"""surf_f2_prereg.py -> lab/out/surf_f2_prereg.json : PHASE 2 of SURFISTI_ESATTI, written BEFORE the production runs.
Development done before (declared): the Bloch analysis of the valley spin frames (U2(q+d) = V U2(R d) V^dag, found exactly
at 5e-16: q=(pi/2,pi/2,0): V = sigma_x with L<->R swap, R = (-,-,+); (pi/2,0,pi/2): V = sigma_x; (0,pi/2,pi/2): V = sigma_y,
R = (-,+,+)); reading docs/BATTIMENTO_1S.md (master 83485bc). Nothing of V or S below has been measured yet.
"""
import json
import os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
# continuum reference (docs/EQUILIBRIO_QUANTISTICO.md, |psi|^4 start, 16x16 cells): tau(M=16) = 4.3, tau(M=25) = 2.4;
# rms spread of the mode frequencies E = (n^2+m^2)/2: 4.016 (M=16), 6.116 (M=25)
cont = {16: dict(tau=4.3, sigma=4.0156), 25: dict(tau=2.4, sigma=6.1156)}
for M in cont:
    cont[M]["tau_sigma"] = cont[M]["tau"] * cont[M]["sigma"]
P = dict(
    study="SURFISTI_ESATTI fase 2", date="2026-10-13",
    part1=dict(
        title="Valentini on the lattice rule with the exact flow (no fake mixing)",
        setup=dict(walk="pre_core.LatticeWalk (2-D zig-zag, mirror walls), fixed axis order x,y -> period 1 tick",
                   lattices=[dict(n=48, theta=0.25, cell=4, ticks=4000), dict(n=96, theta=0.125, cell=8, ticks=8000)],
                   states="psi = sum of the M eigenmodes of the one-tick U closest to exp(-i omega0), random phases "
                          "(pre_core.modes_state, fixed ARPACK start vector, 1 thread), omega0 = 2 theta",
                   M=[1, 4, 16, 25], starts=["born (|psi|^2)", "rho2 (|psi|^4)", "uniform", "half box (x < n/2)"],
                   surfers=50000, mover="S4 exact flow (net link flux of the tick, Pollock), deterministic"),
        exams=dict(
            V1=dict(what="start at |psi|^2 stays: max H <= 3 F (F = (cells-1)/(2N))", thr=3.0),
            V2=dict(what="M = 1 does not relax (|psi|^4 and half box): H_final/H_0 >= 0.9", thr=0.9),
            V3=dict(what="M = 16 and 25 relax: H_final - F < 0.1 (H_0 - F) for |psi|^4 and uniform starts", thr=0.1),
            V4=dict(what="faster with more modes: tau(1) > tau(4) > tau(16) >= tau(25) (|psi|^4 start, n = 48)"),
            V5=dict(what="law: log(H - F) linear in t until 5 % (R^2 >= 0.9), M = 16, 25", thr=0.9),
            V6=dict(what="continuum comparison, dimensionless tau * sigma_omega (rms spread of the mode frequencies, "
                         "per unit time) for M = 16, 25, |psi|^4 start: within a factor 2 of the continuum",
                    continuum=cont),
            V7=dict(what="lattice refinement: tau in physical units (ticks * a, a = theta) for M = 16, |psi|^4: n = 48 vs "
                         "96 agree within 30 %", thr=0.3)),
        report=["power-law alternative fit H - F ~ t^-p", "half-box start", "uniform start"]),
    part2=dict(
        title="Why <Sigma_z> = 0.51 and v_phi = 0.21 Za sin(theta) for the 3-D 1s (E5): the 4 valleys",
        hypothesis="each exact level of the 48^3 atom is 25 % on each of the 4 valleys (BATTIMENTO_1S); in valley v the "
                   "envelope is the k = 0 physics seen through the spin frame V_v (Bloch, above). Lab spin = sum over "
                   "valleys of w_v <Sigma_z>_v; the circulation of the lab current partly cancels.",
        exams=dict(
            S1=dict(what="Ritz (filter diagonalisation, even ticks, 1200 ticks, 25 frequencies 0.82..0.94) reproduces the "
                         "multiplet 0.8404 / 0.8863 / 0.8993 / 0.8996", thr=1e-3),
            S2=dict(what="valley weights of the exact 1s (0.8996) and of the E5 state: each valley 0.25 +- 0.05 for the exact "
                         "level", thr=0.05),
            S3=dict(what="valley decomposition explains the lab spin: |sum_v <Sigma_z>_v (lab) - <Sigma_z>| <= 0.05 and "
                         "|<Sigma_z>_v in its own frame (V_v^dag Sigma_z V_v)| >= 0.9 for every valley with weight > 0.1",
                    thr=dict(sum=0.05, frame=0.9)),
            S4=dict(what="control with valley mixing switched off: nucleus smeared to -Za erf(r/R)/r, R = 1 lambda_C "
                         "(2 cells; a DIAGNOSTIC, not physics), 1s listened from the radial seed (2 passes of 4096 ticks): "
                         "k = 0 valley weight >= 0.9, <Sigma_z> >= 0.9, two-tick-flow surfer v_phi/(Za sin) in [0.75, 1.25] "
                         "(the E5 window), measured as in E5",
                    thr=dict(valley0=0.9, sigma=0.9, vphi=[0.75, 1.25])),
            S5=dict(what="E5 redone on the exact level 0.8996 (two-tick flow, 2e5 surfers, 1000 ticks): report v_phi, "
                         "|z|, drift; threshold of E5 kept: v_phi/(Za sin) in [0.75, 1.25]", thr=[0.75, 1.25])),
        note="a smaller lattice spacing (a = 0.25, 96^3) is NOT run: the valley mixing comes from the point nucleus on "
             "the cell corner and does not go away with a (BATTIMENTO_1S); S4 tests the cause directly."),
)
with open(os.path.join(OUT, "surf_f2_prereg.json"), "w") as f:
    json.dump(P, f, indent=1)
print("written lab/out/surf_f2_prereg.json")
