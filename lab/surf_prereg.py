"""surf_prereg.py -> lab/out/surf_prereg.json : pre-registration of the study SURFISTI_ESATTI (written and committed
BEFORE the production runs of lab/surf_run.py). Development done before this file (declared): surf_core.py self tests
(walk with flux == zz_core, 1-D quantile 7e-16), surf_dev.py (random 12^3 state: chi2 ok, reversibility,
Liouville over 1-5 ticks), timings, and SMOKE runs of every mode of surf_run.py with tiny lengths (10-40 ticks, outputs
*_smoke.json, not used for any exam). Seen in the smoke runs: old reader on the free 3-D packet z = 190 after 12 ticks;
1-D lattice-vs-continuum quantile distance ~0.8 a at t = 1 (first order). The 3-D listened states (surf_states.py)
were still being computed when this file was written.
"""
import json
import os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = dict(
    study="SURFISTI_ESATTI", prefix="surf_", date="2026-10-13",
    goal="a surfer rule on the zig-zag lattice that keeps |psi|^2 EXACTLY (machine precision) with the least randomness",
    schemes=dict(
        flow="S4 NEW: deterministic flow in flow time (Pollock / lowest-order Raviart-Thomas) of the NET link flux of the tick",
        flow2="S4 over the two-tick period of the 3-D rule (x,y,z then z,y,x): net flux of the pair, rho(2m) -> rho(2m+2)",
        bellflux="S2b: Bell continuous-time jumps with rates F+/rho(s) (net flux of the tick, flow time), exact master equation",
        bell="S2: Bell minimal discrete jumps per sub-step (pre_core lambda = 0)",
        ballistic="S1: Bell gross-flow jumps per sub-step (pre_core lambda = 1): every surfer moves like the walk",
        cdfasse="S3: monotone CDF map along lines per sub-step (IDEE_RIPRESE R3)",
        old="S0: Rust zz.rs reader (trilinear spinor at the tick boundary, Bohm-Dirac v, midpoint) - not equivariant",
        old_rebirth="S0 + 0.2 % of surfers reborn from |psi|^2 every tick (as zzscena.rs: 40 of 20000)",
        proiettato="S5 (orbitals, continuous time): Dirac current at faces + minimal Poisson correction -> exact flux -> S4",
    ),
    free_parameters=dict(count=0, discrete_choices=[
        "net (minimal) flux instead of gross flows", "aggregation period (1 tick; 2 ticks for flow2)",
        "RT0 field inside a cell (each component linear in its own coordinate)", "rho(s) linear in time inside the tick",
        "S5 only: Poisson (L2-minimal) correction of the face current"]),
    exams=dict(
        E0a=dict(what="Walk3D with flux == zz_core.Walk; rho1 = rho0 - div F", thr=dict(walk=1e-14, div=1e-13), src="surf_core.py"),
        E0b=dict(what="Liouville over one tick on the 3-D production states (finite differences h = 1e-6)",
                 thr=dict(median=1e-5, q99=1e-3), runs=["a3d 1s", "a3d 1s2p", "a3d free"]),
        E0c=dict(what="reversibility of S4: 20 ticks forward then 20 back (with -F), cells", thr=dict(median=1e-12, q99=1e-9)),
        E1=dict(what="old F6 test: 1-D free walk theta = 0.1, R-Gaussian sigma = 5, 40 hbar/mc^2, 2000 surfers; max quantile error",
                thr=dict(flow=1e-12), report="old reader (expected ~1.26e-2 as in ZIGZAG F6b)", run="f6"),
        E2=dict(what="continuum limit 1-D: max_u |Q_lattice - Q_continuum| over t <= 40 for a = theta = 0.2, 0.1, 0.05",
                thr=dict(monotone=True, order_min=0.8, bohm_identity=1e-4),
                report="Bell (S2b) rms distance from S4 at t = 40 vs a (Vink: expected to shrink)", run="cont"),
        E3=dict(what="hydrogen Z = 1 (alpha) 1s of the radial rule, a = 0.5, 1e4 ticks, 2e4 surfers: stationary real state",
                thr=dict(flow_max_dq=1e-9, jumps_per_surfer_per_tick_bell_and_bellflux=1e-6, fraction_ever_jumped=1e-2),
                report="ballistic, old", run="h1s"),
        E4=dict(what="radial Z = 60, a = 0.05, (1s + 2s)/sqrt2, 20 beat periods: equivariance",
                thr=dict(flow_quantile_err=1e-10, max_abs_z=4.0, schemes=["flow", "bellflux", "bell", "ballistic"]),
                report="old", run="beat"),
        E5=dict(what="3-D zig-zag atom Z = 60, n = 48, a = 0.5 (the live atom of zz.rs), listened 1s, 1000 ticks",
                thr=dict(max_abs_z=4.0, tv_ratio=1.2, schemes=["flow2", "flow", "bellflux", "bell", "ballistic", "cdfasse"],
                         flow2_drift_growth_1000_over_100=1.5, flow2_neq_H_ratio_min=0.9,
                         flow2_vphi_over_Za_sin=[0.75, 1.25], flow2_vr_rms_over_Za_max=0.25),
                note="continuum Bohm-Dirac for the exact 1s: v_phi = Z alpha c sin(theta), v_r = 0", run="a3d 1s"),
        E6=dict(what="3-D atom, (1s + 2p1/2)/sqrt2 (listened), 1670 ticks (~10 beat periods)",
                thr=dict(max_abs_z=4.0, tv_ratio=1.2, schemes=["flow2", "flow", "bellflux"]), report="others", run="a3d 1s2p"),
        E7=dict(what="3-D free packet (sigma = 3, k = (0.6, 0.3, 0), theta = 0.5, periodic 48^3), 300 ticks",
                thr=dict(max_abs_z=4.0, tv_ratio=1.2, schemes=["flow2", "flow", "bellflux"]), report="others", run="a3d free"),
        E8=dict(what="fake mixing (c): random bits / jumps per surfer per tick, decorrelation of r, relaxation of a |psi|^4 start "
                     "in states where continuum de Broglie-Bohm does NOT mix (radial beat: 1-D order frozen; 3-D 1s: r conserved)",
                thr=dict(flow_neq_H_ratio_radial_min=0.999), report="all schemes", runs=["beat", "a3d 1s"]),
        E9=dict(what="analytic orbitals (as orbitali.rs): 1s + 2p1/2 continuous time on 48^3 at a = 0.5, dt = period/200, 2000 steps",
                thr=dict(max_abs_z_proiettato=4.0), report="old orbital reader with and without 0.5 %/step rebirth", run="orb"),
    ),
    statistics=dict(chi2="surfer histogram vs |psi|^2 on 2^3 blocks (3-D) / 4 or 8 cells (radial), bins with expected < 5 merged; "
                         "z = (chi2 - dof)/sqrt(2 dof); 41 checkpoints per run",
                    tv="total variation distance divided by the mean TV of an ideal iid sampler of the same N"),
    seeds=dict(base=20261013),
)
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "surf_prereg.json"), "w") as f:
    json.dump(P, f, indent=1)
print("written lab/out/surf_prereg.json")
