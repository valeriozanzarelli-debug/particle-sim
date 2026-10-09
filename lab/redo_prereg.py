"""redo_prereg.py - PRE-REGISTRATION of the three re-tests R1-R3 (written and committed BEFORE any run).
Writes lab/out/redo_prereg.json.

DEV (declared): nothing measured yet. Only the existing results were read (pre_p1p5.json, flow_a.json).
"""
import json, os, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

R1 = {
    "idea": "Valerio: the vacuum is a sea of perfectly mixed pairs, invisible until a strong field breaks them (Schwinger).",
    "rule": "1+1 D zig-zag walk (zz_core restricted to one axis): R moves +1, L moves -1, then coin "
            "L' = e^{-iV}(cos th L - i sin th R), R' = e^{-iV}(cos th R - i sin th L). Dispersion cos w = cos th cos k. "
            "Mass m = th (rest quasi-energy). Units hbar = c = a = dt = 1, charge e = 1.",
    "vacuum": "filled lower Floquet band (quasi-energy in (-pi+th, -th)) = Dirac sea.",
    "R1a_uniform_field": {
        "how": "Peierls phase A(t) = E t on every hop (R picks e^{-iA}, L picks e^{+iA}): every sea plane wave stays a "
               "plane wave, k(t) = k + E t. Exact single-particle evolution of each sea state (2x2 per k). "
               "Pair probability per passage of the k = 0 gap: P0 = |<upper(k_end)|psi>|^2 for a state started in the "
               "lower band at k = -pi/2 and stopped at k = +pi/2 (adiabatic Floquet basis). Rate per length per time "
               "Gamma_lat = (E/2pi) P0 (E/2pi = number of k-modes per unit length crossing k = 0 per unit time).",
        "reference": "1+1 D spinor QED, MEAN NUMBER of pairs per length per time (Nikishov): Gamma = (eE/2pi) exp(-pi m^2/(eE)). "
                     "(The vacuum-decay rate 2 Im L = (eE/2pi) sum_n n^-1 exp(-n pi m^2/eE) is a different quantity.) "
                     "E_crit = m^2/e. Also a lattice Dykhne prediction P_D = exp(-(4/E) int_0^{s_c} arccos(cos th cosh s) ds), "
                     "cosh s_c = 1/cos th (the exact lattice dispersion continued to imaginary k), reported.",
        "scan": "th = 0.1 (main) and 0.05 (convergence); E/E_crit in {0.25, 0.5, 1, 2, 4} (pass range) plus "
                "{8, 16, 32, 64} (to find where the lattice breaks; no pass criterion).",
        "PASS": "|Gamma_lat/Gamma_Schwinger - 1| <= 0.10 for every E/E_crit in {0.25,0.5,1,2,4} at th = 0.1.",
        "real_space_check": "ring of N = 128 sites, Peierls phases in real space, the whole sea (128 states) evolved for one "
                            "half Bloch period: upper-band population per k equal to the k-space result to <= 1e-10.",
        "valid_window": "reported: the E/E_crit range where |ratio - 1| <= 0.10, and the time window t < pi/E (half Bloch period).",
    },
    "R1b_doublers_prediction": "PREDICTION (artefact): the walk has a second gap at quasi-energy pi (the pi-modes), with the "
                               "same size 2 th. A sea state crossing k = pi tunnels with P_pi = P0 (exact symmetry k->k+pi, "
                               "w->w+pi). So the RAW count of pairs over a full Bloch period is 2x Schwinger. CONFIRMED if "
                               "|P_pi/P0 - 1| <= 0.01 at every E of the pass range.",
    "R1c_sauter_step": {
        "how": "static potential V(x) = (V0/2)(1 + tanh(x/ell)) as the phase e^{-iV(x)} of the coin. Exact Floquet scattering "
               "at quasi-energy w by the lattice transfer matrix. Pairs (Klein zone) = transmission N(w) of a sea state "
               "through the step; pair rate Gamma_step = int_{Klein} dw/2pi N(w).",
        "params": "th = 0.1, V0 = 0.6 (= 6 m > 2 m), ell in {5, 10, 20, 40} sites (E_max/E_crit = V0/(2 ell th^2) = 6, 3, 1.5, 0.75).",
        "reference": "continuum 1+1 Dirac on the same tanh step, solved numerically (ODE, rtol 1e-10); the Sauter / "
                     "Gavrilov-Gitman closed form is quoted only if it matches the ODE to 1e-4.",
        "PASS": "|Gamma_step_lattice/Gamma_step_continuum - 1| <= 0.10 for ell >= 10.",
        "wavepacket_check": "one time-dependent run (lattice, gaussian packet of the sea in the Klein zone): transmitted "
                            "fraction within 0.02 of the transfer-matrix N averaged over the packet spectrum.",
    },
    "R1d_unruh": "optional; skipped unless cheap (declared now: probably skipped).",
    "parameters": "th (mass), E or (V0, ell): physical inputs. No fitted parameter.",
}

R2 = {
    "idea": "flux through the holes keeps the holes from collapsing (moduli stabilisation)",
    "start": "the SAME 'crumpled torus' of OCCHIO_E_FLUSSO case (a): hex torus L = 15, w = 1 + 0.4 U(-1,1), default_rng(20261050)",
    "flow": "flow_core unchanged: Ollivier kappa (alpha 0.5, exp(-w^2)), eps 0.5, tol 1e-3, maxit 100, normalised total length",
    "flux_term": "for each vertex v and each torus cycle class a in {1, 2} (lattice directions a1, a2): l_a(v) = length of the "
                 "shortest closed loop through v with ODD winding along a (Z2 double cover, Dijkstra on current lengths), "
                 "gamma_a(v) its edges. Flux energy E_flux = (g n^2 / Z0) sum_{a,v} L / l_a(v)  (energy ~ n^2/cycle length, "
                 "L = 15 the flat cycle length). Pressure on edge e: P_e = -L dE/dw_e = (g n^2/Z0) sum_{a,v: e in gamma_a(v)} "
                 "(L/l_a(v))^2, Z0 = mean number of loops per edge on the flat lattice (w = 1), so <P> ~ g n^2 when flat. "
                 "Step: w <- w - eps * w * (kappa - P), then the same normalisation. Flat (all w = 1) is still a fixed point.",
    "why_not_area": "an energy ~ n^2/area would be (almost) inert: the normalised flow keeps the total length fixed.",
    "coupling": "g = 0.1 fixed now (same order as the initial kappa spread 0.089). n in {0, 1, 2} -> g n^2 = 0, 0.1, 0.4. "
                "Extra (declared tuning probe, not an exam): n = 4 (g n^2 = 1.6).",
    "new_parameters": 1,
    "exams": {
        "one_piece_loops": "final MI read by the L1 eye (flow_run.read, THR2 from eye_calib): pieces == 1 AND loops == 2",
        "spread": "std(kappa_final) <= std(kappa_initial)/10",
        "spinor": "Regge holonomy N_inv == 2 on the final lengths",
    },
    "PASS": "the brake WORKS if n = 1 or n = 2 passes all three exams. PARTIAL if it passes one-piece+loops only. "
            "Also reported: fraction of edges with w > 2, length spread w_cv, versus n = 0.",
    "code_check": "n = 0 must reproduce flow_a.json (spread ratio 3.83, eye loops 0) up to LP round-off.",
    "diagnostic_D2": "superficial-error check of debt 9 (eps too big?): flat hex torus + 1% noise, 20 steps at eps = 0.5 and "
                     "100 steps at eps = 0.1 (same flow time 10), no flux. The instability is REAL (not a step-size "
                     "artefact) if the length spread w_cv grows in both.",
}

R3 = {
    "fix": "replace the random lattice jumps of pre_core (Markov chain -> fake mixing) by a DETERMINISTIC surfer: continuous "
           "position; during each axis sub-step of the walk the density along every lattice line changes from rho to "
           "rho' (the coin does not change rho). The surfer is moved by the monotone CDF map along its line, with the "
           "density taken uniform inside each site: x' = C'^-1(C(x)). In 1-D (per line, box walls anchor C) this map is "
           "exactly the flow with velocity j/rho of the Bell minimal net current (lam = 0), integrated over the sub-step: "
           "deterministic, no crossing, EXACTLY equivariant (a |psi|^2 crowd with uniform jitter stays |psi|^2).",
    "setup": "IDENTICAL to pre_p1p5 P1: pre_core.LatticeWalk box n = 24, th = 0.5, modes_state(omega0 = 0.9), 20000 surfers, "
             "cells 4x4, T = 400 ticks, seed 20261009; starts psi2 / uniform / half; M = 16 and M = 1.",
    "exams": {
        "R3a_equivariance": "start |psi|^2, M = 16: max_t H <= 3F",
        "R3b_relax": "M = 16, uniform and half starts: H(T)/H(0) < 0.10",
        "R3c_control": "M = 1 must NOT relax (uniform start is primary, half reported): tau(M=1)/tau(M=16) >= 3 "
                       "OR (H_M1(T)/H_M1(0) > 0.10 AND H_M16(T)/H_M16(0) < 0.03)",
    },
    "PASS": "R3a and R3b and R3c -> the old fake mixing (P1c) was a superficial artefact of the random jumps.",
    "secondary": "same exams at th = 0.125, n = 96 (cells 16x16, same physical box), reported, only if each run fits in 4 min.",
    "parameters": "0 new (the map has no knob).",
}

if __name__ == "__main__":
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=HERE, capture_output=True, text=True).stdout.strip()
    meta = dict(written=datetime.datetime.now().isoformat(timespec="seconds"), git_head_before=head)
    os.makedirs(OUT, exist_ok=True)
    json.dump(dict(meta=meta, R1=R1, R2=R2, R3=R3), open(os.path.join(OUT, "redo_prereg.json"), "w"), indent=1)
    print("written", os.path.join(OUT, "redo_prereg.json"))
