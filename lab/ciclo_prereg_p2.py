"""ciclo_prereg_p2.py - phase 2 entry of the CICLO_LIVELLI pre-registration (imported by ciclo_prereg.py).
Committed BEFORE the production runs of phase 2. Development runs (declared, numbers seen):
  * small 1-D chains (base 233-729 cells): first design (slow set = electron band up to half Nyquist, des Cloizeaux
    on the n C slowest modes) gave nonlocal rules (errors 1e-2 at w = 6); cause: the Dirac regime (cells 0.2 lambda_C:
    the positive band alone is nonlocal below lambda_C) and the coarse-zone edge. Changes, all before this file:
    common BASE level in the Pauli regime (B0 = 32 fine cells of theta = 0.05: 1.6 lambda_C, stencil w = 8 with
    G(0) = 0 and the kinetic mass transported EXACTLY, phase-1 lesson); sector = the n C modes the coarse variables
    capture best; relative least squares (as learn_h); check window = >= 24 coarse cells per wavelength (project HD
    rule) and >= 9 modes; slow set (variables, weight 1) = 12 cells per wavelength.
  * seen in development: periodic b = 2 and b = 3 chains pass every level with C = 1, w = 1-4, accumulated error
    1e-3 .. 4e-3 on 512/729 cells; the Fibonacci chain passes level 1 (w = 4-6) but NO candidate at level 2 (w <= 6);
    the exact (untruncated) Fibonacci coarse operator has irregular tails ~0.04 at 6-8 cells. Therefore wmax is raised
    to 8 for ALL schemes (declared before production)."""

PHASES = dict(fase2=dict(
    title="Fase 2: variabili grosse automatiche + fattore di blocco 2 / 3 / phi (Fibonacci)",
    date="2026-10-13",
    config_1d=dict(
        fine="1-D zig-zag walk theta = 0.05 (cell 0.05 lambda_C)",
        base="B0 = 32 fine cells (1.6 lambda_C), C = 1, real even stencil w = 8, G(0) = 0 and c2 exact (constraints), "
             "weights 1 (|kappa| <= pi/2), 0.03 beyond",
        schemes={"2": "periodic b = 2, 2048 base cells", "3": "periodic b = 3, 2187 base cells",
                 "fib": "Fibonacci blocks (A <- AB, B <- A), 2584 base cells (F_18), mean factor phi per level"},
        loop="variables: block density matrix of the slow set (>= 12 coarse cells per wavelength), C = fewest numbers "
             "(1..2) that pass; sector = the n C modes best captured; rule: Hermitian generator tied by local pattern "
             "of block types, |j| <= w, w = smallest in 1..8 that passes; relative least squares (floor 1e-2 of the "
             "check scale), fast modes weight 0.03; check: lowest max(9, #modes with >= 24 cells/wavelength) eigenvalues "
             "vs the level below, budget 1e-3 of the check kinetic scale, captured weight of the check set >= 0.98; "
             "repeat from the LEARNED rule; stop when the next level has < 48 cells",
        observables="free: base eigenvalues with kappa_base <= 2 pi 4 / 2048 (9 / 9 / 11 modes), accumulated error "
                    "max|d omega| / (range); coulomb: soft Coulomb V = -g/sqrt(x^2 + s^2), g = 0.015, s = 10 lambda_C, "
                    "base V = block mean of the fine V, coarse V = Galerkin W V W^dag; 3 lowest bound states (binding "
                    "energies) against the exact base",
        fine_anchor="bound states of the fine zig-zag (N = 32 x 2048 fine cells, sparse shift-invert) vs the base"),
    config_3d=dict(
        start="existing chain-A level B = 16 (lab/out/agg_zz_3d_A.pkl, scalar part; spin part 1e-10)",
        schemes="B_{n+1} = b B_n with b = 2, 3, phi (phi: non-integer rescaling, possible in k-space because the "
                "generator is transported as b G(kappa / b); NO quasi-periodic structure in 3-D: declared)",
        learn="real even stencil, cube |j|_inf <= w, w = smallest in 1..5 that passes; samples 75 % like the 1s "
              "(k_s = alpha m) + 25 % uniform with weight 0.03, relative least squares; G(0) = 0 and the axis c2 "
              "transported exactly from the level below (phase-1 lesson)",
        check="energy error weighted by the 1s momentum density vs the level below <= 3e-4; bz_rel90 <= 0.2; no negative "
              "motion energy (<= 1 %) -- the chain-A budgets",
        observables="accumulated error vs the fine walk's exact positive band; 1s energy shift at first order "
                    "dE = <(G - G_fine)/tau>_1s relative to E_1s = alpha^2/2; chain until B >= 600"),
    exams=dict(
        A1=dict(what="1-D: at every admitted level of every scheme the automatic block variables lie in the hand "
                     "descriptor subspace (block mean of the large component / of the level-below amplitude): "
                     "overlap >= 0.99; C = 1 everywhere", tol=0.99, prediction="PASS"),
        A2=dict(what="3-D: the 2 leading eigenvectors of the positive-band projector of the fine walk averaged over the "
                     "1s momentum density span the large component (L + R)/sqrt 2 (x) spin: overlap >= 0.999 (the "
                     "machine rediscovers 'pauli2' without being told)", tol=0.999, prediction="PASS"),
        A3=dict(what="1-D free, periodic b = 2 and b = 3: every level admitted down to the stop, last B >= 27, "
                     "accumulated error at the last level <= 1e-2", tol=1e-2, prediction="PASS"),
        A4=dict(what="1-D free, Fibonacci: the same (every level admitted, last B >= 27, accumulated <= 1e-2)",
                tol=1e-2, prediction="FAIL (development: level 2 not admitted at w <= 6)"),
        A5=dict(what="1-D free, blocking score = accumulated error x taps (2w+1) C^2 at B = 8 and B = 16 (log "
                     "interpolation between levels): phi WINS if its score is lower than the best of {2, 3} by >= 20 % "
                     "at both B, LOSES if higher by >= 20 % at either B or if its chain does not reach B, else TIE",
                prediction="LOSES (known physics gives no reason for phi; quasi-periodic rules are less local)"),
        A6=dict(what="1-D soft Coulomb: binding energies of the 3 lowest bound states at the level nearest B = 4 within "
                     "1e-2 of the exact base, for each scheme that reaches it; score as A5 at B = 4", tol=1e-2,
                prediction="PASS for b = 2, 3; phi as A4"),
        A7=dict(what="3-D free: each scheme reaches B >= 400 with every level admitted, and the first-order 1s shift "
                     "|dE|/E_1s at B = 400 (log interpolation) <= 1e-3; score = |dE|/E_1s x taps; same win/lose rule",
                tol=1e-3, prediction="PASS; phi TIE or LOSES"),
        A8=dict(what="anchor: base bound states vs the fine zig-zag bound states within 1e-2 (binding energies)",
                tol=1e-2, prediction="PASS"),
    ),
    noise_note="step 3 of the loop asks for a fluctuation-dissipation noise term: for ONE wave (no temperature) the "
               "discarded detail is a deterministic memory, measured here as the leak 1 - captured weight; FD noise "
               "needs a thermal many-body state (phase 4, the droplet): not tested in phase 2 (declared).",
))
