"""zstar_amend.py - AMENDMENTS to lab/out/zstar_prereg.json, declared after seeing data. Writes lab/out/zstar_amend.json."""
import json, os
A = {
    "A1": {
        "when": "after the production runs of S0.25 and R+i at a = 3e-4, 1e-4",
        "seen": "the pre-registered tracker (DZMAX 0.5, start Z*-20) returned Z_cr = 146.23 / 144.43 (S0.25) and "
                "139.25 (R+i) instead of ~103.7 / ~95.2: for 2 gamma* > 1 the dive happens in a window "
                "|Z - Z*| ~ a^(2 gamma*) (0.01 in Z at a = 1e-4) much narrower than the step; the tracker jumped "
                "over it and followed the next level to its own (zeta > 1) dive",
        "change": "extra runs 'VARIANT:fine' = same tracker with DZMAX = 0.01 started at Z0 = Z*_formula - 1 for "
                  "S0.25 and R+i at all a. The pre-registered values are kept in zstar_runs.json and reported; the "
                  "exams are judged on BOTH (pre-registered result = FAIL of the tracker, declared)",
    },
    "A2": {
        "when": "after the production runs of K-2 (p3/2) and K+2 (d3/2)",
        "seen": "K-2: the tracker (start 180, DZMAX 0.5) returned 285.50, 282.07, 279.79, 278.46 (a = 3e-3 ... 1e-4): "
                "the p3/2 level's own dive just above zeta = 2 (Z = 274, the kappa = -2 fall to the centre), i.e. it "
                "jumped over the predicted contact crossing at 206.02. K+2: the Z grid of step 0.1 found the crossing "
                "only at a = 1e-3 (42.0748); at 3e-4 and 1e-4 the crossing window is narrower than the grid",
        "change": "K-2:fine = same tracker with DZMAX = 0.002 from Z*-0.3; K+2:fine = grid step 0.002 in [41.85, 42.35]. "
                  "Pre-registered results kept and reported",
    },
    "A3": {
        "when": "after A1 failed: the fine tracker at a = 3e-4 hit the 10-minute limit (S0.25:fine killed, no result)",
        "change": "new method lab/zstar_run_sign.py: at omega = -1 the walk's own regular solution from r = 0 (same rule, "
                  "exact recursion) changes the sign of Im g(r_end = 20) exactly when a level crosses -mc^2; scan in Z "
                  "with step h + bisection to 1e-7. Validated BEFORE use on the base rule: 132.7834570, 132.0840766, "
                  "131.7123161, 131.5451387 vs the giro-1 tracker 132.7834570, 132.0840766, 131.7123161, 131.5451387 "
                  "(a = 3e-3 ... 1e-4). Applied to every variant (validation where the tracker worked) and extended to "
                  "a = 3e-5, 1e-5. K+2 (A2 replaced): crossing of omega = -1 and of omega = 0 by the same method. "
                  "Exams F2/F3/F4/K2/K3 are judged on these values; the pre-registered tracker values are reported too",
    },
}
json.dump(A, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "zstar_amend.json"), "w"), indent=1)
print(json.dumps(A, indent=1))
