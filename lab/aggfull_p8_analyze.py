"""aggfull_p8_analyze.py - verdicts of PHASE 8 -> lab/out/aggfull_p8_verdict.json"""
import os, json
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "aggfull_p8_prereg.json")))
E = json.load(open(os.path.join(OUT, "aggfull_p8_err.json")))
J = P["config"]["J0_map"]
V = []


def check(n, ok, msg):
    V.append((n, bool(ok), msg)); print(f"[{'PASS' if ok else 'FAIL'}] {n}: {msg}", flush=True)


def row(Z, js):
    return "; ".join(f"Z{Z}/{2 * int(Z)} J{j}: 1s {100 * E[f'p8_P{Z}_J{j}_1s']:+.2f} %, 2s {100 * E[f'p8_P{Z}_J{j}_2s']:+.2f} %" for j in js)


ok = lambda Z, js: all(abs(E[f"p8_P{Z}_J{j}_{n}s"]) <= 0.01 for j in js for n in (1, 2))
check("Q1 pair (1,2), J0 = 32, 48, 64 within +-1 %", ok("1", J["1"][1:]), row("1", J["1"]))
check("Q2 pairs (3,6), (5,10) at R >= 4.06 a0/Z", ok("3", J["3"][1:]) and ok("5", J["5"][1:]), row("3", J["3"]) + " | " + row("5", J["5"]))
check("Q3 pairs (3,6), (5,10) at R >= 6.06 a0/Z", ok("3", J["3"][2:]) and ok("5", J["5"][2:]), "see Q2")
q4 = {Z: E[f"th1_P{Z}_J{J[Z][-1]}_1s"] / (-0.1187 / int(Z)) for Z in J}
check("Q4 unpatched chain: 1s error = -11.87 %/Z (30 %)", all(abs(r - 1) <= 0.30 for r in q4.values()),
      "; ".join(f"Z{Z}: {100 * E[f'th1_P{Z}_J{J[Z][-1]}_1s']:+.2f} % (pred {-11.87 / int(Z):+.2f}), 2s {100 * E[f'th1_P{Z}_J{J[Z][-1]}_2s']:+.2f} %" for Z in J))
json.dump(dict(err=E, verdicts=[dict(exam=n, passed=o, detail=m) for n, o, m in V]), open(os.path.join(OUT, "aggfull_p8_verdict.json"), "w"), indent=1)
print(f"{sum(o for _, o, _ in V)} PASS / {sum(not o for _, o, _ in V)} FAIL")
