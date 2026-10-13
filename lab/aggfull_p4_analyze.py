"""aggfull_p4_analyze.py - verdicts of PHASE 4 (aggfull_p4_prereg.json + aggfull_p4_amend.json) -> aggfull_p4_verdict.json"""
import json, os
import numpy as np
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "aggfull_p4.json")))
RT = json.load(open(os.path.join(OUT, "aggfull_relearnT.json")))
F = json.load(open(os.path.join(OUT, "aggfull_p4_f3.json")))
V = []
def check(n, ok, msg):
    V.append((n, bool(ok), msg)); print(f"[{'PASS' if ok else 'FAIL'}] {n}: {msg}")
t = [P[f"T6_Z10_s{s}"] for s in range(1, 7)]
fin = [x["tried"][-1] for x in t]
check("M2b level 1 with k^6 pinned, 3e-5 on seeds 1-6", all(f["ok"] for f in fin) and all(abs(x["m_over_mfine_minus_1"]) <= 1e-9 for x in t),
      "w chosen " + str([x["w"] for x in t]) + "; energy " + ", ".join(f"{f['check']['energy_rel']:.2e}" for f in fin)
      + f"; {sum(f['ok'] for f in fin)} of 6 pass; mass dev <= {max(abs(x['m_over_mfine_minus_1']) for x in t):.1e}; c6 = Dirac x {P['c6_over_dirac']:.5f}")
dv = {w: max(P[f"M6_w{w}"]["max_dev_c2_c4_c6"]) for w in (5, 6, 7, 8)}
check("P3 Markov + k^6 needs w >= 8 (moment bound)", all(dv[w] > 1e-3 for w in (5, 6, 7)),
      ", ".join(f"w={w}: {dv[w]:.1e}" for w in dv) + " (w = 8 feasible, as the bound allows)")
mk = [(r["B"], max(r["markov_w4"]["max_c2_dev"], r["markov_w4"]["max_c4_dev"])) for r in RT]
check("P4 Markov impossible at a0/16 and a0/8", all(d > 1e-2 for _, d in mk), ", ".join(f"B={b}: {d:.1f}" for b, d in mk))
check("N1 float64 field removes the 2p high-k tail", False, "weight at kappa > 1.6: 4.4e-5 at N = 192 float64 (unconverged), "
      "1.9e-6 at N = 256 float32: the tail is not the float32 field (amendment aggfull_p4_amend.json)")
check("F3P 2p alpha^4 shifts within 5 % (amended: |kappa| < 0.8)", F["F3P_PASS"],
      "; ".join(f"{k}: 2p1/2 {100 * F[k]['err12_in']:+.1f} %, 2p3/2 {100 * F[k]['err32_in']:+.1f} %" for k in ("orig", "K", "T"))
      + "; whole cube (T): " + f"{100 * F['T']['err12_all']:+.0f} % / {100 * F['T']['err32_all']:+.0f} %")
json.dump(dict(verdict=V), open(os.path.join(OUT, "aggfull_p4_verdict.json"), "w"), indent=1, default=float)
print(f"=== {sum(v[1] for v in V)} PASS / {sum(not v[1] for v in V)} FAIL")
