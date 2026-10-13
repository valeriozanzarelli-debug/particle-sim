"""aggfull_p6_analyze.py - verdicts of PHASE 6 (lab/out/aggfull_p6_prereg.json) -> lab/out/aggfull_p6_verdict.json"""
import os, json
import numpy as np
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "aggfull_p6_prereg.json")))
E = P["exams"]
L = json.load(open(os.path.join(OUT, "aggfull_p6_learn.json")))
A = json.load(open(os.path.join(OUT, "aggfull_p6_atom.json")))
G = json.load(open(os.path.join(OUT, "aggfull_p6_glue_err.json")))
V, rep = [], {}


def check(n, ok, msg):
    V.append((n, bool(ok), msg)); print(f"[{'PASS' if ok else 'FAIL'}] {n}: {msg}", flush=True)


lv = L["levels"]
ph = {l["B"]: {k: v[0] if isinstance(v, list) else v for k, v in l["physics"].items()} for l in lv}
phi = {l["B"]: {k: v[1] if isinstance(v, list) else 0 for k, v in l["physics"].items()} for l in lv}
ok = all(abs(p["a"] + 0.25) <= E["L1"]["tol"] and abs(p["b"]) <= E["L1"]["tol"] and abs(p["c"] - 0.25) <= E["L1"]["tol"]
         and abs(p["gso"] - 0.25) <= E["L1"]["tol"] for p in ph.values())
check("L1 LC coefficients a, b, c, gso vs Dirac (0.005)", ok, "; ".join(
    f"B={B}: a {p['a']:+.5f} b {p['b']:+.5f} c {p['c']:+.5f} gso {p['gso']:+.5f}" for B, p in ph.items()))
ok = all(abs(l["N_of_k"]["n2"] / -0.125 - 1) <= E["L2"]["tol"] and l["N_of_k"]["nonscalar_max"] <= E["L2"]["nonscalar"] for l in lv)
check("L2 large-component amplitude N(k)", ok, "; ".join(
    f"B={l['B']}: n2 {l['N_of_k']['n2']:+.6f} (-0.125), non-scalar {l['N_of_k']['nonscalar_max']:.1e}, vs sqrt((E+1)/2E) {l['N_of_k']['dirac_rel_dev']:.1e}" for l in lv))
t4 = E["L3"]["theta_over_4"]
ok = all(abs(phi[B]["a"] / t4 - 1) <= E["L3"]["tol"] and abs(-phi[B]["b"] / t4 - 1) <= E["L3"]["tol"] for B in phi)
check("L3 imaginary parts = +-theta/4", ok, "; ".join(f"B={B}: Im a {phi[B]['a']:+.6f}, Im b {phi[B]['b']:+.6f} (theta/4 {t4:.6f})" for B in phi))
rep["parity_odd_rel"] = {l["B"]: l["parity_odd_rel"] for l in lv}

q = A["eT8_2p_c"]["LC"]
check("A1 2p alpha^4 (a0/8, T + LC) within 5 %", abs(q["err12"]) <= 0.05 and abs(q["err32"]) <= 0.05,
      f"2p1/2 {100 * q['err12']:+.2f} %, 2p3/2 {100 * q['err32']:+.2f} %")
s = {j: A[j]["LC"] for j in ("eT16_1s", "eT8_1s", "eT8_2s")}
check("A2 s alpha^4 (T + LC) within 5 %", all(abs(x["err"]) <= 0.05 for x in s.values()),
      "; ".join(f"{j}: {100 * x['err']:+.1f} % (Darwin part {100 * (x['darwin_part_eV'] / x['darwin_dirac_eV'] - 1):+.1f} %, kinetic "
                f"{100 * (A[j]['kinetic_remainder_eV'] / x['kinetic_dirac_eV'] - 1):+.1f} %)" for j, x in s.items()))
diffs = {}
for j, a in A.items():
    if a["state"] == "2p":
        mv = -(7.2973525693e-3 ** 4 * 510998.95 / 128) * (8 / 1.5 - 3)
        d = max(abs(a["LC"]["shift12_eV"] - a["FW"]["shift12_eV"]), abs(a["LC"]["shift32_eV"] - a["FW"]["shift32_eV"]))
    else:
        mv = a["LC"]["kinetic_dirac_eV"]
        d = abs(a["LC"]["shift_eV"] - a["FW"]["shift_eV"])
    diffs[j] = d / abs(mv)
check("A3 language independence LC vs FW (1 % of mass-velocity)", all(v <= E["A3"]["tol"] for v in diffs.values()),
      "; ".join(f"{j}: {100 * v:.2f} %" for j, v in diffs.items()))

h_over = 1 / 8                                                      # cell in a0/Z for both pairs
pred, ratios = {}, {}
for n, js in ((1, (5, 8, 12, 16)), (2, (5, 8, 12))):
    for J in js:
        x = (J - 0.5) * h_over
        p = 4 * np.exp(-2 * x) if n == 1 else 1.6 * (1 - x / 2) ** 2 * np.exp(-x)
        d = G[f"lc_J{J}_{n}s"] - G[f"base_J{J}_{n}s"]
        ratios[f"{n}s_J{J}"] = d / p
check("G1 learned outer term = outer spread Darwin (30 %)", all(abs(r - 1) <= E["G1"]["tol"] for r in ratios.values()),
      "ratios measured/predicted " + ", ".join(f"{k} {v:.3f}" for k, v in ratios.items()))
rep["G1"] = ratios
big = (24, 32, 48, 64)
check("G2 plateau cured (|err| <= 5 % at J0 >= 24)", all(abs(G[f"lc_J{J}_{n}s"]) <= 0.05 for J in big for n in (1, 2)),
      "; ".join(f"J{J}: 1s {100 * G[f'lc_J{J}_1s']:+.1f} % (base {100 * G[f'base_J{J}_1s']:+.1f}), 2s {100 * G[f'lc_J{J}_2s']:+.1f} % "
                f"(base {100 * G[f'base_J{J}_2s']:+.1f})" for J in big))
check("G3 J0 = 3 with the learned term (5 %)", abs(G["lc_J3_1s"]) <= 0.05 and abs(G["lc_J3_2s"]) <= 0.05,
      f"1s {100 * G['lc_J3_1s']:+.1f} % (base {100 * G['base_J3_1s']:+.1f}), 2s {100 * G['lc_J3_2s']:+.1f} % (base {100 * G['base_J3_2s']:+.1f})")
rep["glue_table"] = G
rep["verdicts"] = [dict(exam=n, passed=o, detail=m) for n, o, m in V]
json.dump(rep, open(os.path.join(OUT, "aggfull_p6_verdict.json"), "w"), indent=1, default=float)
print(f"{sum(o for _, o, _ in V)} PASS / {sum(not o for _, o, _ in V)} FAIL")
