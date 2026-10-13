"""aggmol_h2_gate_posthoc.py - POST HOC (after the G1 FAIL; not pre-registered): the same gate with larger bases, to tell
a basis-size limit from a bug. Writes lab/out/aggmol_h2_gate_posthoc.json"""
import json, os, aggmol_ecg as E, aggmol_exact as ex
M = E.M_KIN
out = {}
for K in (14, 24, 34):
    e, B, t = E.build(1, 0.0, K, ntrial=40, nref=8, budget=150)
    out[f"H_K{K}"] = e + M / 2
key = M * ex.energy(M * 2.0)[0]
for K in (40, 80):
    e, B, t = E.build(1, 2.0, K, ntrial=40, nref=8, budget=200)
    out[f"H2p_K{K}"] = e - key
print(out)
json.dump(out, open(os.path.join(E.OUT, "aggmol_h2_gate_posthoc.json"), "w"), indent=1)
