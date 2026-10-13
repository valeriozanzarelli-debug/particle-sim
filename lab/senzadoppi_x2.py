"""senzadoppi_x2.py - exploration after the pre-registration (declared): the two lines of FW1 at theta = 0,5
(r4: 0,8844 and 0,8719) - valley content and spin, with the r7 machinery on 40^3, a = 0,5, 1600 periods, spin up and
spin along x.  Output lab/out/senzadoppi_x2.json"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import senzadoppi_r7 as r7

r7.A, r7.N, r7.NPER = 0.5, 40, 1600
r7.SCR = os.path.join(r7.SCR, "x2")
os.makedirs(r7.SCR, exist_ok=True)
for _ in range(10):
    try:
        r7.run("FW1")
    except SystemExit:
        continue
    break
st = json.load(open(os.path.join(r7.SCR, "senzadoppi_r7_FW1_ck.json")))
json.dump(st, open(os.path.join(os.path.dirname(__file__), "out", "senzadoppi_x2.json"), "w"), indent=1)
print(json.dumps(st))
