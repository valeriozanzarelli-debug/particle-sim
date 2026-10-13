"""motore_amend.py - DECLARED post-hoc check (amendment A1, written after the B1 run): the two flagged B1 packets
(N = 64, sigma = 2 lambda) were run in a box of +-4.5 sigma = +-9 lambda, but over 100 lambda the packet spreads by
diffraction to ~ +-9 lambda and its tails wrap around the periodic box (the same artefact seen in development).
Here the same two cases in a box of +-7 sigma. The B1 verdict is NOT changed: this file only says whether the flag
sub-criterion failed because of the engine's model or because of the measurement. -> lab/out/motore_amend.json"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import motore_core as mc
from motore_m1 import run_pair, screen_meas
rows = []
for p in (0.1, 0.2):
    r = run_pair(64, 30.0, 100.0, 2.0, box_sig=7.0)
    e, env, fr, res_ = screen_meas(r, p)
    r.update(footprint=p, meas_px=e, env_px=env, fringe_px=fr)
    rows.append(r)
    print(p, r["n"], r["meas_dx"], r["meas_phase"], "px", e, flush=True)
json.dump(dict(amendment="A1 post-hoc box check for the flagged B1 packets", rows=rows), open(os.path.join(mc.OUT, "motore_amend.json"), "w"), indent=1)
