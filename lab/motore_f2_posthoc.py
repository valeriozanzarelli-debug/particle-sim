"""motore_f2_posthoc.py - DECLARED post-hoc check (after the judge of phase 2): is the +2.3 % (Ultra) / +3.2 % (Intel)
of P(2p0) after the kick a box effect? Intel level (a0/8) in boxes +-12, +-18, +-24 a0 (octant 96, 144, 192).
The verdict of G1b is NOT changed. -> lab/out/motore_f2_posthoc.json"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import motore_f2_core as f2
lev = f2.load_level3(512)
rows = []
for n in (96, 144, 192):
    ae, ao = f2.OctAtom(lev, n, "even"), f2.OctAtom(lev, n, "odd")
    v1, _, _ = f2.solve(ae, "1s"); v2, _, _ = f2.solve(ao, "2p")
    E1, E2 = ae.energy(v1)[0] * f2.MC2, ao.energy(v2)[0] * f2.MC2
    X, Y, Z = ae.grid()
    amp = np.vdot(v2, (np.sin(0.5 * Z / f2.A0) * v1.reshape((n,) * 3)).ravel())
    rows.append(dict(octant=n, box_a0=n / 8, E1s=E1, E2p=E2, line=E2 - E1, P2p=float(abs(amp) ** 2), P2p_rel=float(abs(amp) ** 2 / 0.073728 - 1)))
    print(rows[-1], flush=True)
json.dump(rows, open(os.path.join(f2.OUT, "motore_f2_posthoc.json"), "w"), indent=1)
