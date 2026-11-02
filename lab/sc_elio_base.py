"""sc_elio_base.py - SCENE_CALCOLATE: rebuild the helium ECG basis (the .npz was never in git). Same machine, same seed, same
size as multi_e_ecg_prod.py (job he1S_Z2): the energy must come back (exam S2). Writes lab/out/sc_he_basis.npz in the same
form as the H2 bases of h2v_wf.py (A, s = 0, nrm, c, E_el), lengths already scaled to the phase-5 mass m0.

    cd lab && python sc_elio_base.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("OMP_NUM_THREADS", "1")
from multi_e_core import OUT  # noqa: E402
from multi_e_ecg import System, run_svm  # noqa: E402

m_A = json.load(open(os.path.join(OUT, "multi_e_kernel.json")))["R3"]["m_kin_theta_A_used"]
m0 = json.load(open(os.path.join(OUT, "multi_e_v5.json")))["checks"]["R3v5"]["m0"]
ref_A = json.load(open(os.path.join(OUT, "multi_e_ecg.json")))["prod"]["nr"]["he1S_Z2"]["E"]
s = System(2, 1, 2.0, m_A, name="he1S_Z2")
sv, _, _ = run_svm(s, 200, 2.0, 2026, Lmax=20.0, log=lambda *a: None)
E_A = sv.rebuild()
E, c = sv.state()
E_m0 = E_A * m0 / m_A
print(f"E(m_A) {E_A:.10f} registrata {ref_A:.10f} diff {E_A - ref_A:.2e}; E(m0) {E_m0:.8f} esportazione -2.90372128 diff {E_m0 + 2.9037212829955945:.2e}", flush=True)
np.savez(os.path.join(OUT, "sc_he_basis.npz"), A=sv.A * (m0 / m_A) ** 2, s=np.zeros((len(sv.A), 2, 3)), nrm=sv.nset, c=c,
         E_el=E_m0, E_A=E_A, ref_A=ref_A)
print("ALL DONE", flush=True)
