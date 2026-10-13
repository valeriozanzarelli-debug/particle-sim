"""emiss_x5b.py - declared extra diagnosis for X5 (lab/out/emiss_amend3.json): the true frequency of the weak dipole
lines.  The run records the dipole once per light tick (dt = 4: alias band 2 pi/4 = 1.571).  Here: final state of
P4_E1pi evolved freely, <z> on every EVEN tick (dt = 1), Blackman-Harris spectrum on [0, pi].
Run: python lab/emiss_x5b.py -> lab/out/emiss_x5b.json"""
import json, os, time
import numpy as np
from emiss_core import CK, OUT
from emiss_v4_dev import make_atom, S_LEARNED

t0 = time.time()
d = np.load(os.path.join(CK, "run_P4_E1pi.npz"))
psi = d["psi"].copy(); tk = int(d["t_ticks"])
D = make_atom(n=48, s=S_LEARNED)
NT = 6000
T, Zs = [], []
for t in range(NT):
    if (tk + t) % 2 == 0:
        T.append(t * D.a); Zs.append(D.position(psi)[2])
    D.tick(psi, tk + t)
T, Zs = np.array(T), np.array(Zs)
x = (T - T[0]) / (T[-1] - T[0])
h = 0.35875 - 0.48829 * np.cos(2 * np.pi * x) + 0.14128 * np.cos(4 * np.pi * x) - 0.01168 * np.cos(6 * np.pi * x)
om = np.arange(0.01, np.pi, 1e-4)
S = np.abs(np.exp(1j * np.outer(om, T)) @ (h * (Zs - Zs.mean()))) * 2 / np.sum(h)
pk = [(float(om[i]), float(S[i])) for i in range(1, len(om) - 1) if S[i] > S[i - 1] and S[i] > S[i + 1] and S[i] > 0.01 * S.max()]
pk.sort(key=lambda p: -p[1])
alias = [(w, a, float(abs((w + np.pi / 4) % (np.pi / 2) - np.pi / 4))) for w, a in pk]
res = dict(peaks=pk[:12], aliased_at_dt4=[(round(w, 5), round(f, 5)) for w, a, f in alias[:12]], seconds=time.time() - t0)
json.dump(res, open(os.path.join(OUT, "emiss_x5b.json"), "w"), indent=1)
print("dipole lines (omega, amplitude):", [(round(w, 5), round(a, 4)) for w, a in pk[:10]])
print("where they land when sampled every 4 time units:", res["aliased_at_dt4"][:10], f"({time.time() - t0:.0f} s)")
