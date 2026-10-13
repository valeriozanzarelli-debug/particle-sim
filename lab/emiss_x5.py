"""emiss_x5.py - test X5 (lab/out/emiss_prereg3.json): where do the weak E4 lines at +0.0131 come from?
The final state of P4_E1pi (after the pulse, learned nucleus) evolved freely; spectrum of <psi_end|psi(t)> (even ticks, Hann).
Run: python lab/emiss_x5.py -> lab/out/emiss_x5.json"""
import json, os, time
import numpy as np
from emiss_core import CK, OUT, ovl, norm2
from emiss_v4_dev import make_atom, S_LEARNED

t0 = time.time()
d = np.load(os.path.join(CK, "run_P4_E1pi.npz"))
psi = d["psi"].copy()
tk = int(d["t_ticks"])
D = make_atom(n=48, s=S_LEARNED)
L4 = json.load(open(os.path.join(OUT, "emiss_levels4_48.json")))
w1s, w12, w32 = L4["1s+"]["omega"], L4["2p1/2+"]["omega"], L4["2p3/2+"]["omega"]
import sys
NT = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
WIN = sys.argv[2] if len(sys.argv) > 2 else "hann"
ref = psi.copy()
C = []
ck = os.path.join(CK, f"x5_{NT}.npz")
t_start = 0
if os.path.exists(ck):
    z = np.load(ck); psi = z["psi"]; C = [tuple(x) for x in z["C"]]; t_start = int(z["t"]); C = [(int(a.real), b) for a, b in C]
for t in range(t_start, NT):
    if time.time() - t0 > 520:
        np.savez(ck, psi=psi, C=np.array([(a, b) for a, b in C], complex), t=t)
        print("checkpoint", t, flush=True); sys.exit(3)
    if (tk + t) % 2 == 0:
        C.append((t, ovl(ref, psi)))
    D.tick(psi, tk + t)
tt = np.array([c[0] for c in C]) * D.a
cc = np.array([c[1] for c in C])
x_ = np.arange(len(cc)) / (len(cc) - 1)
if WIN == "bh":   # Blackman-Harris 4-term: side lobes -92 dB
    h = 0.35875 - 0.48829 * np.cos(2 * np.pi * x_) + 0.14128 * np.cos(4 * np.pi * x_) - 0.01168 * np.cos(6 * np.pi * x_)
else:
    h = np.sin(np.pi * x_)**2
om = np.arange(0.80, 1.06, 2e-5)
S = np.abs(np.exp(1j * np.outer(om, tt)) @ (h * cc)) / np.sum(h)          # amplitude ~ weight of each level
pk = [(float(om[i]), float(S[i])) for i in range(1, len(om) - 1) if S[i] > S[i - 1] and S[i] > S[i + 1] and S[i] > 1e-4]
pk.sort(key=lambda p: -p[1])
win1 = [p for p in pk if 0.880 <= p[0] <= 0.892]
win2 = [p for p in pk if 0.985 <= p[0] <= 0.996]
W1 = max([p[1] for p in win1], default=0.0)
W2 = max([p[1] for p in win2], default=0.0)
res = dict(peaks=pk[:15], levels=dict(w1s=w1s, w12=w12, w32=w32), H1_weight=W1, H2_weight=W2,
           H1_peaks=win1[:3], H2_peaks=win2[:3], decision="H1 (1s-type level near w1s - 0.0131)" if W1 > W2 else "H2 (2p-type level near w2p + 0.0131)",
           predicted_lines=dict(from_H1=[w32 - (w1s - 0.0131), w12 - (w1s - 0.0131)]), seconds=time.time() - t0)
res["NT"] = NT; res["window"] = WIN
json.dump(res, open(os.path.join(OUT, "emiss_x5.json" if NT == 4000 else f"emiss_x5_{NT}_{WIN}.json"), "w"), indent=1)
if os.path.exists(ck): os.remove(ck)
print("peaks (omega, weight):", [(round(a, 5), round(b, 4)) for a, b in pk[:12]])
print(f"H1 window [0.880, 0.892]: {win1[:3]};  H2 window [0.985, 0.996]: {win2[:3]}  -> {res['decision']}  ({time.time() - t0:.0f} s)")
