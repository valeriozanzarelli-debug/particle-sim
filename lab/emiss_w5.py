"""emiss_w5.py - test W5 (lab/out/emiss_prereg3.json): what splits the 2p3/2 doublet 0.97872/0.97887?
Valley weights of each member (Fourier, k reduced mod pi per component), the radial start for reference, and the
overlap of the members found from the m = 3/2 and m = 1/2 starts.   Run: python lab/emiss_w5.py -> lab/out/emiss_w5.json"""
import json, os
import numpy as np
from emiss_core import CK, OUT, norm2, ovl
from emiss_v4_dev import make_atom, start_state, S_LEARNED

L = json.load(open(os.path.join(OUT, "emiss_levels5.json")))


def valleys(psi):
    n = psi.shape[-1]
    F = np.fft.fftn(psi, axes=(1, 2, 3))
    P = np.sum(np.abs(F)**2, 0)
    k = 2 * np.pi * np.fft.fftfreq(n)
    km = np.abs((k + np.pi / 2) % np.pi - np.pi / 2)          # distance from 0 mod pi, in [0, pi/2]
    near0 = km < np.pi / 4
    a, b, c = np.meshgrid(near0, near0, near0, indexing="ij")
    out = {"0": P[a & b & c].sum(), "xy": P[~a & ~b & c].sum(), "xz": P[~a & b & ~c].sum(), "yz": P[a & ~b & ~c].sum()}
    tot = P.sum()
    out = {k_: float(v / tot) for k_, v in out.items()}
    out["other"] = 1 - sum(out.values())
    return out


D = make_atom(n=48, s=S_LEARNED)
res = {}
mem = {}
for st in ("2p3/2+3", "2p3/2+"):
    tag = f"48_{st.replace('/', '_')}_2400_0.0025"
    top = L[tag]["top"][:2]
    mem[st] = []
    for t_ in top:
        v = np.load(os.path.join(CK, f"ritzvec5_{tag}_{t_['index']}.npy"))
        v = v / np.sqrt(norm2(v))
        mem[st].append(v)
        vw = valleys(v)
        res[f"{st} member {t_['omega']:.5f}"] = dict(omega=t_["omega"], weight_in_start=t_["weight"], valleys=vw)
        print(f"{st} member {t_['omega']:.5f} (weight {t_['weight']:.3f}): valleys", {k: round(x, 3) for k, x in vw.items()}, flush=True)
    s0 = start_state(D, st)
    res[f"{st} radial start"] = dict(valleys=valleys(s0 / np.sqrt(norm2(s0))))
    print(f"{st} radial start: valleys", {k: round(x, 3) for k, x in res[f'{st} radial start']['valleys'].items()})
# decision rule (pre-registered)
mvals = [v["valleys"] for k, v in res.items() if "member" in k]
is_valley = all(sum(1 for key in ("0", "xy", "xz", "yz") if mv[key] >= 0.15) >= 3 for mv in mvals)
not_valley = all(mv["0"] >= 0.80 for mv in mvals)
res["W5_decision"] = "VALLEYS" if is_valley else ("NOT VALLEYS" if not_valley else "UNDECIDED")
print("W5 decision:", res["W5_decision"])
json.dump(res, open(os.path.join(OUT, "emiss_w5.json"), "w"), indent=1)
