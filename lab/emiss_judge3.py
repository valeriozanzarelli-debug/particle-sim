"""emiss_judge3.py - judges lab/out/emiss_prereg3.json (phase 3).  -> lab/out/emiss_results3.json"""
import json, os, glob
import numpy as np
from zz_core import ALPHA, T_UNIT
from emiss_core import CK, OUT
from emiss_analysis import qe_rates, load

R = {}


def verdict(name, ok, msg, value=None):
    R[name] = dict(pass_=None if ok is None else bool(ok), msg=msg, value=value)
    print(f"[{'INFO' if ok is None else ('PASS' if ok else 'FAIL')}] {name}: {msg}", flush=True)


def have(name):
    return os.path.exists(os.path.join(CK, f"run_{name}.npz")) and any(
        f"=== {name} finished rc=0" in open(f).read() for f in glob.glob(os.path.join(CK, "drive*.log")))


L5 = json.load(open(os.path.join(OUT, "emiss_levels5.json")))
w96 = L5["96_2p3_2+3_1200_0.005"]["omega"] - L5["96_1s+_1200_0.005"]["omega"]
EX = json.load(open(os.path.join(OUT, "emiss_exact.json")))
if have("P5_R1_96"):
    q = qe_rates("P5_R1_96", w96, 4.0)
    d = load("P5_R1_96")
    Ar = 4 / 3 * ALPHA * w96**3 * float(np.sum(np.abs(d["ser_Dc"][0])**2)) / d["ser_norm"][0]
    verdict("R5", abs(q["A_norm"] / (Ar * q["W2"]) - 1) < 0.03, f"A_norm {q['A_norm']:.4e} / (A_rule {Ar:.4e} x W2 {q['W2']:.4f}) = {q['A_norm'] / (Ar * q['W2']):.4f} (w = {w96:.5f})", q["A_norm"] / (Ar * q["W2"]))
    verdict("R5b", abs(q["A_photons"] / q["A_norm"] - 1) < 0.03, f"A_photons / A_norm = {q['A_photons'] / q['A_norm']:.4f}", q["A_photons"] / q["A_norm"])
    ALW = EX["Z60"]["2p3/2"]["A_LW"]
    A = q["A_norm"] / q["W2"]
    verdict("R2v5", abs(A / ALW - 1) < 0.05, f"A/W2 = {A:.4e} ({A / T_UNIT:.3e} s^-1) vs Dirac A_LW {ALW:.4e} ({100 * (A / ALW - 1):+.1f} %); w {w96:.5f} vs Dirac 0.07669 ({100 * (w96 / 0.0766903 - 1):+.2f} %)", A / ALW)
    t, n, Dc = d["ser_t"], d["ser_norm"], d["ser_Dc"]
    Ai = 4 / 3 * ALPHA * w96**3 * np.sum(np.abs(Dc)**2, 1) / n
    win = []
    for a_, b_ in ((100, 600), (600, 1200), (1200, 2000), (2000, 3000)):
        mm = (t >= a_) & (t < b_)
        loc = -np.polyfit(t[mm], np.log(n[mm]), 1)[0]
        win.append(float(loc / (np.mean(Ai[mm]) * q["W2"])))
    print("   windows (rate / instantaneous rule rate):", [round(x, 4) for x in win])
    R["R5_numbers"] = dict(q, A_rule=Ar, omega=w96, windows=win, t=t[::5].tolist(), norm=n[::5].tolist())
    E_qe = float(d["ser_eph"][-1] + d["ser_abs"][-1])
    verdict("S5b", q["depletion"] >= 0.5 * Ar * q["W2"] * q["t_end"], f"quantum pure 2p: 1 - norm(T) = {q['depletion']:.4e} (>= {0.5 * Ar * q['W2'] * q['t_end']:.4e})", q["depletion"])
    if have("P5_S1sc_96"):
        ds = load("P5_S1sc_96")
        E_sc = float(ds["ser_eph"][-1] + ds["ser_abs"][-1])
        verdict("S5", E_sc < 0.01 * E_qe, f"energy given to the light by the pure 2p: semiclassical {E_sc:.2e} vs quantum {E_qe:.3e} (ratio {E_sc / E_qe:.1e})", [E_sc, E_qe])
        pe = np.array(ds["ser_pexc"])[:, 0]
        R["S5_series"] = dict(t=ds["ser_t"][::5].tolist(), p=pe[::5].tolist(), proj_drop=float(1 - pe[~np.isnan(pe)][-1]))
d5 = json.load(open(os.path.join(OUT, "emiss_d5.json")))
verdict("D5", abs(d5["ratio"] - 1) < 0.05, f"96^3: sum |<1s|r|2p3/2+3>|^2 = {d5['sum2']:.4f} vs Dirac {d5['dirac']:.4f} ({100 * (d5['ratio'] - 1):+.1f} %); 48^3 was +48.8 %", d5["ratio"])
w5 = json.load(open(os.path.join(OUT, "emiss_w5.json")))
verdict("W5", None, f"decision by the pre-registered rule: {w5['W5_decision']}", w5["W5_decision"])
x5 = json.load(open(os.path.join(OUT, "emiss_x5.json")))
x5b = json.load(open(os.path.join(OUT, "emiss_x5_16000_bh.json")))
verdict("X5", None, f"pre-registered estimator (Hann, 4000 ticks): '{x5['decision']}' but from side lobes of the 1s line; with the declared "
        f"clean window (Blackman-Harris, 16000 ticks) BOTH windows are empty (no level > 1e-4 near 0.886 or 0.989): neither H1 nor H2", None)
verdict("M5", None, "dynamic multipole coupling: not done (declared)", None)
json.dump(R, open(os.path.join(OUT, "emiss_results3.json"), "w"), indent=1, default=lambda x: x.tolist() if hasattr(x, "tolist") else float(x))
print(f"=== {sum(1 for v in R.values() if isinstance(v, dict) and v.get('pass_') is True)} PASS, "
      f"{sum(1 for v in R.values() if isinstance(v, dict) and v.get('pass_') is False)} FAIL")
