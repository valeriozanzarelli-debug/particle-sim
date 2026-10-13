"""emiss_judge2.py - judges the phase-2 tests of lab/out/emiss_prereg2.json (thresholds copied, never changed here).
Uses the outputs of emiss_retard.py (N2), emiss_l1b.py (L1b), emiss_v4_checks.py (D4, V4) and the P4_* runs.
Writes lab/out/emiss_results2.json.  Run: python lab/emiss_judge2.py"""
import json, os, glob
import numpy as np
from zz_core import ALPHA, T_UNIT
from emiss_core import CK, OUT
from emiss_analysis import qe_rates, load
os.environ["EMISS_CFG"] = "none"
import emiss_configs as EC

L4 = json.load(open(os.path.join(OUT, "emiss_levels4_48.json")))
R = {}


def verdict(name, ok, msg, value=None):
    R[name] = dict(pass_=None if ok is None else bool(ok), msg=msg, value=value)
    print(f"[{'INFO' if ok is None else ('PASS' if ok else 'FAIL')}] {name}: {msg}", flush=True)


def have(name):
    return os.path.exists(os.path.join(CK, f"run_{name}.npz")) and any(
        f"=== {name} finished rc=0" in open(f).read() for f in glob.glob(os.path.join(CK, "drive*.log")))


def spectrum_peaks(t, y, lo=0.06, hi=0.10):
    h = np.sin(np.pi * (t - t[0]) / (t[-1] - t[0]))**2
    om = np.arange(lo, hi, 1e-5)
    S = np.abs(np.exp(1j * np.outer(om, t)) @ (h * y)) * 2 / np.sum(h)
    pk = [(om[i], S[i]) for i in range(1, len(om) - 1) if S[i] > S[i - 1] and S[i] > S[i + 1]]
    pk.sort(key=lambda p: -p[1])
    return pk, om, S


w1 = L4["1s+"]["omega"]
w32 = L4["2p3/2+"]["omega"] - w1
w12 = L4["2p1/2+"]["omega"] - w1
D = EC.atom4()
s1 = EC.level4("1s+")
for name, pol in (("P4_E1pi", "z"), ("P4_E1sig", "sigma-")):
    if not have(name):
        continue
    dr = EC.CONFIGS[name]["drive"]
    d = load(name)
    t = d["ser_t"]
    y = np.array(d["ser_probe"])[:, 2 if pol == "z" else 1]
    m = t >= dr["t0"] + 4 * dr["tau"] + 160.0
    pk, om, S = spectrum_peaks(t[m], y[m])
    Sw = lambda w: np.exp(-((w - dr["omega"]) * dr["tau"])**2 / 4)
    if pol == "z":
        d32 = D.dip_matrix(s1, EC.level4("2p3/2+"))[2]
        d12 = D.dip_matrix(s1, EC.level4("2p1/2+"))[2]
        ratio_p = abs(d32)**2 * Sw(w32) * w32**2 / (abs(d12)**2 * Sw(w12) * w12**2)
        verdict("E4a", abs(pk[0][0] / w32 - 1) < 0.005, f"strongest line {pk[0][0]:.5f} vs w(1s->2p3/2) {w32:.5f} ({100 * (pk[0][0] / w32 - 1):+.2f} %)", pk[0][0])
        verdict("E4b", abs(pk[1][0] / w12 - 1) < 0.005, f"second line {pk[1][0]:.5f} vs w(1s->2p1/2) {w12:.5f} ({100 * (pk[1][0] / w12 - 1):+.2f} %)", pk[1][0])
        rm = pk[0][1] / pk[1][1]
        verdict("E4c", abs(rm / ratio_p - 1) < 0.20, f"line ratio 3/2 : 1/2 = {rm:.3f} vs predicted {ratio_p:.3f} ({100 * (rm / ratio_p - 1):+.1f} %; |d|^2 {abs(d32)**2:.4f} vs {abs(d12)**2:.4f})", [rm, ratio_p])
        R["E4_spectrum_pi"] = dict(om=om[::20].tolist(), S=S[::20].tolist(), w32=w32, w12=w12, peaks=pk[:6])
    else:
        c32 = D.dip_matrix(EC.level4("2p3/2-"), s1)
        c12 = D.dip_matrix(EC.level4("2p1/2-"), s1)
        a32 = abs(c32[0] - 1j * c32[1]) * abs(D.dip_matrix(s1, EC.level4("2p3/2-"))[1]) * Sw(w32) * w32**2
        a12 = abs(c12[0] - 1j * c12[1]) * abs(D.dip_matrix(s1, EC.level4("2p1/2-"))[1]) * Sw(w12) * w12**2
        verdict("E4s", abs(pk[0][0] / w12 - 1) < 0.005, f"sigma-: strongest line {pk[0][0]:.5f} vs w(1s->2p1/2) {w12:.5f} ({100 * (pk[0][0] / w12 - 1):+.2f} %); "
                f"second {pk[1][0]:.5f}; ratio measured {pk[0][1] / pk[1][1]:.3f}, lattice prediction (1/2 : 3/2) {a12 / a32:.3f}", pk[0][0])
        R["E4_spectrum_sigma"] = dict(om=om[::20].tolist(), S=S[::20].tolist(), peaks=pk[:6])

w48 = L4["2p3/2+3"]["omega"] - w1
if have("P4_R1"):
    q = qe_rates("P4_R1", w48, 4.0)
    d = load("P4_R1")
    Ar = 4 / 3 * ALPHA * w48**3 * float(np.sum(np.abs(d["ser_Dc"][0])**2)) / d["ser_norm"][0]
    verdict("R4", abs(q["A_norm"] / (Ar * q["W2"]) - 1) < 0.03, f"A_norm {q['A_norm']:.4e} / (A_rule {Ar:.4e} x W2 {q['W2']:.4f}) = {q['A_norm'] / (Ar * q['W2']):.4f}", q["A_norm"] / (Ar * q["W2"]))
    verdict("R4b", abs(q["A_photons"] / q["A_norm"] - 1) < 0.03, f"A_photons / A_norm = {q['A_photons'] / q['A_norm']:.4f}", q["A_photons"] / q["A_norm"])
    t, n, Dc = d["ser_t"], d["ser_norm"], d["ser_Dc"]
    Ai = 4 / 3 * ALPHA * w48**3 * np.sum(np.abs(Dc)**2, 1) / n
    win = []
    for a_, b_ in ((100, 600), (600, 1200), (1200, 2000), (2000, 3000)):
        mm = (t >= a_) & (t < b_)
        loc = -np.polyfit(t[mm], np.log(n[mm]), 1)[0]
        win.append(dict(window=[a_, b_], ratio=float(loc / (np.mean(Ai[mm]) * q["W2"]))))
    R["R4_numbers"] = dict(q, A_rule=Ar, omega=w48, A_over_W2_s=q["A_norm"] / q["W2"] / T_UNIT, windows=win,
                           D_start=np.abs(Dc[0]).tolist(), D_end=np.abs(Dc[-1]).tolist(), t=t[::5].tolist(), norm=n[::5].tolist())
    print("   windows (rate / instantaneous rule rate):", [round(w_["ratio"], 4) for w_ in win])
    if have("P4_S1sc"):
        ds = load("P4_S1sc")
        pe = np.array(ds["ser_pexc"])[:, 0]
        pe = pe[~np.isnan(pe)]
        dsc = 1 - pe[-1]
        dep = q["depletion"]
        e_sc = float(ds["ser_eph"][-1] + ds["ser_abs"][-1])
        verdict("S4", dep >= 0.5 * Ar * q["W2"] * q["t_end"] and dsc < 0.01 * dep,
                f"qe 1 - norm(T) = {dep:.4e} (>= {0.5 * Ar * q['W2'] * q['t_end']:.4e}); semiclassical 1 - p_2p(T) = {dsc:.2e} ({dsc / dep:.2f} of it); "
                f"energy to light: semiclassical {e_sc:.1e}", [dep, dsc, e_sc])
        R["S4_series"] = dict(t=ds["ser_t"][::5].tolist(), p=np.array(ds["ser_pexc"])[::5, 0].tolist())

for f, keys in (("emiss_retard.json", None), ("emiss_l1b.json", None), ("emiss_v4_checks.json", None)):
    p = os.path.join(OUT, f)
    if not os.path.exists(p):
        continue
    j = json.load(open(p))
    if f == "emiss_retard.json":
        verdict("N2", j["N2"]["pass_"], f"(A_N1/W2) x F_ret(rule) = {j['N2']['value']:.5e} vs Dirac A_EX {j['N2']['dirac_A_EX']:.5e} ({100 * (j['N2']['ratio'] - 1):+.2f} %); "
                f"velocity k->0 / length on the rule's states {j['rule']['vel0_over_len']:.5f}", j["N2"]["ratio"])
        verdict("N2x", None, f"F_ret rule {j['N2x']['F_rule']:.5f} vs Dirac {j['N2x']['F_dirac']:.5f}", j["N2x"]["F_rule"])
    elif f == "emiss_l1b.json":
        verdict("L1b", j["pass_"], f"new export: amplitude at the old beat {j['new']['amplitude']:.4f} (< 0.01); old 1s free run {j['old_free_run']['amplitude']:.4f}", j["new"]["amplitude"])
    else:
        verdict("D4", None, "; ".join(f"{k}: {v['sum2']:.3f} ({v['ratio']:.2f} x Dirac)" for k, v in j["D4"].items()), j["D4"])
        verdict("V4", None, "; ".join(f"{k}: max drift {v['max']:.4f}" for k, v in j["V4"].items()), {k: v["max"] for k, v in j["V4"].items()})

json.dump(R, open(os.path.join(OUT, "emiss_results2.json"), "w"), indent=1, default=lambda x: x.tolist() if hasattr(x, "tolist") else float(x))
print(f"=== {sum(1 for v in R.values() if isinstance(v, dict) and v.get('pass_') is True)} PASS, "
      f"{sum(1 for v in R.values() if isinstance(v, dict) and v.get('pass_') is False)} FAIL")
