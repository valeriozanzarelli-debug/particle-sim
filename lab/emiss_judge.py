"""emiss_judge.py - judges the production runs against lab/out/emiss_prereg.json (thresholds copied from there,
never changed here).  Writes lab/out/emiss_results.json.  Run: python lab/emiss_judge.py"""
import json, os
import numpy as np
from zz_core import ALPHA, T_UNIT
from emiss_core import Dirac3D, Light3D, CK, OUT
from emiss_analysis import qe_rates, load

os.environ["EMISS_CFG"] = "none"
import emiss_configs as EC

EX = json.load(open(os.path.join(OUT, "emiss_exact.json")))
L48 = json.load(open(os.path.join(OUT, "emiss_levels_48.json")))
R = {}


def verdict(name, ok, msg, value=None):
    R[name] = dict(pass_=bool(ok) if ok is not None else None, msg=msg, value=value)
    print(f"[{'PASS' if ok else ('INFO' if ok is None else 'FAIL')}] {name}: {msg}", flush=True)


def have(name):
    """only runs that finished (the drive logs say '=== NAME finished rc=0')"""
    if not os.path.exists(os.path.join(CK, f"run_{name}.npz")):
        return False
    import glob
    return any(f"=== {name} finished rc=0" in open(f).read() for f in glob.glob(os.path.join(CK, "drive*.log")))


def spectrum_peaks(t, y, lo=0.06, hi=0.10):
    h = np.sin(np.pi * (t - t[0]) / (t[-1] - t[0]))**2
    om = np.arange(lo, hi, 1e-5)
    S = np.abs(np.exp(1j * np.outer(om, t)) @ (h * y)) * 2 / np.sum(h)
    pk = [(om[i], S[i]) for i in range(1, len(om) - 1) if S[i] > S[i - 1] and S[i] > S[i + 1]]
    pk.sort(key=lambda p: -p[1])
    return pk, om, S


# ------------------------------------------------------------------ E1v3
D48 = Dirac3D(48, 0.5)
w1 = L48["1s+"]["omega"]
w32 = L48["2p3/2+"]["omega"] - w1
w12 = L48["2p1/2+"]["omega"] - w1
for name, pol in (("P_E1pi", "z"), ("P_E1sig", "sigma-")):
    if not have(name):
        continue
    cfg = EC.CONFIGS[name]
    dr = cfg["drive"]
    d = load(name)
    t = d["ser_t"]
    pr = np.array(d["ser_probe"])
    comp = 2 if pol == "z" else 1
    y = pr[:, comp]
    t_start = dr["t0"] + 4 * dr["tau"] + 160.0
    m = t >= t_start
    pk, om, S = spectrum_peaks(t[m], y[m])
    s1 = EC.level(48, "1s+")
    Sw = lambda w: np.exp(-((w - dr["omega"]) * dr["tau"])**2 / 4)
    if pol == "z":
        d32 = D48.dip_matrix(s1, EC.level(48, "2p3/2+"))[2]
        d12 = D48.dip_matrix(s1, EC.level(48, "2p1/2+"))[2]
        amp32 = abs(d32)**2 * Sw(w32) * w32**2
        amp12 = abs(d12)**2 * Sw(w12) * w12**2
        verdict("E1v3a", abs(pk[0][0] / w32 - 1) < 0.005, f"strongest line {pk[0][0]:.5f} vs w(1s->2p3/2) {w32:.5f} ({100 * (pk[0][0] / w32 - 1):+.2f} %)", pk[0][0])
        verdict("E1v3b", abs(pk[1][0] / w12 - 1) < 0.005, f"second line {pk[1][0]:.5f} vs w(1s->2p1/2) {w12:.5f} ({100 * (pk[1][0] / w12 - 1):+.2f} %)", pk[1][0])
        ratio_m = pk[0][1] / pk[1][1]
        ratio_p = amp32 / amp12
        verdict("E1v3c", abs(ratio_m / ratio_p - 1) < 0.20, f"line ratio 3/2 : 1/2 = {ratio_m:.3f} vs predicted {ratio_p:.3f} "
                f"(|d|^2: {abs(d32)**2:.4f} vs {abs(d12)**2:.4f})", [ratio_m, ratio_p])
        # the v2 estimator
        zc = []
        mm = (t >= dr["t0"] + 4 * dr["tau"] + 400) & (t <= dr["t0"] + 4 * dr["tau"] + 1200)
        tt, yy = t[mm], y[mm]
        for i in range(len(yy) - 1):
            if (yy[i] <= 0 < yy[i + 1]) or (yy[i] >= 0 > yy[i + 1]):
                zc.append(tt[i] + (tt[i + 1] - tt[i]) * yy[i] / (yy[i] - yy[i + 1]))
        wzc = np.pi * (len(zc) - 1) / (zc[-1] - zc[0])
        verdict("E1v3zc", None, f"zero crossings {wzc:.5f}: vs 1s->2p3/2 {100 * (wzc / w32 - 1):+.2f} %, vs 1s->2p1/2 {100 * (wzc / w12 - 1):+.2f} %", wzc)
        R["E1_spectrum_pi"] = dict(om=om[::20].tolist(), S=S[::20].tolist(), w32=w32, w12=w12, peaks=pk[:5])
    else:
        sm = EC.level(48, "1s+")
        X = D48.x
        c32 = D48.dip_matrix(EC.level(48, "2p3/2-"), sm)
        c12 = D48.dip_matrix(EC.level(48, "2p1/2-"), sm)
        a32 = abs(c32[0] - 1j * c32[1]) * abs(D48.dip_matrix(sm, EC.level(48, "2p3/2-"))[1]) * Sw(w32) * w32**2
        a12 = abs(c12[0] - 1j * c12[1]) * abs(D48.dip_matrix(sm, EC.level(48, "2p1/2-"))[1]) * Sw(w12) * w12**2
        verdict("E1v3s", abs(pk[0][0] / w12 - 1) < 0.005, f"sigma-: strongest line {pk[0][0]:.5f} vs w(1s->2p1/2) {w12:.5f} "
                f"({100 * (pk[0][0] / w12 - 1):+.2f} %); second {pk[1][0]:.5f}; ratio 1/2 : 3/2 measured "
                f"{pk[0][1] / pk[1][1]:.3f}, predicted {a12 / a32:.3f}", pk[0][0])
        R["E1_spectrum_sigma"] = dict(om=om[::20].tolist(), S=S[::20].tolist(), peaks=pk[:5])

# ------------------------------------------------------------------ R0
if have("P_R0"):
    c = EC.CONFIGS["P_R0"]
    d = load("P_R0")
    t, E = d["ser_t"], d["ser_eph"] + d["ser_abs"]
    m = t >= 400
    P = np.polyfit(t[m], E[m], 1)[0]
    W2 = Light3D(128, 4.0, sponge=0, sigma=1.0).form_factor2(c["omega"])
    PL = c["I0"]**2 * c["omega"]**2 / (12 * np.pi)
    verdict("R0", abs(P / (PL * W2) - 1) < 0.03, f"P / (Larmor W2) = {P / (PL * W2):.4f} (Bloch prediction of the rule 1.006)", P / (PL * W2))


def rule_rate(name, omega):
    d = load(name)
    D0 = d["ser_Dc"][0]
    n0 = d["ser_norm"][0]
    return 4 / 3 * ALPHA * omega**3 * float(np.sum(np.abs(D0)**2)) / n0


# ------------------------------------------------------------------ R1 / S1 / K10 / SC10 (48^3)
w48 = L48["2p3/2+3"]["omega"] - w1
if have("P_R1_48"):
    q = qe_rates("P_R1_48", w48, 4.0)
    Ar = rule_rate("P_R1_48", w48)
    verdict("R1_48", abs(q["A_norm"] / (Ar * q["W2"]) - 1) < 0.03, f"A_norm {q['A_norm']:.4e} / (A_rule {Ar:.4e} x W2 {q['W2']:.4f}) = {q['A_norm'] / (Ar * q['W2']):.4f}", q["A_norm"] / (Ar * q["W2"]))
    verdict("R1b_48", abs(q["A_photons"] / q["A_norm"] - 1) < 0.03, f"A_photons / A_norm = {q['A_photons'] / q['A_norm']:.4f}", q["A_photons"] / q["A_norm"])
    R["R1_48_numbers"] = dict(q, A_rule=Ar, omega=w48, A_over_W2_s=q["A_norm"] / q["W2"] / T_UNIT)
    dep = q["depletion"]
    if have("P_S1sc"):
        ds = load("P_S1sc")
        pe = np.array(ds["ser_pexc"])[:, 0]
        pe = pe[~np.isnan(pe)]
        dsc = 1 - pe[-1]
        verdict("S1", dep >= 0.5 * Ar * q["W2"] * q["t_end"] and dsc < 0.01 * dep,
                f"pure 2p, quantum: 1 - norm(T) = {dep:.4e} (>= {0.5 * Ar * q['W2'] * q['t_end']:.4e}); semiclassical: 1 - p_2p(T) = {dsc:.2e} ({dsc / dep:.1e} of it)", [dep, dsc])
    if have("P_K10"):
        qk = qe_rates("P_K10", w48, 4.0)
        dk = load("P_K10")
        tk, nk = dk["ser_t"], dk["ser_norm"]
        dev = float(np.max(np.abs(nk - np.exp(-qk["A_norm"] * tk))))
        verdict("K10", dev < 0.02 and abs(qk["A_norm"] / 100 / q["A_norm"] - 1) < 0.05,
                f"max |norm - exp(-A t)| = {dev:.4f}; A(K=10)/100 / A(K=1) = {qk['A_norm'] / 100 / q['A_norm']:.4f}; norm(T) = {nk[-1]:.4f}", [dev, qk["A_norm"]])
        R["K10_series"] = dict(t=tk[::5].tolist(), norm=nk[::5].tolist(), A=qk["A_norm"])
        if have("P_SC10"):
            dsc = load("P_SC10")
            ts = dsc["ser_t"]
            ps = np.array(dsc["ser_pexc"])[:, 0]
            ok = ~np.isnan(ps)
            p0 = ps[0]
            A = qk["A_norm"]
            pl = p0 * np.exp(-A * ts) / (1 - p0 + p0 * np.exp(-A * ts))
            dev = float(np.max(np.abs(ps[ok] - pl[ok])))
            verdict("SC10", dev < 0.05, f"semiclassical p(t) vs logistic with A of K10 (no fit): max deviation {dev:.4f}; p(T) = {ps[ok][-1]:.4f} vs {pl[ok][-1]:.4f}", dev)
            R["SC10_series"] = dict(t=ts[ok][::5].tolist(), p=ps[ok][::5].tolist(), logistic=pl[ok][::5].tolist())

# ------------------------------------------------------------------ R1 / R2 96^3
if have("P_R1_96"):
    L96 = json.load(open(os.path.join(OUT, "emiss_levels_96.json")))
    w96 = L96["2p3/2+3"]["omega"] - L96["1s+"]["omega"]
    q = qe_rates("P_R1_96", w96, 4.0)
    Ar = rule_rate("P_R1_96", w96)
    verdict("R1_96", abs(q["A_norm"] / (Ar * q["W2"]) - 1) < 0.03, f"A_norm/(A_rule W2) = {q['A_norm'] / (Ar * q['W2']):.4f} (w = {w96:.5f})", q["A_norm"] / (Ar * q["W2"]))
    ALW = EX["Z60"]["2p3/2"]["A_LW"]
    verdict("R2_96", abs(q["A_norm"] / q["W2"] / ALW - 1) < 0.05, f"A/W2 = {q['A_norm'] / q['W2']:.4e} vs Dirac A_LW {ALW:.4e} ({100 * (q['A_norm'] / q['W2'] / ALW - 1):+.1f} %)", q["A_norm"] / q["W2"] / ALW)
    R["R1_96_numbers"] = dict(q, A_rule=Ar, omega=w96)

# ------------------------------------------------------------------ radial
for name, Z in (("P_N1", 60), ("P_H1a", 1), ("P_H1b", 1)):
    if not have(name):
        continue
    os.environ["EMISS_CFG"] = name
    import importlib
    importlib.reload(EC)
    c = EC.CONFIGS[name]
    q = qe_rates(name, c["omega"], c["b"])
    Ar = rule_rate(name, c["omega"])
    R[name + "_numbers"] = dict(q, A_rule=Ar, omega=c["omega"], b=c["b"], A_over_W2=q["A_norm"] / q["W2"])
    verdict(name + "_R1", abs(q["A_norm"] / (Ar * q["W2"]) - 1) < 0.03, f"A_norm/(A_rule W2) = {q['A_norm'] / (Ar * q['W2']):.4f}; photons/norm {q['A_photons'] / q['A_norm']:.4f}", q["A_norm"] / (Ar * q["W2"]))
    if name == "P_N1":
        A = q["A_norm"] / q["W2"]
        ALW, AEX = EX["Z60"]["2p3/2"]["A_LW"], EX["Z60"]["2p3/2"]["A_EX"]
        verdict("N1", abs(A / ALW - 1) < 0.03, f"A/W2 = {A:.5e} ({A / T_UNIT:.4e} s^-1) vs Dirac A_LW {ALW:.5e} ({100 * (A / ALW - 1):+.2f} %); Z^4 NR {EX['Z60']['A_NR_Z4']:.4e} ({100 * (A / EX['Z60']['A_NR_Z4'] - 1):+.1f} %)", A / ALW)
        verdict("N1x", None, f"vs Dirac A_EX (retardation) {AEX:.5e}: {100 * (A / AEX - 1):+.2f} %", A / AEX)
if "P_H1a_numbers" in R and "P_H1b_numbers" in R:
    Aa, Ab = R["P_H1a_numbers"]["A_over_W2"], R["P_H1b_numbers"]["A_over_W2"]
    A0 = (4 * Ab - Aa) / 3
    ref = EX["Z1"]["2p3/2"]["A_EX"]
    verdict("H1", abs(A0 / ref - 1) < 0.03, f"hydrogen: A(a=0.5) {Aa / T_UNIT:.4e}, A(a=0.25) {Ab / T_UNIT:.4e}, Richardson {A0 / T_UNIT:.4e} s^-1 vs Dirac "
            f"{ref / T_UNIT:.4e} ({100 * (A0 / ref - 1):+.2f} %); x mu/m = {A0 / T_UNIT * (1 - 1 / 1837.15267343):.4e} vs NIST 6.2649e8", A0 / ref)

# ------------------------------------------------------------------ toy and Laboratory replay
f = os.path.join(OUT, "emiss_twa.json")
if os.path.exists(f):
    tw = json.load(open(f))
    verdict("T1", abs(tw["T1_htwa_initial_rate_over_A"] - 0.5) < 0.1, f"hybrid TWA initial rate / A = {tw['T1_htwa_initial_rate_over_A']:.3f} (predicted 0.5)", tw["T1_htwa_initial_rate_over_A"])
    verdict("T2", abs(tw["T2_htwa_ground_p_end"] - 0.293) < 0.05, f"hybrid TWA from the ground state: p(3000) = {tw['T2_htwa_ground_p_end']:.3f} (predicted 0.293)", tw["T2_htwa_ground_p_end"])
    verdict("T3", abs(tw["T3_stwa_rate_over_A"] - 1) < 0.15, f"spin TWA rate / A = {tw['T3_stwa_rate_over_A']:.3f}", tw["T3_stwa_rate_over_A"])
    verdict("T4", tw["T4_sc_pure_e_dev"] < 1e-6, f"semiclassical pure e: max |p - 1| = {tw['T4_sc_pure_e_dev']:.1e}", tw["T4_sc_pure_e_dev"])
f = os.path.join(OUT, "emiss_lab1s.json")
if os.path.exists(f):
    lb = json.load(open(f))
    pp, pl = lb["free"]["period_pred"], lb["lab_series"]["beat_period"]
    verdict("L1", abs(pp / pl - 1) < 0.05, f"beat period predicted {pp:.1f} vs Laboratory series {pl:.1f} ({100 * (pp / pl - 1):+.2f} %)", pp / pl)

json.dump(R, open(os.path.join(OUT, "emiss_results.json"), "w"), indent=1, default=lambda x: x.tolist() if hasattr(x, "tolist") else float(x))
npass = sum(1 for v in R.values() if isinstance(v, dict) and v.get("pass_") is True)
nfail = sum(1 for v in R.values() if isinstance(v, dict) and v.get("pass_") is False)
print(f"=== {npass} PASS, {nfail} FAIL")
