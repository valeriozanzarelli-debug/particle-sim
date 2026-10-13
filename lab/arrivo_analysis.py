"""arrivo_analysis.py - judges the exams E1-E18 of lab/out/arrivo_prereg.json from the production runs.
Reads lab/out/arrivo_r1_*.json, lab/out/arrivo_r3.json and the scratch npz of r1, r2, r4.
Writes lab/out/arrivo_results.json. Amendments (declared): lab/out/arrivo_emendamenti.json."""
import json
import os
import time
import numpy as np
import arrivo_nr as A
import arrivo_rule1d as R1

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SCR = os.environ.get("ARRIVO_SCRATCH", "/tmp/claude-0/-home-claude-particle-sim/9de435bc-26c2-5d88-b70d-07efb8033e5b/scratchpad")
t0 = time.time()
R = {}
EX = {}


def ld(p):
    return np.load(os.path.join(SCR, p), allow_pickle=True)


def rich(f_a, f_a2):
    """Richardson a^2 from f(a), f(a/2)"""
    return (4 * f_a2 - f_a) / 3


def quadfit(avals, fvals):
    p = np.polyfit(np.asarray(avals), np.asarray(fvals), 2)
    return float(p[-1])


def st(T):
    f = np.isfinite(T)
    Tf = T[f]
    return dict(arrived=float(f.mean()), mean=float(Tf.mean()), q50=float(np.median(Tf)),
                q90=float(np.percentile(Tf, 90)), q999=float(np.percentile(Tf, 99.9)), max=float(Tf.max()))


def verdict(ok):
    return "PASS" if ok else "FAIL"


def refine_max(name, L, om, te, T, X0, Y0, Z0, ntop=50):
    """amendment A1: the latest ntop (and all not-arrived) NR trajectories re-integrated with dt/4"""
    Tf = np.where(np.isfinite(T), T, np.inf)
    idx = np.argsort(-Tf)[:ntop]
    f = 0.25
    Tn = A.bohm_arrivals(len(idx), L, om, np.pi / 2, te, None,
                         dt_fn=lambda t: f * (2e-4 if t < 0.05 else (2e-3 if t < 2 else 1e-2)),
                         XYZ0=(X0[idx], Y0[idx], Z0[idx]))
    T2 = T.copy()
    T2[idx] = Tn
    return T2


# ------------------------------------------------------------------ NR reference
nr = {}
for name, L, om, te in (("NR_A", 100.0, 1000.0, 60.0), ("NR_B", 10.0, 1000.0, 40.0), ("NR_C", 10.0, 20.0, 40.0)):
    js = json.load(open(os.path.join(OUT, "arrivo_r1_%s.json" % name)))
    d = ld("arrivo_r1_%s.npz" % name)
    Tx_raw = d["T_alpha=1.5708"]
    Tx = refine_max(name, L, om, te, Tx_raw, d["X0"], d["Y0"], d["Z0"])
    nr[name] = dict(Tx=Tx, Tz=d["T_alpha=0.0000"], js=js, X0=d["X0"], Y0=d["Y0"], Z0=d["Z0"], d=d)
    R[name] = dict(x_raw=st(Tx_raw), x_refined=st(Tx), z=st(d["T_alpha=0.0000"]), KS_up_flux=js["KS_spin_up_vs_flux"],
                   dt_half=js["dt_half_check"])
    if name == "NR_C":
        R[name]["alpha_scan_mean"] = {k: v["mean"] for k, v in js["runs"].items()}
    print(name, R[name]["x_refined"], flush=True)

taumax_A = R["NR_A"]["x_refined"]["max"]
EX["E1"] = dict(value=taumax_A, raw_value=R["NR_A"]["x_raw"]["max"], ref=42.9, rel=taumax_A / 42.9 - 1,
                verdict=verdict(abs(taumax_A / 42.9 - 1) < 0.02),
                note="raw sample max 43.69 passes only because of 2 badly integrated trajectories (amendment A1)")
EX["E2"] = dict(value=R["NR_A"]["KS_up_flux"], verdict=verdict(R["NR_A"]["KS_up_flux"] < 0.015))
Tz = nr["NR_A"]["Tz"]
p_after = float(np.mean(~np.isfinite(Tz) | (Tz > taumax_A)))
EX["E3"] = dict(value=p_after, verdict=verdict(p_after > 0.05))

# ------------------------------------------------------------------ rule 1-D: flux
tg = np.linspace(1e-4, 60, 600001)


def L1_vs_qf(L, tau_c, f_c, t_lo, t_hi):
    sel = (tau_c >= t_lo) & (tau_c <= t_hi)
    ex = A.flux(L, tau_c[sel])
    dt = np.gradient(tau_c[sel])
    return float(np.sum(np.abs(f_c[sel] - ex) * dt))


def flux_curve(tag):
    d = ld("r2/%s.npz" % tag)
    return d["tau"], d["flux"]


def common(t1, f1, t2, f2, tmax):
    tt = np.linspace(0.02, tmax, 20000)
    return tt, np.interp(tt, t1, f1), np.interp(tt, t2, f2)


fl = {}
for D in (20, 40, 80):
    tmax = 10.0
    curves = {}
    for a in (0.5, 0.25, 0.125, 0.0625):
        try:
            curves[a] = flux_curve("z_L10_D%d_a%g" % (D, a))
        except FileNotFoundError:
            pass
    rows = {}
    for a, (t, f) in curves.items():
        rows[str(a)] = L1_vs_qf(10.0, t, f, 0.02, tmax)
    if 0.125 in curves and 0.0625 in curves:
        tt, f1, f2 = common(*curves[0.125], *curves[0.0625], tmax)
        f0 = rich(f1, f2)
        rows["extrap"] = L1_vs_qf(10.0, tt, f0, 0.02, tmax)
    fl[D] = rows
R["flux_L10_L1"] = fl
EX["E4"] = dict(value=fl[80].get("extrap"), by_D={D: fl[D].get("extrap") for D in fl},
                verdict=verdict(fl[80].get("extrap", 9) < 0.01 and fl[20]["extrap"] > fl[40]["extrap"] > fl[80]["extrap"]))

# L = 100 flux from the x runs (flux is spin independent)
c100 = {}
for D, alist in ((20, (0.5, 0.25, 0.125, 0.0625)), (40, (0.5, 0.25, 0.125))):
    rows = {}
    cur = {}
    for a in alist:
        try:
            cur[a] = flux_curve("x_L100_D%d_a%g" % (D, a))
            rows[str(a)] = L1_vs_qf(100.0, *cur[a], 0.05, 50.0)
        except FileNotFoundError:
            pass
    pairs = [(0.125, 0.0625)] if D == 20 else [(0.25, 0.125)]
    for (p, q) in pairs:
        if p in cur and q in cur:
            tt, f1, f2 = common(*cur[p], *cur[q], 50.0)
            rows["extrap"] = L1_vs_qf(100.0, tt, rich(f1, f2), 0.05, 50.0)
    c100[D] = rows
R["flux_L100_L1"] = c100
EX["E5"] = dict(value=c100[40].get("extrap"), D20=c100[20].get("extrap"),
                verdict=verdict(c100[40].get("extrap", 9) < 0.01))

# tail ratio
tails = {}
for D in (20, 40):
    vals = {}
    cur = {}
    for a in (0.5, 0.25, 0.125, 0.0625):
        try:
            t, f = flux_curve("z_L10_D%d_a%g" % (D, a))
        except FileNotFoundError:
            continue
        sel = (t >= 20) & (t <= 40)
        vals[str(a)] = float(np.mean(f[sel] / A.flux(10.0, t[sel])))
        cur[a] = vals[str(a)]
    if 0.125 in cur and 0.0625 in cur:
        vals["extrap"] = rich(cur[0.125], cur[0.0625])
    vals["Mstar_cubed_a0.0625"] = (np.tan(0.0625) / 0.0625)**3
    tails[D] = vals
R["tail_ratio_L10"] = tails
EX["E6"] = dict(value={D: tails[D].get("extrap") for D in tails},
                verdict=verdict(all(abs(tails[D].get("extrap", 9) - 1) < 0.01 for D in tails)))
# edge check (ramp vs cone)
try:
    tc, fc = flux_curve("z_L10_D20_a0.25")
    trp, frp = flux_curve("zramp_L10_D20_a0.25")
    tt, f1, f2 = common(tc, fc, trp, frp, 40.0)
    R["edge_check_L1_ramp_vs_cone"] = float(np.sum(np.abs(f1 - f2)) * (tt[1] - tt[0]))
    R["edge_check_T"] = dict(cone=st(ld("r2/z_L10_D20_a0.25.npz")["T"]), ramp=st(ld("r2/zramp_L10_D20_a0.25.npz")["T"]))
except FileNotFoundError:
    pass

# ------------------------------------------------------------------ hybrid surfers
def hyb(tag):
    return ld("r2/%s.npz" % tag)["T"]


hy = {}
for D in (20, 40):
    rows = {}
    for a in (0.5, 0.25, 0.125, 0.0625):
        try:
            rows[str(a)] = st(hyb("x_L10_D%d_a%g" % (D, a)))
        except FileNotFoundError:
            pass
    if "0.125" in rows and "0.0625" in rows:
        rows["extrap"] = {k: rich(rows["0.125"][k], rows["0.0625"][k]) for k in ("mean", "q50", "q90", "max", "q999")}
        rows["quadfit"] = {k: quadfit([0.25, 0.125, 0.0625], [rows[s][k] for s in ("0.25", "0.125", "0.0625")])
                           for k in ("mean", "q50", "q90", "max", "q999")}
    hy[D] = rows
ref = R["NR_B"]["x_refined"]
R["hybrid_L10"] = hy
rels = {k: hy[40]["extrap"][k] / ref[k] - 1 for k in ("mean", "q50", "q90", "max")}
EX["E7"] = dict(rel=rels, rel_q999=hy[40]["extrap"]["q999"] / ref["q999"] - 1,
                verdict=verdict(all(abs(v) < 0.015 for v in rels.values())))

hy100 = {}
nA = 5000
refA = st(nr["NR_A"]["Tx"][:nA])
for D, alist in ((20, (0.5, 0.25, 0.125, 0.0625)), (40, (0.5, 0.25, 0.125))):
    rows = {}
    for a in alist:
        try:
            rows[str(a)] = st(hyb("x_L100_D%d_a%g" % (D, a)))
        except FileNotFoundError:
            pass
    p, q = ("0.125", "0.0625") if D == 20 else ("0.25", "0.125")
    if p in rows and q in rows:
        rows["extrap"] = {k: rich(rows[p][k], rows[q][k]) for k in ("mean", "q50", "q90", "max", "q999")}
    hy100[D] = rows
R["hybrid_L100"] = hy100
R["NR_A_first5000"] = refA
rels = {k: hy100[40]["extrap"][k] / refA[k] - 1 for k in ("mean", "max")}
EX["E8"] = dict(rel=rels, rel_D20={k: hy100[20]["extrap"][k] / refA[k] - 1 for k in ("mean", "max")},
                verdict=verdict(all(abs(v) < 0.02 for v in rels.values())))


def ks2(T1, T2, tmax):
    """two-sample KS on [0, tmax] with the full sample sizes as denominators (not-arrived count as > tmax)"""
    g = np.linspace(0, tmax, 4001)
    c1 = np.array([np.sum(T1 <= x) for x in g]) / T1.size
    c2 = np.array([np.sum(T2 <= x) for x in g]) / T2.size
    return float(np.max(np.abs(c1 - c2)))


Tz_rule = hyb("z_L10_D40_a0.0625")
Tz_nr = nr["NR_B"]["Tz"][:Tz_rule.size]
EX["E9"] = dict(value=ks2(np.nan_to_num(Tz_rule, nan=1e9), np.nan_to_num(Tz_nr, nan=1e9), 40.0),
                stats_rule=st(Tz_rule), stats_nr=st(Tz_nr))
EX["E9"]["verdict"] = verdict(EX["E9"]["value"] < 0.02)

# ------------------------------------------------------------------ 2-D rule surfers
two = {}
for s in ("x", "y", "z"):
    for a in (0.5, 0.25, 0.125):
        try:
            d = ld("r4/free_%s_a%g.npz" % (s, a))
            two["%s_%g" % (s, a)] = dict(stats=st(d["T"]), spin0=d["spin0"].tolist(), prep_resid=float(d["prep_resid"]))
        except FileNotFoundError:
            pass
R["rule2d_free"] = two
refC = st(nr["NR_C"]["Tx"][:2000])
refCz = nr["NR_C"]["Tz"][:2000]
R["NR_C_first2000"] = dict(x=refC, z=st(refCz))
if "x_0.25" in two and "x_0.125" in two:
    ex2 = {k: rich(two["x_0.25"]["stats"][k], two["x_0.125"]["stats"][k]) for k in ("mean", "max", "q50", "q90", "q999")}
    R["rule2d_x_extrap"] = ex2
    rels = {k: ex2[k] / refC[k] - 1 for k in ("mean", "max")}
    EX["E10"] = dict(rel=rels, verdict=verdict(all(abs(v) < 0.03 for v in rels.values())))
ks = {}
for s in ("z", "y"):
    try:
        T = ld("r4/free_%s_a0.125.npz" % s)["T"]
        ks[s] = ks2(np.nan_to_num(T, nan=1e9), np.nan_to_num(refCz, nan=1e9), 8.0)
    except FileNotFoundError:
        pass
if ks:
    EX["E11"] = dict(value=ks, verdict=verdict(all(v < 0.03 for v in ks.values()) and len(ks) == 2))
try:
    Tx2 = ld("r4/free_x_a0.125.npz")["T"]
    Tz2 = ld("r4/free_z_a0.125.npz")["T"]
    tmx = np.nanmax(Tx2)
    frac_z = float(np.mean(~np.isfinite(Tz2) | (Tz2 > tmx)))
    n_x_late = int(np.sum(~np.isfinite(Tx2) | (Tx2 > 1.03 * refC["max"])))
    EX["E12"] = dict(frac_z_after=frac_z, taumax_x_rule=float(tmx), n_x_late=n_x_late,
                     verdict=verdict(frac_z > 0.03 and n_x_late == 0))
except FileNotFoundError:
    pass

# ------------------------------------------------------------------ detectors
r3 = json.load(open(os.path.join(OUT, "arrivo_r3.json")))
e13 = r3["E13"]
EX["E13"] = dict(rows=e13, verdict=verdict(all(abs(r["R_lattice"] / r["R_dirac"] - 1) < 1e-3 for r in e13)),
                 rel_vs_dirac=[abs(r["R_lattice"] / r["R_dirac"] - 1) for r in e13])
det = {}
fF = r3["flux_L10"]
for D in (20, 40, 80):
    cur = {}
    for a in (0.5, 0.25, 0.125, 0.0625):
        try:
            d = ld("r2/det_L10_D%d_a%g.npz" % (D, a))
            cur[a] = (d["tau"], d["F"])
        except FileNotFoundError:
            pass
    tmax = 100.0 if D < 80 else 60.0
    tt = np.linspace(0.05, tmax, 5000)
    abr = r3["ABR_L10"]["2D_%d" % D]
    Fabr = np.interp(tt, abr["tau"], abr["F"])
    Ffl = np.interp(tt, fF["tau"], fF["F"])
    rows = {}
    for a, (t, F) in cur.items():
        rows[str(a)] = float(np.max(np.abs(np.interp(tt, t, F) - Fabr)))
    if 0.125 in cur and 0.0625 in cur:
        F0 = rich(np.interp(tt, *cur[0.125]), np.interp(tt, *cur[0.0625]))
        rows["extrap"] = float(np.max(np.abs(F0 - Fabr)))
        rows["vs_flux_extrap"] = float(np.max(np.abs(F0 - Ffl)))
        rows["F20_extrap"] = float(np.interp(20, tt, F0))
        rows["F20_abr"] = float(np.interp(20, tt, Fabr))
        rows["F20_flux"] = float(np.interp(20, tt, Ffl))
    det[D] = rows
R["det1d_vs_abr"] = det
EX["E14"] = dict(value={D: det[D].get("extrap") for D in det},
                 verdict=verdict(det[80].get("extrap", 9) < 0.01 and det[20]["extrap"] > det[40]["extrap"] > det[80]["extrap"]))
EX["E15"] = dict(value=det[80].get("vs_flux_extrap"), verdict=verdict(det[80].get("vs_flux_extrap", 0) > 0.1))

# 2-D detector: spin independence and comparison with the 1-D detector at the same D, a
d2 = {}
for a in (0.5, 0.25, 0.125):
    Fs = {}
    for s in ("x", "y", "z"):
        try:
            d = ld("r4/det_%s_a%g.npz" % (s, a))
            Fs[s] = (d["tau"], d["F"])
        except FileNotFoundError:
            pass
    if "x" in Fs and "z" in Fs:
        row = dict(xz=float(np.max(np.abs(Fs["x"][1] - Fs["z"][1]))))
        if "y" in Fs:
            row["yz"] = float(np.max(np.abs(Fs["y"][1] - Fs["z"][1])))
        # 1-D rule detector at D = 30, same a, same tau_end
        tau_end = Fs["z"][0][-1]
        run = R1.Run1D(30, a, 10 * 30, "detector", tau_end * 900)
        run.advance(Fs["z"][1].size)
        F1 = np.cumsum(np.concatenate(run.absorbed))
        row["vs1d"] = float(np.max(np.abs(F1 - Fs["z"][1])))
        row["F_end_2d"] = float(Fs["z"][1][-1])
        row["F_end_1d"] = float(F1[-1])
        d2[str(a)] = row
R["det2d"] = d2
if d2:
    EX["E16"] = dict(value=d2, verdict=verdict(all(max(r["xz"], r.get("yz", 0)) < 1e-3 for r in d2.values())))
    EX["E17"] = dict(value={k: r["vs1d"] for k, r in d2.items()}, verdict=verdict(all(r["vs1d"] < 0.02 for r in d2.values())))
F20pi = r3["ABR_L10"]["pi"]["F20"]
EX["E18"] = dict(value=F20pi, ref=0.96, verdict=verdict(abs(F20pi - 0.96) < 0.02))

R["exams"] = EX
R["summary"] = {k: v["verdict"] for k, v in EX.items()}
R["secs"] = time.time() - t0


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


json.dump(clean(R), open(os.path.join(OUT, "arrivo_results.json"), "w"), indent=1)
print(json.dumps(clean(R["summary"]), indent=0))
for k, v in EX.items():
    print(k, {kk: vv for kk, vv in clean(v).items() if kk != "rows"})
