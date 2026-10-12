"""press_ode.py - P1, P2, P3, P5 of the pressure study (see lab/out/press_prereg.json).
Energy-momentum tensor of the exact B=1 hedgehog (knot_core.radial_bvp), Skyrme units:
  e2 = f'^2 + 2 S/r^2, e4 = 2 S f'^2/r^2 + S^2/r^4  (S = sin^2 f),  T00 = e2 + e4
  p = (e4 - e2)/3,  s = 2 (f'^2 - S/r^2)(1 + S/r^2),  T^mu_mu = 2 e2
Writes lab/out/press_ode.json and lab/out/press_profiles.npz (for the figure)."""
import sys, os, json, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot_core as kc

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
PRE = json.load(open(os.path.join(OUT, "press_prereg.json")))
HBARC = 0.1973269804            # GeV fm
GEVFM3_PA = 1.602176634e-10 / 1e-45   # 1 GeV/fm^3 in Pa

# ------------------------------------------------------------------ profile (inner solve + analytic 1/r^2 tail)
R = 60.0
P = kc.radial_bvp(R)
C = P["C_tail"]
n_in = len(P["fp"])
r_in = P["r"][2:1 + n_in]; f_in = P["f"][2:1 + n_in]; fp_in = P["fp"][1:]   # drop r=1e-4: collocation BC node has a bad f'
rt = np.geomspace(R, 2e4, 6000)[1:]
r = np.concatenate([r_in, rt]); f = np.concatenate([f_in, C / rt**2]); fp = np.concatenate([fp_in, -2 * C / rt**3])
S = np.sin(f) ** 2
e2 = fp**2 + 2 * S / r**2
e4 = 2 * S * fp**2 / r**2 + S**2 / r**4
T00 = e2 + e4
p = (e4 - e2) / 3
p2, p4 = -e2 / 3, e4 / 3
s = 2 * (fp**2 - S / r**2) * (1 + S / r**2)
trace = 2 * e2
a = -float(np.polyfit(r_in[3:60]**2, fp_in[3:60], 1)[1])    # f = pi - a r + O(r^3) near 0
I = lambda g: np.trapezoid(g, r)
Rmax = r[-1]
# analytic tails beyond Rmax: T00 ~ 6C^2/r^6, p ~ -2C^2/r^6, s ~ 6C^2/r^6
E2 = 4 * np.pi * (I(r**2 * e2) + 2 * C**2 / Rmax**3)
E4 = 4 * np.pi * I(r**2 * e4)
M = E2 + E4
Ip0 = I(r**2 * p) - (2 / 3) * C**2 / Rmax**3
Ipabs = I(r**2 * np.abs(p))
Ip4 = 4 * np.pi * (I(r**4 * p) - 2 * C**2 / Rmax)        # int d^3r r^2 p
Is4 = 4 * np.pi * (I(r**4 * s) + 6 * C**2 / Rmax)        # int d^3r r^2 s
Itr = 4 * np.pi * (I(r**2 * trace) + 4 * C**2 / Rmax**3)
Itol = 4 * np.pi * (I(r**2 * (T00 + 3 * p)))             # tail of T00+3p is O(r^-10)

# sign structure
sg = np.sign(p[r < 50])
changes = int(np.sum(sg[1:] != sg[:-1]))
i0 = int(np.argmax(sg < 0))
r0 = float(np.interp(0.0, [p[i0], p[i0 - 1]], [r[i0], r[i0 - 1]]))
tail_neg = bool(np.all(p[(r > r0 + 1e-9) & (r < 50)] < 0))
p0 = a**2 * (a**2 - 1); T00_0 = 3 * a**2 + 3 * a**4
r2p = r**2 * p
ipk = int(np.argmax(r2p)); imn = int(np.argmin(r2p[r < 20]))


def D_of(e):
    return M * Ip4 / (4 * e**4), -(4 / 15) * M * Is4 / (4 * e**4)


def phys(e, F):            # F in MeV
    eps0 = F / (4 * e) / 1000.0                      # GeV
    ell = 2 * HBARC * 1000.0 / (e * F)               # fm
    return eps0, ell, eps0 / ell**3                  # GeV, fm, GeV/fm^3


res = {"profile": {"R_bvp": R, "C_tail": C, "a_slope": a, "bvp_status": P["bvp_status"]}}
# ---------------------------------------------------------------- P1
Dp, Ds = D_of(1.0)          # e = 1 -> pure number M I / 4 ; scale 1/e^4
P1 = {"M": M, "E_bvp": P["E"], "a_rel": abs(M - P["E"]) / P["E"],
      "vonLaue_rel": abs(Ip0) / Ipabs, "p0": p0, "T00_0": T00_0, "sign_changes_r<50": changes, "r0": r0,
      "tail_negative": tail_neg, "D_times_e4_from_p": Dp, "D_times_e4_from_s": Ds, "D_rel_diff": abs(Dp - Ds) / abs(Dp),
      "r2p_max_at": float(r[ipk]), "r2p_min_at": float(r[imn])}
P1["pass_a"] = bool(P1["a_rel"] < 1e-6)
P1["pass_b"] = bool(P1["vonLaue_rel"] < 1e-3)
P1["pass_c"] = bool(bool(p0 > 0 and changes == 1 and tail_neg))
P1["pass_d"] = bool(P1["D_rel_diff"] < 0.02)
# P1e Cebulla chiral limit
cb = PRE["references"]["Cebulla2007"]; ec, Fc = cb["params"]["e"], cb["params"]["F_pi_MeV"]
eps0, ell, pu = phys(ec, Fc)
mine = {"T00_0_GeV_fm3": T00_0 * pu, "p0_GeV_fm3": p0 * pu, "r0_fm": r0 * ell, "d1": 1.25 * Dp / ec**4,
        "M_MeV": M * eps0 * 1000, "int_r2p4_MeV": E4 / (12 * np.pi) * eps0 * 1000}
ref = cb["chiral_limit_table2"]
dev = {k: mine[k] / ref[k] - 1 for k in ref}
P1["e_cebulla"] = {"mine": mine, "ref": ref, "rel_dev": dev}
P1["pass_e"] = bool(bool(all(abs(v) < 0.02 for v in dev.values())))
res["P1"] = P1

# ---------------------------------------------------------------- P2
e, F = PRE["references"]["ANW"]["e"], PRE["references"]["ANW"]["F_pi_MeV"]
eps0, ell, pu = phys(e, F)
D_anw = Dp / e**4
beg = dict(d1=-2.04, M2=1.02, alpha=3, mN=0.938272)
D0 = 0.8 * beg["d1"]; Mt = math.sqrt(beg["M2"])
pbeg0 = -D0 * Mt**5 / (64 * np.pi * beg["mN"]) / HBARC**3           # GeV/fm^3
r0beg = 3 / Mt * HBARC
P2 = {"eps0_MeV": eps0 * 1000, "ell_fm": ell, "M_sol_MeV": M * eps0 * 1000,
      "p0_GeV_fm3": p0 * pu, "p0_Pa": p0 * pu * GEVFM3_PA, "T00_0_GeV_fm3": T00_0 * pu, "r0_fm": r0 * ell,
      "r2p_max_at_fm": float(r[ipk]) * ell, "r2p_min_at_fm": float(r[imn]) * ell, "D": D_anw, "d1": 1.25 * D_anw,
      "BEG": {"D0": D0, "p0_GeV_fm3": pbeg0, "p0_Pa": pbeg0 * GEVFM3_PA, "r0_fm": r0beg,
              "r2p_max_at_fm": (3 - math.sqrt(3)) / Mt * HBARC, "r2p_min_at_fm": (3 + math.sqrt(3)) / Mt * HBARC}}
P2["ratio_r0"] = P2["r0_fm"] / 0.6
P2["ratio_p0_vs_1e35"] = P2["p0_Pa"] / 1e35
P2["ratio_p0_vs_BEG"] = P2["p0_Pa"] / P2["BEG"]["p0_Pa"]
P2["pass_a"] = bool(abs(P2["ratio_r0"] - 1) < 0.25)
P2["pass_b1"] = bool(0.5 < P2["ratio_p0_vs_1e35"] < 2)
P2["pass_b2"] = bool(0.5 < P2["ratio_p0_vs_BEG"] < 2)
P2["pass_c"] = bool(-5.28 <= D_anw <= -3.58)
P2["D_compare"] = {"BEG_quark_mu2_1.4": D0, "lattice_total_2024": [-3.87, 0.97],
                   "lattice_pull_sigma": (D_anw + 3.87) / 0.97}
res["P2"] = P2

# ---------------------------------------------------------------- P3
P3 = {"E2": E2, "E4": E4, "virial_rel": abs(E2 - E4) / M, "trace_int": Itr, "trace_rel": abs(Itr - M) / M,
      "tolman_int": Itol, "tolman_rel": abs(Itol - M) / M, "frac_E2": E2 / M, "frac_E4": E4 / M,
      "skyrme_trace_pointwise_max": float(np.max(np.abs(e4 - 3 * p4))),
      "int_r2_p4_over_M": I(r**2 * p4) / (M / (4 * np.pi)), "MeV_ANW": {"E2": E2 * eps0 * 1e3, "E4": E4 * eps0 * 1e3}}
P3["pass_a"] = bool(P3["virial_rel"] < 1e-4)
P3["pass_b"] = bool(P3["trace_rel"] < 1e-4 and P3["skyrme_trace_pointwise_max"] < 1e-12)
P3["pass_c"] = bool(P3["tolman_rel"] < 1e-4)
res["P3"] = P3

# ---------------------------------------------------------------- P5
G, c = 6.67430e-11, 2.99792458e8
Mkg = 939e6 * 1.602176634e-19 / c**2
ratio = G * Mkg / (0.6e-15 * c**2)
rad = 1.0 + 3 * (1 / 3)       # rho + 3p/c^2 for p = rho c^2/3
res["P5"] = {"GM2_over_R_over_Mc2": ratio, "pass_a": 1e-40 <= ratio <= 1e-38, "tolman_radiation_factor": rad,
             "pass_b": rad == 2.0, "M_kg": Mkg,
             "note": "Skyrme term alone has p4 = e4/3 (radiation-like, Tolman weight 2), sigma term p2 = -e2/3 (weight 0); "
                     "summed over the knot von Laue gives total Tolman mass = M exactly"}

json.dump(res, open(os.path.join(OUT, "press_ode.json"), "w"), indent=1, default=float)
# profiles for figure / lattice comparison (Skyrme units + scales)
np.savez(os.path.join(OUT, "press_profiles.npz"), r=r, f=f, fp=fp, T00=T00, p=p, p2=p2, p4=p4, s=s, ell=ell, pu=pu,
         eps0=eps0, r0=r0, D0=D0, Mt=Mt, mN=beg["mN"])

def pf(k, d):
    print(k, {x: (round(v, 6) if isinstance(v, float) else v) for x, v in d.items() if not isinstance(v, dict)})
for k in ("P1", "P2", "P3", "P5"):
    pf(k, res[k])
print("P1e", P1["e_cebulla"])
print("BEG", P2["BEG"])
