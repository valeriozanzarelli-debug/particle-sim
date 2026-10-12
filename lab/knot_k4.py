"""knot_k4.py - K4: moment of inertia, collective quantization (Adkins-Nappi-Witten 1983), what is a prediction and
what is a fit. Uses the radial ODE and the lattice runs of knot_k1.py (lab/out/knot_k1_*.json).
Spin 1/2: Finkelstein-Rubinstein (pi_1 of the configuration space = Z2) is CITED, not derived here.
"""
import os, sys, json, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot_core as kc

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
P0 = kc.radial_profile()          # shooting (used for the prereg numbers 103.90)
P = kc.radial_bvp()               # collocation, exact tail (reference)
ld = lambda tag: json.load(open(os.path.join(OUT, f"knot_k1_{tag}.json")))
a, b = ld("h0.2_L48"), ld("h0.15_L64")           # same half-box 4.8
rich = lambda x2, x1: x1 + (x1 - x2) * 0.15**2 / (0.2**2 - 0.15**2)
lat = dict(E=rich(a["E"], b["E"]), Lam=rich(a["Lambda"], b["Lambda"]), rB2=rich(a["baryon"]["r2"], b["baryon"]["r2"]),
           Lam_h02=a["Lambda"], Lam_h015=b["Lambda"])

MN, MD, hbarc = 938.92, 1232.0, 197.3269804
expt = dict(F_pi_MeV=184.2, mu_p=2.7928, mu_n=-1.9130, mu_ratio=2.7928 / 1.9130, g_A=1.2754, r_I0_fm=0.77)


def anw(E, Lam, rB2, gA_D=None):
    """E, Lam, <r^2>_B in Skyrme units. Returns ANW-style numbers with M_N, M_Delta as the two inputs."""
    lam = 3 / (2 * (MD - MN))                     # 1/MeV (rotor: M_Delta - M_N = 3/(2 lambda))
    M = MN - 3 / (8 * lam)
    e = (E * Lam / (4 * M * lam)) ** 0.25        # M = (E/4) F_pi/e,  lambda = Lam/(e^3 F_pi)
    Fpi = 4 * M * e / E
    r2 = rB2 * 4 / (e * Fpi) ** 2                 # MeV^-2
    mu_p = r2 * MN * (MD - MN) / 9 + MN / (2 * (MD - MN))
    mu_n = r2 * MN * (MD - MN) / 9 - MN / (2 * (MD - MN))
    out = dict(e=e, F_pi_MeV=Fpi, M_classical=M, r_I0_fm=math.sqrt(r2) * hbarc, mu_p=mu_p, mu_n=mu_n,
               mu_ratio=abs(mu_p / mu_n), mu_NDelta=(mu_p - mu_n) / math.sqrt(2),
               F_pi_over_MN=Fpi / MN)
    if gA_D is not None:
        out["g_A"] = -math.pi / (3 * e**2) * gA_D
    return out


# g_A integral (ANW units x = 2 r)
x = 2 * P["r"][1:]; F = P["f"][1:]; Fp = np.gradient(F, x); s2 = np.sin(F) ** 2; s2F = np.sin(2 * F)
D = np.trapezoid(x**2 * Fp + x * s2F + 4 * (x * s2F * Fp**2 + 2 * s2 * Fp + s2 * s2F / x), x)
rad = anw(P["E"], P["Lam"], P["rB2"], D)
latp = anw(lat["E"], lat["Lam"], lat["rB2"], D)

# e from ONE ratio only: M_Delta/M_N fixes e (dimensionless) -> mu_p/mu_n and g_A do not need F_pi
# (check: they are identical to the 2-input numbers above, as they must be)
zero_input = dict(E_over_12pi2=P["E"] / (12 * math.pi**2), Lambda=P["Lam"], rE_over_rB=math.sqrt(P["rE2"] / P["rB2"]),
                  rotor_ratio=(35 / 4 - 15 / 4) / (15 / 4 - 3 / 4), muNDelta_over_mup_minus_mun=1 / math.sqrt(2))
rad_box = {}
for tag in ("h0.2_L32", "h0.2_L48", "h0.2_L64", "h0.15_L32", "h0.15_L48", "h0.15_L64"):
    if os.path.exists(os.path.join(OUT, f"knot_k1_{tag}.json")):
        d = ld(tag); rad_box[tag] = dict(half_box=d["half_box"], Lambda=d["Lambda"], E=d["E"])
res = dict(radial=dict(E=P["E"], E2=P["E2"], E4=P["E4"], Lam=P["Lam"], rB2=P["rB2"], rE2=P["rE2"], C_tail=P["C_tail"], D=D),
           radial_shooting_prereg=dict(E=P0["E"], Lam=P0["Lam"], E2=P0["E2"], E4=P0["E4"]), lattice_boxes=rad_box, lattice=lat, anw_radial=rad, anw_lattice=latp,
           zero_input=zero_input, experiment=expt,
           literature=dict(ANW_mu_p=1.87, ANW_mu_n=-1.31, ANW_ratio=1.43, ANW_g_A=0.61, ANW_F_pi=129, ANW_e=5.45,
                           note="ANW numbers quoted from memory/standard reviews; our ODE reproduces them (see doc)"),
           exams=dict(K4a=abs(lat["Lam"] / P["Lam"] - 1), K4a_vs_prereg_103p90=abs(lat["Lam"] / 103.90 - 1),
                      K4b_radial=abs(rad["mu_ratio"] / expt["mu_ratio"] - 1), K4b=abs(latp["mu_ratio"] / expt["mu_ratio"] - 1),
                      K4c_gA=abs(latp.get("g_A", rad["g_A"]) / expt["g_A"] - 1), K4c_Fpi=abs(latp["F_pi_MeV"] / expt["F_pi_MeV"] - 1)))
res["pass"] = {k: bool(res["exams"][k] < lim) for k, lim in (("K4a", 0.02), ("K4b", 0.03), ("K4b_radial", 0.03), ("K4c_gA", 0.30), ("K4c_Fpi", 0.30))}
# informative (post-hoc): Lambda has a 1/R box tail (r^2 sin^2 f ~ C^2/r^2); fit Lambda(R) = Lam_inf - c/R from L = 48, 64 at each h
ext = {}
for hh in ("0.2", "0.15"):
    p1, p2 = rad_box[f"h{hh}_L48"], rad_box[f"h{hh}_L64"]
    c = (p2["Lambda"] - p1["Lambda"]) / (1 / p1["half_box"] - 1 / p2["half_box"])
    ext[hh] = dict(Lam_inf=p2["Lambda"] + c / p2["half_box"], c=c, c_expected_sphere=16 * math.pi / 3 * P["C_tail"]**2)
ext["richardson"] = rich(ext["0.2"]["Lam_inf"], ext["0.15"]["Lam_inf"])
ext["dev_from_radial"] = ext["richardson"] / P["Lam"] - 1
res["lattice_Lambda_box_extrapolated_posthoc"] = ext
res["anw_lattice_boxextrap_posthoc"] = anw(lat["E"], ext["richardson"], lat["rB2"], D)
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(OUT, "knot_k4.json"), "w"), indent=1)
