"""kk_shape_leptons.py - TEST 3: the lepton exam. Can the lowest three rest masses of one hidden shape be
m_e : m_mu : m_tau = 1 : 206.768 : 3477.2 (Koide Q = 0.6667)? Families and parameter counts are fixed in
lab/out/kk_shape_prereg.json (T3). A ratio counts as a PREDICTION only if no parameter was fixed with it.
Writes lab/out/kk_shape_leptons.json and kk_shape_leptons.npz (scan curves for the figure)."""
import json, time, sys
import numpy as np
from scipy.optimize import brentq
import kk_shape_core as S

t0 = time.time()
R2, R3, QO = S.R2_OBS, S.R3_OBS, S.Q_OBS
out = dict(targets=dict(r2=R2, r3=R3, Q=QO))
save = {}
print(f"targets: m_mu/m_e = {R2:.3f}, m_tau/m_e = {R3:.1f}, Koide Q = {QO:.6f}")


def distinct3(E, tol=1e-7):
    E = np.sort(E)
    d = [E[0]]
    for e in E[1:]:
        if e - d[-1] > tol * max(1.0, e) and len(d) < 3:
            d.append(e)
    return np.array(d[:3])


# ---------------- E0: modes of the 2-D shape itself (0 parameters per shape; ellipse aspect = 1) ----------------
sp = json.load(open("out/kk_shape_spectra.json"))
out["E0_fixed_shapes"] = {k: dict(m2_m1=v["m2_m1"], m3_m1=v["m3_m1"], Q=v["koide"]) for k, v in sp["shapes"].items()}
asp = [1.0, 1.5, 2, 3, 4, 6, 8]
r2s, r3s, Qs = [], [], []
for a_ in asp:
    mask, info = S.shape_mask("ellipse", np.pi * 20**2, a_)
    E = S.pair_up(S.rest_masses_sparse(np.where(mask, 0.0, S.WALL), k=8))[0]
    r2s.append(E[1] / E[0]); r3s.append(E[2] / E[0]); Qs.append(S.koide(E[:3]))
out["E0_ellipse_scan"] = dict(aspect=asp, r2=r2s, r3=r3s, Q=Qs, max_r2=max(r2s), max_r3=max(r3s))
save["e0_aspect"] = np.array(asp); save["e0_r2"] = np.array(r2s); save["e0_r3"] = np.array(r3s)
print(f"[E0] walk shape modes: fixed shapes m2/m1 in "
      f"[{min(v['m2_m1'] for v in out['E0_fixed_shapes'].values()):.2f}, {max(v['m2_m1'] for v in out['E0_fixed_shapes'].values()):.2f}], "
      f"Q in [{min(v['Q'] for v in out['E0_fixed_shapes'].values()):.3f}, {max(v['Q'] for v in out['E0_fixed_shapes'].values()):.3f}]; "
      f"ellipse aspect 1..8: max m2/m1 = {max(r2s):.2f}, max m3/m1 = {max(r3s):.2f}  -> FAIL (need 206.8)")
# Ashbaugh-Benguria (1991): lambda_2/lambda_1 <= (j_{1,1}/j_{0,1})^2 for ANY 2-D Dirichlet drum
from scipy.special import jn_zeros
ab = (jn_zeros(1, 1)[0] / jn_zeros(0, 1)[0])**2
out["E0_theorem"] = dict(AB_lambda2_over_lambda1_max=ab, mass_ratio_max=np.sqrt(ab))
print(f"     theorem (Laplacian, any 2-D shape): m2/m1 = sqrt(lambda2/lambda1) <= {np.sqrt(ab):.3f}")
# charge sector n != 0 only compresses: E_i = sqrt(k0^2 + e_i^2)
e = np.array(sp["shapes"]["disk"]["R40"]["x"][:3])
k0s = np.array([0, 1, 2, 5])
out["E0_charge_sector_disk"] = {int(k): (np.sqrt(k**2 + e**2) / np.sqrt(k**2 + e[0]**2)).tolist() for k in k0s}
print(f"     charge sector (k0 R = 0,1,2,5): m2/m1 = {[round(out['E0_charge_sector_disk'][int(k)][1], 3) for k in k0s]}")

# ---------------- E1: 'mass = curvature' on a hidden loop: theta(s) = mu kappa(s)/<kappa> (1 parameter) ----------------
N = 200
loops = [("ellipse", 2.0), ("trefoil", 0.2), ("stadium", 1.0), ("egg", (0.15, 0.10))]
mus = np.geomspace(1e-4, 1.0, 30)
out["E1"] = {}
for kind, par in loops:
    kap, _ = S.loop_sample(kind, par, N)
    rel = kap / kap.mean()
    rows = []
    for mu in mus:
        m = distinct3(S.loop_masses(mu * rel, Emax=1.2))
        rows.append((m[1] / m[0], m[2] / m[0], S.koide(m)))
    rows = np.array(rows)
    res = dict(mu=mus.tolist(), r2=rows[:, 0].tolist(), r3=rows[:, 1].tolist(), Q=rows[:, 2].tolist())
    # fit mu to r2 (if reachable) and predict r3, Q
    f = lambda lm: (lambda m: m[1] / m[0])(distinct3(S.loop_masses(np.exp(lm) * rel, Emax=1.2))) - R2
    idx = np.where(np.diff(np.sign(rows[:, 0] - R2)) != 0)[0]
    res["fits"] = []
    for ii in idx:                       # EVERY crossing is reported (no choosing the nicest one)
        lm = brentq(f, np.log(mus[ii]), np.log(mus[ii + 1]), xtol=1e-12)
        m = distinct3(S.loop_masses(np.exp(lm) * rel, Emax=1.2))
        res["fits"].append(dict(mu=float(np.exp(lm)), mu_N=float(np.exp(lm) * N), pred_r3=float(m[2] / m[0]),
                                pred_Q=S.koide(m), min_kappa_rel=float(rel.min())))
        print(f"[E1 mass=curvature] {kind:8s}: fit mu = {np.exp(lm):.3e} (mu N = {np.exp(lm) * N:.3f}) to m2/m1 = 206.77 -> "
              f"PREDICTS m3/m1 = {m[2] / m[0]:.1f} (obs 3477.2), Q = {S.koide(m):.4f} (obs 0.6667)"
              f"{'  [curvature < 0 somewhere: theta changes sign]' if rel.min() < 0 else ''}")
    if not len(idx):
        print(f"[E1 mass=curvature] {kind:8s}: m2/m1 never reaches 206.77 (range {rows[:, 0].min():.2f}..{rows[:, 0].max():.2f})")
    out["E1"][kind] = res
    save[f"e1_{kind}_r2"] = rows[:, 0]; save[f"e1_{kind}_r3"] = rows[:, 1]
    sys.stdout.flush()
save["e1_mu"] = mus

# ---------------- E1 +-: signed mass theta = mu (kappa/<kappa> - 1) -> domain walls ----------------
out["E1pm"] = {}
for kind, par in loops:
    kap, _ = S.loop_sample(kind, par, N)
    rel = kap / kap.mean() - 1
    row = {}
    for mu in (0.05, 0.2, 0.5):
        m = np.sort(S.loop_masses(mu * rel, Emax=1.2))[:6]
        row[mu] = m.tolist()
    out["E1pm"][kind] = row
    print(f"[E1 signed] {kind:8s}: lowest masses at mu = 0.2: {np.array2string(np.array(row[0.2]), precision=3)}  "
          f"(sum of theta = {np.sum(rel):.1e})")

# ---------------- E2: exponential hierarchy from localization (domain walls on the walk) ----------------
# loop: background theta = +M; three segments with theta = -M of lengths d, 2d, 3d (ratio fixed in advance),
# separated by gaps of 50 sites. Each segment carries one light mode m_i ~ 2M exp(-M L_i).
d, gap = 10, 50
lengths = [3 * d, 2 * d, d]                  # e, mu, tau
Nl = sum(lengths) + 3 * gap
def theta_E2(M, Ls):
    th = np.full(int(sum(Ls) + 3 * gap), M)
    pos = 0
    for L in Ls:
        pos += gap // 2
        th[pos:pos + int(L)] = -M
        pos += int(L) + gap - gap // 2
    return th
def masses_E2(M, Ls=lengths):
    return distinct3(S.loop_masses(theta_E2(M, Ls), Emax=1.2), tol=1e-12)
Ms = np.linspace(0.04, 0.8, 39)
rows = []
for M in Ms:
    m = masses_E2(M)
    rows.append((m[1] / m[0], m[2] / m[0], S.koide(m), m[0]))
rows = np.array(rows)
save["e2_Md"] = Ms * d; save["e2_r2"] = rows[:, 0]; save["e2_r3"] = rows[:, 1]
slope = np.polyfit(Ms[8:30] * d, np.log(rows[8:30, 0]), 1)[0]
fM = lambda M: masses_E2(M)[1] / masses_E2(M)[0] - R2
i = np.where(np.diff(np.sign(rows[:, 0] - R2)) != 0)[0][0]
Mfit = brentq(fM, Ms[i], Ms[i + 1], xtol=1e-12)
m = masses_E2(Mfit)
out["E2"] = dict(d=d, gap=gap, lengths=lengths, slope_dlnr2_dMd=slope, fit_M=Mfit, fit_Md=Mfit * d,
                 pred_r3=float(m[2] / m[0]), pred_Q=S.koide(m), masses=m.tolist(),
                 scan=dict(Md=(Ms * d).tolist(), r2=rows[:, 0].tolist(), r3=rows[:, 1].tolist(), Q=rows[:, 2].tolist()))
print(f"[E2 localization, 1:2:3 segments] ln(m2/m1) grows linearly with M d (slope {slope:.3f}; ideal 1): EXPONENTIAL hierarchy.")
print(f"     fit M d = {Mfit * d:.3f} to m2/m1 = 206.77 -> PREDICTS m3/m1 = {m[2] / m[0]:.0f} (obs 3477; geometric (m2/m1)^2 = {R2**2:.0f}), "
      f"Q = {S.koide(m):.4f} (obs 0.6667)")
# Koide alone with the same 1-parameter family
fQ = lambda M: S.koide(masses_E2(M)) - 2 / 3
iq = np.where(np.diff(np.sign(rows[:, 2] - 2 / 3)) != 0)[0]
if len(iq):
    MQ = brentq(fQ, Ms[iq[0]], Ms[iq[0] + 1], xtol=1e-12); mq = masses_E2(MQ)
    out["E2"]["koide_fit"] = dict(Md=MQ * d, r2=float(mq[1] / mq[0]), r3=float(mq[2] / mq[0]))
    print(f"     (instead fitting Q = 2/3: M d = {MQ * d:.3f} gives m2/m1 = {mq[1] / mq[0]:.1f}, m3/m1 = {mq[2] / mq[0]:.0f}: wrong ratios)")
# 2 free parameters (segment lengths L_mu, L_tau at fixed M): ln r_i = M (L_e - L_i) + O(1) fits any two ratios
out["E2_free_lengths"] = "2 parameters for 2 ratios: exact fit by construction, no prediction"
# continuous: use the measured slope model ln r_i = M (L_e - L_i) + c_i  -> 2 lengths for 2 ratios, Q follows from the ratios
Qfit = S.koide([1, R2, R3])
out["E2_free_lengths_Q"] = Qfit
print(f"     with free segment lengths (2 parameters) both ratios fit by construction; Q is then a function of the "
      f"two fitted ratios (= {Qfit:.6f}), NOT a prediction.")

# analytic one-parameter families (for the doc)
x = np.sqrt(R2)
out["analytic"] = dict(
    geometric_r3=R2**2, geometric_Q=(1 - x + x * x) / (1 + x + x * x),
    gaussian_split_fermions_equal_spacing_r3=R2**(8 / 3),
    koide_geometric_q=((5 + np.sqrt(21)) / 2)**2,
    required_exponent_ratio=np.log(R3) / np.log(R2))
print(f"[analytic] geometric family m_i ~ q^i: r3 = r2^2 = {R2**2:.0f}, Q(q = 206.77) = {out['analytic']['geometric_Q']:.4f}; "
      f"Q = 2/3 needs q = {out['analytic']['koide_geometric_q']:.2f}; observed ln r3/ln r2 = {out['analytic']['required_exponent_ratio']:.4f}")

out["seconds"] = time.time() - t0
with open("out/kk_shape_leptons.json", "w") as fh:
    json.dump(out, fh, indent=1, default=float)
np.savez_compressed("out/kk_shape_leptons.npz", **save)
print(f"done {time.time() - t0:.0f}s")
