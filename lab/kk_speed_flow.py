"""kk_speed_flow.py - TASK 3: can different limiting speeds flow to a common one ("algorithm back to
normal speed")? Chadha-Nielsen (Nucl. Phys. B217 (1983) 125), Anber-Donoghue (PRD 83 (2011) 105027).

(i) 1+1 D, two chains (speeds c1, c2), exact (Gaussian) coarse-graining: block size b keeps k <= pi/b.
    A: gradient coupling g (u'v') = forward-scattering density-density coupling of 1-D fermions after
       bosonization (Luttinger): speeds at every block size.
    B: on-site mixing mu^2 (u - v)^2 / 2: the low band at every block size.
(ii) 3+1 D, Yukawa g phi psibar psi, fermion speed c_f, scalar speed c_s, one loop.
    Scalar: the fermion loop is Lorentz invariant with speed c_f (rescale q -> q / c_f), so its log term
    is exactly (g^2 / (4 pi^2 c_f^3)) (p0^2 + c_f^2 p^2) ln:   d c_s^2 / dl = a (c_f^2 - c_s^2), a = g^2/(4 pi^2 c_f^3).
    Fermion: log coefficients of gamma0 p0 and c_f gamma.p computed NUMERICALLY on the unit 3-sphere
    (integrand homogeneous of degree -4):  d c_f / dl = c_f g^2 (I1 - I0) / (2 pi)^4.
    Check: c_f = c_s -> I0 = I1 = pi^2 (gamma_psi = g^2/(32 pi^2), the textbook value).
    Coupling: Lorentz-invariant one-loop running dg^2/dl = -5 g^4 / (8 pi^2) (APPROXIMATION, flagged).
    l = ln(Lambda / mu) grows towards the infrared. Planck -> electron mass: l = 51.5.
Writes out/kk_speed_flow.json. Run from lab/.
"""
import json, time
import numpy as np
from scipy.integrate import solve_ivp, quad

t0 = time.time()
out = {}

# ---------------- (i) 1+1 D chains ----------------
c1, c2 = 1.0, 0.5
bs = 2.0 ** np.arange(0, 13)
rowsA, rowsB = [], []
for g in (0.3, 0.8):
    for b in bs:
        k = np.pi / b
        M = 4 * np.sin(k / 2)**2 * np.array([[c1**2, g * c1 * c2], [g * c1 * c2, c2**2]])
        w = np.sqrt(np.linalg.eigvalsh(M))
        rowsA.append(dict(g=g, block=b, v_slow=float(w[0] / k), v_fast=float(w[1] / k)))
for mu in (0.05,):
    for b in bs:
        k = np.pi / b
        M = 4 * np.sin(k / 2)**2 * np.diag([c1**2, c2**2]) + mu**2 * np.array([[1, -1], [-1, 1]])
        w = np.sqrt(np.linalg.eigvalsh(M))
        rowsB.append(dict(mu=mu, block=b, v_low=float(w[0] / k), gap_band=float(w[1])))
vA_exact = {g: np.sqrt(np.linalg.eigvalsh(np.array([[c1**2, g * c1 * c2], [g * c1 * c2, c2**2]]))).tolist() for g in (0.3, 0.8)}
out["i_chains"] = dict(c1=c1, c2=c2, A_gradient=rowsA, A_exact_k0=vA_exact, B_mass_mixing=rowsB,
                       B_common_speed=float(np.sqrt((c1**2 + c2**2) / 2)))
for r in rowsA:
    if r["block"] in (1, 64, 4096):
        print(f"(i-A) g={r['g']}: block {int(r['block']):5d}: speeds {r['v_slow']:.4f}, {r['v_fast']:.4f}")
for r in rowsB:
    if r["block"] in (1, 16, 256, 4096):
        print(f"(i-B) mu={r['mu']}: block {int(r['block']):5d}: low-band speed {r['v_low']:.4f} (sqrt((c1^2+c2^2)/2) = {np.sqrt((c1**2+c2**2)/2):.4f}); other band gapped at {r['gap_band']:.3f}")


# ---------------- (ii) 3+1 D Yukawa one-loop ----------------
def I_fermion(cf, cs):
    """log coefficients (unit 3-sphere integrals) of gamma0 p0 and cf gamma.p in Sigma/(-i g^2)"""
    def F(chi, which):
        q0, qs = np.cos(chi), np.sin(chi)
        u = q0**2 + cf**2 * qs**2
        D = 1.0 / (q0**2 + cs**2 * qs**2)
        f = (1 / u - 2 * q0**2 / u**2) if which == 0 else (1 / u - 2 * cf**2 * (qs**2 / 3) / u**2)
        return 4 * np.pi * np.sin(chi)**2 * f * D
    return (quad(F, 0, np.pi, args=(0,), epsabs=1e-13, epsrel=1e-12)[0],
            quad(F, 0, np.pi, args=(1,), epsabs=1e-13, epsrel=1e-12)[0])


I0, I1 = I_fermion(1.0, 1.0)
out["ii_check_lorentz"] = dict(I0=I0, I1=I1, pi2=np.pi**2)
print(f"(ii) check c_f = c_s = 1: I0 = {I0:.12f}, I1 = {I1:.12f}, pi^2 = {np.pi**2:.12f}")
for r in (0.5, 0.9, 1.1, 2.0):
    a, b = I_fermion(1.0, r)
    print(f"(ii) c_f = 1, c_s = {r}: I1 - I0 = {b - a:+.5f}  -> fermion speed {'rises' if b > a else 'falls'} (towards c_s)")


def rhs(l, y):
    cf, cs, g2 = y
    i0, i1 = I_fermion(cf, cs)
    return [cf * g2 * (i1 - i0) / (2 * np.pi)**4,
            (g2 / (4 * np.pi**2 * cf**3)) * (cf**2 - cs**2) / (2 * cs),
            -5 * g2**2 / (8 * np.pi**2)]


lP = 51.5
flows = []
for cs0 in (0.5, 2.0, 1.1):
    for g20 in (1.0, 5.0, 20.0):
        ls = np.concatenate([[0], np.geomspace(0.1, 1e4, 120)])
        sol = solve_ivp(rhs, (0, 1e4), [1.0, cs0, g20], t_eval=ls, rtol=1e-8, atol=1e-12)
        d = sol.y[1] / sol.y[0] - 1
        dP = float(np.interp(lP, sol.t, d))
        big = sol.t > 1e3
        p = -np.polyfit(np.log(sol.t[big]), np.log(np.abs(d[big])), 1)[0]
        flows.append(dict(cs0=cs0, g2_0=g20, l=sol.t.tolist(), delta=d.tolist(), g2=sol.y[2].tolist(),
                          delta_at_planck_to_me=dP, reduction=float(d[0] / dP), exponent_large_l=float(p)))
        print(f"(ii) start c_s/c_f = {cs0}, g^2 = {g20:4.1f}: delta = c_s/c_f - 1: {d[0]:+.3f} -> {dP:+.4f} after 51.5 e-folds"
              f" (Planck -> m_e; x{d[0]/dP:.1f} smaller); large-l law delta ~ l^-{p:.3f}")
out["ii_flows"] = flows
# linearized asymptotic law: g^2 -> 8 pi^2 / (5 l), d delta/dl = -(1/(4 pi^2) + kappa/(2 pi)^4) g^2 delta
h = 1e-4
dIp = I_fermion(1.0, 1.0 + h); dIm = I_fermion(1.0, 1.0 - h)
kap = ((dIp[1] - dIp[0]) - (dIm[1] - dIm[0])) / (2 * h)          # d(I1 - I0)/d delta at delta = 0
p_lin = (8 * np.pi**2 / 5) * (1 / (4 * np.pi**2) + kap / (2 * np.pi)**4)
out["ii_linear_exponent"] = dict(kappa_dI_ddelta=float(kap), p=float(p_lin))
out["ii_efolds_needed_1e-15_from_delta1"] = float(10 ** (15 / p_lin))
print(f"(ii) linearized: d(I1-I0)/d delta = {kap:.3f}; delta ~ l^-{p_lin:.3f} (pure power of the log)")
print(f"(ii) from delta = 1 to the bound 1e-15 needs l ~ 10^{15/p_lin:.0f} e-folds, i.e. a scale ratio e^(10^{15/p_lin:.0f}); Planck -> m_e is l = 51.5")

# ---------------- QED running of a mass from Planck to m_e (fixed alpha, one loop) ----------------
alpha = 7.2973525693e-3
fac = float(np.exp(3 * alpha / (2 * np.pi) * lP))
out["qed_mass_running_factor"] = fac
print(f"QED one-loop mass running Planck -> m_e: factor {fac:.3f} (not 1e21)")
out["time_s"] = time.time() - t0
json.dump(out, open("out/kk_speed_flow.json", "w"), indent=1)
print("time %.1f s" % out["time_s"])
