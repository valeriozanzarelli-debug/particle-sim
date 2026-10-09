"""dirac_hydrogen.py - hydrogen (Z = 1) with the Dirac rule.

1) Radial Dirac-Coulomb equations solved numerically (shooting, outward + inward, match of the
   Wronskian G_out F_in - G_in F_out; DOP853 in s = ln r, rtol 1e-13). Units hbar = c = m_e = 1.
       G' = -kappa/r G + (E + 1 - V) F ,  F' = kappa/r F - (E - 1 - V) G ,  V = -Z alpha / r
   The energy formula (Sommerfeld-Dirac) is used ONLY as a check, never as input.
2) Surfer velocity in the eigenstates: v = psi^dag alpha psi / psi^dag psi. For kappa = -1, m_j = +1/2:
       v = -(2 G F / (G^2 + F^2)) sin(theta) phi_hat   (purely azimuthal: circles about the spin axis)
3) A small 3-D surfer ensemble moved with the same RK4 rule (checks |v| <= c, r and z conserved).
Experimental values (cited, not fitted):
   H ionization energy 13.598 434 599 702 eV (NIST ASD, 1s centroid incl. hyperfine average);
   2S1/2-2P1/2 Lamb shift 1057.845(9) MHz (Lundeen & Pipkin 1981; 1057.829(3) Bezginov 2019);
   2S1/2-2P3/2 9911.200(12) MHz (Hagley & Pipkin 1994) -> 2P3/2 - 2P1/2 = 10969.045 MHz;
   1S-2S 2 466 061 413 187.035(10) kHz (Parthey et al. 2011).
Run: timeout 590 python lab/dirac_hydrogen.py
"""
import json, os, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from dirac_core import ALPHA, MC2_EV, LAMBDA_C, T_UNIT, C, SEED

t_start = time.time()
H_EVS = 4.135667696e-15        # Planck constant, eV s
MP_ME = 1836.15267343          # proton / electron mass
MU = 1 / (1 + 1 / MP_ME)       # reduced mass factor
Z = 1
ZA = Z * ALPHA
OUT = os.path.join(os.path.dirname(__file__), "out")
TOL = dict(formula_rel=1e-9, ratio_rel=1e-6, speed_rel=1e-6, vmax=1 + 1e-12)
verdict, res = [], {}


def check(name, ok, msg):
    verdict.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def rhs(s, y, kappa, eps):
    r = np.exp(s)
    G, F = y
    V = -ZA / r
    return [r * (-kappa / r * G + (2 + eps - V) * F), r * (kappa / r * F - (eps - V) * G)]


def shoot(eps, kappa, n_guess):
    gam = np.sqrt(kappa**2 - ZA**2)
    a = 1 / ALPHA                     # Bohr radius in lambda_C
    r0 = 1e-6 * a
    rmax = (40 + 10 * n_guess**2) * a
    rm = 1.5 * n_guess**2 * a
    y0 = [r0**gam, r0**gam * (gam + kappa) / ZA]
    so = solve_ivp(rhs, (np.log(r0), np.log(rm)), y0, args=(kappa, eps), method="DOP853", rtol=1e-13, atol=1e-300)
    E = 1 + eps
    lam = np.sqrt(1 - E**2)
    yi = [1e-200, -1e-200 * np.sqrt((1 - E) / (1 + E))]
    si = solve_ivp(rhs, (np.log(rmax), np.log(rm)), yi, args=(kappa, eps), method="DOP853", rtol=1e-13, atol=1e-300)
    Go, Fo = so.y[:, -1]
    Gi, Fi = si.y[:, -1]
    return (Go * Fi - Gi * Fo) / (np.hypot(Go, Fo) * np.hypot(Gi, Fi))


def eigen(kappa, n_target, nr):
    """find the (nr+1)-th eigenvalue for this kappa by scanning eps (no formula used)"""
    grid = -0.5 * ALPHA**2 * np.geomspace(1.2, 1 / (n_target + 0.6)**2, 120)   # from below 1s to above level n
    vals = [shoot(e, kappa, n_target) for e in grid]
    roots = []
    for i in range(len(grid) - 1):
        if np.sign(vals[i]) != np.sign(vals[i + 1]):
            roots.append(brentq(shoot, grid[i], grid[i + 1], args=(kappa, n_target), xtol=1e-22, rtol=1e-15, maxiter=200))
    # spurious sign flips can only come from the normalization; the bound states are the roots with
    # |W| small on both sides: keep the genuine ones by checking |shoot| at the root
    roots = [r for r in roots if abs(shoot(r, kappa, n_target)) < 1e-6]
    return sorted(roots)[nr]


def formula(n, j):
    k = j + 0.5
    d = k - np.sqrt(k**2 - ZA**2)
    return 1 / np.sqrt(1 + (ZA / (n - d))**2) - 1


levels = {"1s1/2": (1, -1, 0, 0.5), "2s1/2": (2, -1, 1, 0.5), "2p1/2": (2, 1, 0, 0.5), "2p3/2": (2, -2, 0, 1.5),
          "3d5/2": (3, -3, 0, 2.5)}
eps = {}
for nm, (n, kappa, nr, j) in levels.items():
    t0 = time.time()
    e = eigen(kappa, n, nr)
    ef = formula(n, j)
    eps[nm] = e
    res[nm] = dict(eps_mc2=e, eps_eV=e * MC2_EV, formula_eV=ef * MC2_EV, rel_err=abs(e / ef - 1))
    check(f"H1 {nm} numerical vs Dirac formula", abs(e / ef - 1) < TOL["formula_rel"],
          f"E - mc^2 = {e * MC2_EV:.9f} eV (formula {ef * MC2_EV:.9f}), rel diff {abs(e / ef - 1):.1e}  [{time.time() - t0:.1f} s]")

# ---- comparisons with experiment ----
nr1s = -0.5 * ALPHA**2 * MC2_EV
bind_inf = -eps["1s1/2"] * MC2_EV
bind_mu = bind_inf * MU
exp_ion = 13.598434599702
res["1s"] = dict(bind_dirac_inf=bind_inf, bohr=-nr1s, rel_corr=-(bind_inf + nr1s), bind_dirac_mu=bind_mu,
                 exp=exp_ion, diff_eV=bind_mu - exp_ion, diff_GHz=(bind_mu - exp_ion) / H_EVS / 1e9)
print(f"1s binding: Bohr (inf. mass) {-nr1s:.6f} eV; Dirac {bind_inf:.6f} eV (relativistic correction "
      f"{-(bind_inf + nr1s):.3e} eV); Dirac x reduced mass {bind_mu:.6f} eV; exp {exp_ion:.6f} eV; "
      f"Dirac - exp = {bind_mu - exp_ion:.3e} eV = {(bind_mu - exp_ion) / H_EVS / 1e9:.2f} GHz (1s Lamb shift ~ 8.17 GHz)", flush=True)
fs = (eps["2p3/2"] - eps["2p1/2"]) * MC2_EV / H_EVS / 1e9
res["fine_structure_GHz"] = dict(dirac_inf=fs, dirac_mu=fs * MU, exp=10.969045, diff_MHz=(fs * MU - 10.969045) * 1e3)
print(f"2p3/2 - 2p1/2: Dirac {fs:.5f} GHz (inf. mass), {fs * MU:.5f} GHz (reduced mass); exp 10.96905 GHz; "
      f"diff {(fs * MU - 10.969045) * 1e3:.1f} MHz", flush=True)
check("H2 fine structure within 0.5% of experiment", abs(fs * MU / 10.969045 - 1) < 5e-3,
      f"{fs * MU:.5f} vs 10.96905 GHz ({(fs * MU / 10.969045 - 1) * 100:.3f} %)")
lamb = (eps["2s1/2"] - eps["2p1/2"]) * MC2_EV / H_EVS / 1e6
res["lamb_MHz"] = dict(dirac=lamb, exp=1057.845)
print(f"2s1/2 - 2p1/2 (Lamb shift): Dirac {lamb:.3e} MHz (exactly degenerate), exp 1057.845 MHz -> NOT in Dirac", flush=True)
f1s2s = (eps["2s1/2"] - eps["1s1/2"]) * MU * MC2_EV / H_EVS
res["1s2s_Hz"] = dict(dirac_mu=f1s2s, exp=2466061413187035.0, diff_GHz=(f1s2s - 2466061413187035.0) / 1e9)
print(f"1S-2S: Dirac x mu {f1s2s:.6e} Hz, exp 2.466061413187e15 Hz, diff {(f1s2s - 2466061413187035.0) / 1e9:.2f} GHz "
      f"(QED: 1s Lamb 8.17 GHz - 2s Lamb 1.05 GHz ~ 7.1 GHz)", flush=True)


# ---- wavefunctions and surfer velocity ----
def wavefunction(kappa, e, n_guess, r):
    gam = np.sqrt(kappa**2 - ZA**2)
    a = 1 / ALPHA
    r0 = 1e-6 * a
    y0 = [r0**gam, r0**gam * (gam + kappa) / ZA]
    sol = solve_ivp(rhs, (np.log(r0), np.log(r[-1])), y0, args=(kappa, e), method="DOP853", rtol=1e-12,
                    atol=1e-300, t_eval=np.log(r))
    G, F = sol.y
    return G, F


a0 = 1 / ALPHA
r = np.geomspace(1e-4, 8, 2000) * a0
prof = {}
for nm, kappa, ng, rcut in (("1s1/2", -1, 1, 8), ("2s1/2", -1, 2, 8)):
    G, F = wavefunction(kappa, eps[nm], ng, r)
    m = r < rcut * a0 * (0.6 if nm == "1s1/2" else 1.0)   # outward integration of a bound state diverges at large r
    v_eq = 2 * G * F / (G**2 + F**2)          # at theta = 90 deg, |v| = |v_eq| sin(theta)
    prof[nm] = (r[m], G[m], F[m], v_eq[m])
gam = np.sqrt(1 - ZA**2)
ratio_th = -np.sqrt((1 - gam) / (1 + gam))
r1, G1, F1, v1 = prof["1s1/2"]
rat = F1 / G1
check("H3 1s: F/G constant = -sqrt((1-gamma)/(1+gamma))", np.max(np.abs(rat / ratio_th - 1)) < TOL["ratio_rel"],
      f"F/G in [{rat.min():.10f}, {rat.max():.10f}], theory {ratio_th:.10f}")
check("H4 1s surfer speed = Z alpha c sin(theta), independent of r", np.max(np.abs(np.abs(v1) / ZA - 1)) < TOL["speed_rel"],
      f"|v|/(alpha c) at theta=90 deg in [{np.abs(v1).min() / ZA:.9f}, {np.abs(v1).max() / ZA:.9f}] -> {ZA * C:.4e} m/s")
r2, G2, F2, v2 = prof["2s1/2"]
res["surfer_1s"] = dict(speed_over_alpha_min=float(np.abs(v1).min() / ZA), speed_over_alpha_max=float(np.abs(v1).max() / ZA),
                        speed_SI_equator=ZA * C, mean_speed_over_c=ZA * np.pi / 4,
                        omega_at_a0_rad_s=ZA / a0 / T_UNIT, period_at_a0_s=2 * np.pi * a0 / ZA * T_UNIT,
                        direction="sign(G F) < 0 -> v along +phi_hat (counter-clockwise about +z, spin up)")
res["surfer_2s"] = dict(max_speed_over_c=float(np.abs(v2).max()), r_of_max_over_a0=float(r2[np.argmax(np.abs(v2))] / a0))
print(f"1s surfer: circles about the spin axis, speed alpha c sin(theta) = {ZA * C:.4e} m/s * sin(theta); "
      f"angular velocity alpha c / r: at r = a0 {ZA / a0 / T_UNIT:.3e} rad/s (period {2 * np.pi * a0 / ZA * T_UNIT:.3e} s)", flush=True)
print(f"2s surfer: speed at the equator peaks at {np.abs(v2).max():.4f} c at r = {r2[np.argmax(np.abs(v2))] / a0:.4f} a0 "
      f"(near the node of the large component, where |G| = |F|)", flush=True)

# Pauli (non-relativistic + spin term) velocity for 1s: v = (1/2m rho) curl(rho z_hat) -> |dln rho/dr|/2 sin(theta)
v_pauli = 1.0 / a0      # rho ~ exp(-2 r / a0): |d ln rho / dr| / 2 = 1/a0 = alpha
res["pauli_1s_speed_over_alpha"] = v_pauli / ZA
print(f"Pauli limit (Schroedinger 1s + spin current): |v| = {v_pauli / ZA:.6f} alpha c sin(theta); "
      f"Schroedinger without spin term: v = 0 (real wave function)", flush=True)

# ---- 3-D surfer ensemble with the RK4 rule in the 1s state ----
rng = np.random.default_rng(SEED)
N = 2000
# sample r from r^2 (G^2+F^2) on the profile grid, isotropic angles (1s density is spherical)
w = (G1**2 + F1**2) * np.gradient(r1)
w /= w.sum()
rs = rng.choice(r1, size=N, p=w)
ct = rng.uniform(-1, 1, N); ph = rng.uniform(0, 2 * np.pi, N)
P = np.array([rs * np.sqrt(1 - ct**2) * np.cos(ph), rs * np.sqrt(1 - ct**2) * np.sin(ph), rs * ct])
lnr = np.log(r1)


def vel(P):
    rr = np.linalg.norm(P, axis=0)
    vq = np.interp(np.log(rr), lnr, v1)           # 2GF/(G^2+F^2), constant for 1s, interpolated anyway
    # v = -vq * sin(theta) phi_hat ; sin(theta) phi_hat = (-y, x, 0)/r
    return -vq * np.array([-P[1], P[0], np.zeros_like(rr)]) / rr


T = 2 * np.pi * a0 / ZA     # one period at r = a0
nst = 4000
dt = T / nst
vmax = 0.0
P0 = P.copy()
for _ in range(nst):
    k1 = vel(P); k2 = vel(P + dt / 2 * k1); k3 = vel(P + dt / 2 * k2); k4 = vel(P + dt * k3)
    vmax = max(vmax, max(np.max(np.linalg.norm(k, axis=0)) for k in (k1, k2, k3, k4)))
    P = P + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
drz = np.max(np.abs(P[2] - P0[2])) / a0
drr = np.max(np.abs(np.hypot(P[0], P[1]) - np.hypot(P0[0], P0[1]))) / a0
turns = (np.arctan2(P[1], P[0]) - np.arctan2(P0[1], P0[0]))
check("H5 1s surfer ensemble: |v| <= c, circles (z, r_perp conserved)", vmax <= TOL["vmax"] and drz < 1e-9 and drr < 1e-6,
      f"max|v| = {vmax:.6e} c, max dz = {drz:.1e} a0, max d(r_perp) = {drr:.1e} a0 over one period at a0")
res["ensemble_1s"] = dict(vmax=vmax, dz=drz, dr=drr)

# ---- figure: 1s circulation ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2), dpi=85)
xg = np.linspace(-3, 3, 301)
XX, ZZ = np.meshgrid(xg, xg)
RR = np.hypot(XX, ZZ) * a0
dens = np.interp(np.log(np.maximum(RR, r1[0])), lnr, G1**2 + F1**2) / np.maximum(RR, r1[0])**2
vphi = -np.interp(np.log(np.maximum(RR, r1[0])), lnr, v1) * XX / np.maximum(np.hypot(XX, ZZ), 1e-12)
im = ax[0].pcolormesh(XX, ZZ, vphi / ZA, cmap="RdBu_r", vmin=-1, vmax=1, shading="auto")
ax[0].contour(XX, ZZ, np.log10(dens / dens.max() + 1e-30), levels=[-3, -2, -1, -0.3], colors="k", linewidths=0.6)
ax[0].set_aspect("equal"); ax[0].set_xlabel("x [a0]"); ax[0].set_ylabel("z [a0] (asse dello spin)")
ax[0].set_title("1s Dirac: v_phi/(alpha c) nel piano x-z\n(rosso esce dal foglio; contorni: densita)")
fig.colorbar(im, ax=ax[0], fraction=0.046)
ax[1].semilogx(r1 / a0, np.abs(v1) / ZA, label="1s Dirac: |v|/(alpha c)")
ax[1].semilogx(r2 / a0, np.abs(v2) / ZA, label="2s Dirac: |v|/(alpha c)")
ax[1].axhline(1, color="k", lw=0.5, ls=":")
ax[1].set_xlabel("r [a0]"); ax[1].set_ylabel("velocita all'equatore / (alpha c)")
ax[1].set_title("velocita del surfista (theta = 90 gradi)\nSchroedinger senza spin: 0"); ax[1].legend(fontsize=8)
ax[1].set_ylim(0, 1.05 / ZA if np.abs(v2).max() > 2 * ZA else 2)
ax[1].set_yscale("log"); ax[1].set_ylim(1e-3, 300)
Om = np.abs(v1) / r1 / T_UNIT
ax[2].loglog(r1 / a0, Om, label="1s: Omega = alpha c / r")
ax[2].axhline(2 / T_UNIT, color="r", ls="--", lw=0.8, label="Zitterbewegung 2mc^2/hbar")
ax[2].axvline(0.5 * ALPHA, color="gray", ls=":", lw=0.8, label="r = hbar/2mc")
ax[2].set_xlabel("r [a0]"); ax[2].set_ylabel("velocita angolare [rad/s]")
ax[2].set_title("giri del surfista nello stato 1s"); ax[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "dirac_h1s_circulation.png")); plt.close(fig)

res["verdict"] = verdict
res["runtime_s"] = time.time() - t_start
with open(os.path.join(OUT, "dirac_hydrogen.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
print(f"=== {sum(v for _, v in verdict)}/{len(verdict)} PASS, runtime {time.time() - t_start:.1f} s ===", flush=True)
