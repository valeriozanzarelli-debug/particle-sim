"""zz_hydrogen.py - T3: HYDROGEN in the zig-zag rule (radial form), levels, lines, orbitals, surfer.

The electron walker (mass flip theta = m c a / hbar, the INPUT m_e) in the static Coulomb phase of a point
nucleus (INPUT Z = 1, alpha), written for one Dirac partial wave kappa (lab/zz_core.py RadialWalk): the
outward walker chi+ and the inward walker chi- move one cell per tick at c, then rotate by
exp(-i dt (m sigma_z + kappa/r sigma_x - Z alpha / r)).
Why the radial form: the 3-D cubic walk must resolve hbar/mc (theta << 1) AND hold a few Bohr radii
(a0 = 137 hbar/mc): ~ (5 a0 / 0.1 lambda_C)^3 = 3e11 cells. The radial form is the same rule for a state of
given angular momentum (exact in the continuum); the 3-D walk itself is tested on the Coulomb problem at a
larger Z alpha in zz_atom3d.py, and the "listening" dynamics of the balls in zz_listen.py.

How the levels are found: a stationary state is a frequency at which the walk rings. Listening to the
balls for a time T with a lock-in at frequency w selects the states within ~1/T of w; listening for an
infinite time is the resolvent (1 - e^{i w dt} W)^{-1} = sum_n (e^{i w dt} W)^n, which is what the
shift-invert eigen-solver applies (zz_listen.py checks that finite listening converges to it).
Imaginary time / DMC is NOT used: the Dirac spectrum has no bottom (the negative-energy sea), so an
imaginary-time walk would sink into the lattice cutoff at -pi hbar/dt, not into 1s.

Lattice: theta = a / lambda_C in {0.4, 0.2, 0.1, 0.05}; energies extrapolated to a -> 0 with a fit in
theta^2 (the rule is second order). Reduced mass: the walk has an infinitely heavy nucleus; the proton's
mass enters only through mu = m_e m_p/(m_e + m_p) (Dirac-Coulomb energies scale exactly with the mass, so
this equals running the walk with theta_mu). Higher recoil and QED corrections are NOT in the rule.

Tests (tolerances in TOL, fixed before the first run):
  H1 levels (1s,2s,3s,2p1/2,3p1/2,2p3/2,3p3/2,3d3/2,3d5/2): extrapolated walk binding vs Sommerfeld-Dirac
     formula (check only, never input): relative error < 1e-5.
  H2 fine structure 2p3/2 - 2p1/2: walk vs Dirac < 1e-3 relative; vs measurement 10.96905 GHz < 0.5 %
     (the 25 MHz of QED is expected to be missing).
  H3 Lyman-alpha (2p -> 1s, centroid) vs 121.567 nm and Balmer-alpha (n=3 -> 2, centroid) vs 656.46 nm
     (vacuum): relative error < 1e-4.
  H4 orbitals: 1 - |<walk|exact>|^2 < 1e-4 for 1s, 2s, 2p (large components vs hydrogen R_nl).
  H5 1s surfer: |v| = alpha c sin(theta) at every r: density-weighted |v|/(alpha c) at the equator
     within 1e-3 of 1 (theta = 0.05), and the extrapolation within 1e-4.
  H6 Lyman-alpha Einstein coefficient from the walk's radial functions vs 6.2649e8 s^-1: < 1 %.
Run: timeout 590 python lab/zz_hydrogen.py
"""
import json, os, time
import numpy as np
import scipy.sparse.linalg as spl
from zz_core import RadialWalk, dirac_formula, ALPHA, MC2_EV, LAMBDA_C, T_UNIT, H_EVS, HC_EVNM, MP_ME, C

TOL = dict(level=1e-5, fs_dirac=1e-3, fs_meas=5e-3, line=1e-4, overlap=1e-4, surf=1e-3, surf_ex=1e-4, A=0.01)
OUT = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUT, exist_ok=True)
res, verdict = {}, []
t_start = time.time()
MU = 1 / (1 + 1 / MP_ME)
A0 = 1 / ALPHA                    # Bohr radius in lambda_C
THETAS = [0.4, 0.2, 0.1, 0.05]
RMAX = 90 * A0
KAPPAS = {-1: ["1s1/2", "2s1/2", "3s1/2"], 1: ["2p1/2", "3p1/2"], -2: ["2p3/2", "3p3/2"], 2: ["3d3/2"], -3: ["3d5/2"]}
NJ = {"1s1/2": (1, -1), "2s1/2": (2, -1), "3s1/2": (3, -1), "2p1/2": (2, 1), "3p1/2": (3, 1), "2p3/2": (2, -2),
      "3p3/2": (3, -2), "3d3/2": (3, 2), "3d5/2": (3, -3)}
MEAS = dict(fs=10.969045e9, lya=121.567, hal=656.46, A_lya=6.2649e8, ion=13.598434599702, s1s2=2466061413187.035e3)


def check(name, ok, msg):
    verdict.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def nodes(g):
    s = np.sign(g[np.abs(g) > 1e-6 * np.abs(g).max()])
    return int(np.sum(s[1:] != s[:-1]))


# ---------------- levels for every theta ----------------
E = {nm: {} for nm in NJ}
states = {}
guess = 1 - 0.5 * 1.1 * ALPHA**2          # a frequency just below the expected ground state (start of the search)
for th in THETAS:
    for kappa, names in KAPPAS.items():
        t0 = time.time()
        rw = RadialWalk(kappa, th, RMAX)
        Wm = rw.matrix()
        lam0 = np.exp(-1j * guess * th)
        w, v = spl.eigs(Wm, k=len(names) + 1, sigma=lam0, which="LM", tol=0)
        om = -np.angle(w) / th
        o = np.argsort(om)
        found = []
        for q in o:
            chi = rw.palindromic(v[:, q].reshape(rw.N, 2).T)
            g, f = rw.gf(chi)
            if om[q] >= 1:
                continue
            l = kappa if kappa > 0 else -kappa - 1
            n = nodes(g.real) + l + 1
            found.append((n, om[q], g, f))
        for nm in names:
            n_want = NJ[nm][0]
            m = [x for x in found if x[0] == n_want]
            if not m:
                print(f"  theta={th} kappa={kappa}: state {nm} not found ({[x[0] for x in found]})", flush=True)
                continue
            n, omq, g, f = m[0]
            E[nm][th] = 1 - omq                    # binding in m c^2
            if th == THETAS[-1]:
                states[nm] = (rw.r.copy(), g, f)
        print(f"  theta={th:4.2f} kappa={kappa:+d}: N={rw.N} "
              + ", ".join(f"{nm} {E[nm].get(th, np.nan) * MC2_EV:.7f} eV" for nm in names) + f"  ({time.time() - t0:.1f} s)", flush=True)
    guess = 1 - 1.02 * E["1s1/2"][th]


def extrapolate(d):
    th = np.array(sorted(d)); y = np.array([d[t] for t in th])
    # fit y = y0 + c1 th^2 + c2 th^4 on the three smallest, and on the three largest (error estimate)
    f3 = np.polyfit(th[:3]**2, y[:3], 2)
    f3b = np.polyfit(th[1:]**2, y[1:], 2)
    return f3[-1], abs(f3[-1] - f3b[-1])


res["levels"] = {}
print("=== H1 levels (binding energy, infinite nuclear mass) ===", flush=True)
Ex = {}
worst = 0.0
for nm, (n, kappa) in NJ.items():
    y0, err = extrapolate(E[nm])
    Ex[nm] = y0
    ref = 1 - dirac_formula(n, kappa)
    rel = abs(y0 / ref - 1)
    worst = max(worst, rel)
    res["levels"][nm] = dict(theta={str(t): E[nm][t] * MC2_EV for t in E[nm]}, extrapolated_eV=y0 * MC2_EV,
                             extrap_spread_eV=err * MC2_EV, dirac_eV=ref * MC2_EV, rel_err=rel)
    print(f"  {nm}: " + "  ".join(f"th={t}: {E[nm][t] * MC2_EV:.6f}" for t in sorted(E[nm], reverse=True))
          + f"  -> a->0: {y0 * MC2_EV:.7f} eV (+-{err * MC2_EV:.1e})  Dirac {ref * MC2_EV:.7f}  rel {rel:.1e}", flush=True)
check("H1 all 9 levels vs Dirac formula (a -> 0)", worst < TOL["level"], f"worst relative error {worst:.1e}")
print(f"  1s binding: walk {Ex['1s1/2'] * MC2_EV:.6f} eV; Bohr (non-relativistic) 13.605693 eV; Dirac 13.605874 eV;"
      f" with reduced mass {Ex['1s1/2'] * MU * MC2_EV:.6f} eV (measured ionization energy {MEAS['ion']} eV)", flush=True)

# ---------------- H2 fine structure ----------------
fs_inf = (Ex["2p1/2"] - Ex["2p3/2"]) * MC2_EV / H_EVS
fs_mu = fs_inf * MU
fs_dir = ((1 - dirac_formula(2, 1)) - (1 - dirac_formula(2, -2))) * MC2_EV / H_EVS
fs_th = {t: (E["2p1/2"][t] - E["2p3/2"][t]) * MC2_EV / H_EVS for t in THETAS}
res["fine_structure"] = dict(walk_inf_Hz=fs_inf, walk_mu_Hz=fs_mu, dirac_inf_Hz=fs_dir, measured_Hz=MEAS["fs"],
                             by_theta={str(t): v for t, v in fs_th.items()})
print("  fine structure by theta: " + ", ".join(f"{t}: {v / 1e9:.5f} GHz" for t, v in fs_th.items()), flush=True)
check("H2a fine structure 2p3/2-2p1/2 vs Dirac", abs(fs_inf / fs_dir - 1) < TOL["fs_dirac"],
      f"walk {fs_inf / 1e9:.5f} GHz (infinite mass) vs Dirac {fs_dir / 1e9:.5f} GHz")
check("H2b fine structure vs measurement", abs(fs_mu / MEAS["fs"] - 1) < TOL["fs_meas"],
      f"walk {fs_mu / 1e9:.5f} GHz (reduced mass) vs measured {MEAS['fs'] / 1e9:.5f} GHz: {100 * (fs_mu / MEAS['fs'] - 1):+.3f} % "
      f"({(fs_mu - MEAS['fs']) / 1e6:+.1f} MHz = QED, not in the rule)")

# ---------------- radial functions and dipoles ----------------
def overlap_dip(a, b):
    r, ga, fa = states[a]
    _, gb, fb = states[b]
    dr = r[1] - r[0]
    return np.sum((ga.real * gb.real + fa.real * fb.real) * r) * dr


def A_rate(up, lo, l_up, l_lo):
    """E1 Einstein coefficient (s^-1) from the walk's radial functions (non-relativistic angular factor
    max(l,l')/(2 l_up + 1)); omega from the extrapolated levels, reduced mass scaling mu^1."""
    w = (Ex[lo] - Ex[up])                                  # m c^2 / hbar units
    R = overlap_dip(up, lo)                                # lambda_C
    ang = max(l_up, l_lo) / (2 * l_up + 1)
    return 4.0 / 3.0 * ALPHA * w**3 * R**2 * ang / T_UNIT * MU


lines = []
for up, lo, lu, ll, wgt in (("2p3/2", "1s1/2", 1, 0, 2), ("2p1/2", "1s1/2", 1, 0, 1)):
    lines.append((up, lo, wgt))
lya_E = sum(w * (Ex[lo] - Ex[up]) for up, lo, w in lines) / 3 * MU * MC2_EV
lya = HC_EVNM / lya_E
lya_norec = HC_EVNM / (lya_E / MU)
# Balmer alpha: E1 components 3p->2s, 3s->2p, 3d->2p; intensity = (statistical population g_u of the nl level)
# x A(nl -> n'l') x relative line strength inside the multiplet (angular-momentum algebra, not fitted)
A_3p2s = A_rate("3p3/2", "2s1/2", 1, 0)
A_3s2p = A_rate("3s1/2", "2p1/2", 0, 1)
A_3d2p = A_rate("3d5/2", "2p3/2", 2, 1)
comps = [("3p3/2", "2s1/2", 6 * A_3p2s * 2 / 3), ("3p1/2", "2s1/2", 6 * A_3p2s * 1 / 3),
         ("3s1/2", "2p3/2", 2 * A_3s2p * 2 / 3), ("3s1/2", "2p1/2", 2 * A_3s2p * 1 / 3),
         ("3d5/2", "2p3/2", 10 * A_3d2p * 9 / 15), ("3d3/2", "2p3/2", 10 * A_3d2p * 1 / 15), ("3d3/2", "2p1/2", 10 * A_3d2p * 5 / 15)]
wsum = sum(c[2] for c in comps)
hal_E = sum(c[2] * (Ex[c[1]] - Ex[c[0]]) for c in comps) / wsum * MU * MC2_EV
hal = HC_EVNM / hal_E
res["lines"] = dict(lyman_alpha_nm=lya, lyman_alpha_no_recoil_nm=lya_norec, balmer_alpha_nm=hal,
                    balmer_components_nm={f"{u}-{l}": HC_EVNM / ((Ex[l] - Ex[u]) * MU * MC2_EV) for u, l, _ in comps},
                    A_3p2s=A_3p2s, A_3s2p=A_3s2p, A_3d2p=A_3d2p)
print(f"  Lyman-alpha: {lya:.4f} nm (without the proton recoil {lya_norec:.4f} nm); Balmer-alpha centroid {hal:.4f} nm; components "
      + ", ".join(f"{k} {v:.4f}" for k, v in res["lines"]["balmer_components_nm"].items()), flush=True)
check("H3a Lyman-alpha", abs(lya / MEAS["lya"] - 1) < TOL["line"], f"{lya:.4f} nm vs {MEAS['lya']} nm ({(lya / MEAS['lya'] - 1) * 1e6:+.1f} ppm)")
check("H3b Balmer-alpha", abs(hal / MEAS["hal"] - 1) < TOL["line"], f"{hal:.4f} nm vs {MEAS['hal']} nm ({(hal / MEAS['hal'] - 1) * 1e6:+.1f} ppm)")
s1s2 = (Ex["1s1/2"] - Ex["2s1/2"]) * MU * MC2_EV / H_EVS
res["lines"]["1S_2S_Hz"] = s1s2
print(f"  1S-2S: walk {s1s2:.6e} Hz vs measured {MEAS['s1s2']:.6e} Hz ({(s1s2 - MEAS['s1s2']) / 1e9:+.2f} GHz: Lamb shifts, QED)", flush=True)

A_lya = A_rate("2p3/2", "1s1/2", 1, 0)
res["A_lyman_alpha"] = A_lya
check("H6 Lyman-alpha Einstein A from the walk", abs(A_lya / MEAS["A_lya"] - 1) < TOL["A"], f"{A_lya:.4e} s^-1 vs {MEAS['A_lya']:.4e} s^-1 (lifetime {1e9 / A_lya:.3f} ns)")

# ---------------- H4 orbitals ----------------
def R_nl(n, l, r):
    from scipy.special import genlaguerre, factorial
    rho = 2 * r / (n * A0)
    norm = np.sqrt((2 / (n * A0))**3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * np.exp(-rho / 2) * rho**l * genlaguerre(n - l - 1, 2 * l + 1)(rho)


res["orbitals"] = {}
worst_ov = 0
for nm, (n, l) in (("1s1/2", (1, 0)), ("2s1/2", (2, 0)), ("2p1/2", (2, 1)), ("2p3/2", (2, 1))):
    r, g, f = states[nm]
    dr = r[1] - r[0]
    ex = r * R_nl(n, l, r)
    gg = g.real
    ov = (np.sum(gg * ex) * dr)**2 / (np.sum(gg**2) * dr * np.sum(ex**2) * dr)
    small = np.sum(np.abs(f)**2) * dr / (np.sum(np.abs(g)**2 + np.abs(f)**2) * dr)
    res["orbitals"][nm] = dict(one_minus_overlap=1 - ov, small_component_weight=small, r_peak_a0=r[np.argmax(gg**2)] / A0)
    worst_ov = max(worst_ov, 1 - ov)
    print(f"  {nm}: 1 - overlap with hydrogen R_nl = {1 - ov:.1e}; small-component weight {small:.2e}; peak of r^2|psi|^2 at {r[np.argmax(gg**2)] / A0:.3f} a0", flush=True)
check("H4 orbitals 1s, 2s, 2p", worst_ov < TOL["overlap"], f"worst 1 - |<walk|exact>|^2 = {worst_ov:.1e}")

# ---------------- H5 1s surfer ----------------
r, g, f = states["1s1/2"]
rho = np.abs(g)**2 + np.abs(f)**2
v = 2 * np.abs(g.real * f.real) / rho                   # |v| / c at the equator (sin theta = 1)
sel = (r > 0.01 * A0) & (r < 8 * A0)
vw = np.sum(v[sel] * rho[sel]) / np.sum(rho[sel]) / ALPHA
vmaxdev = np.max(np.abs(v[sel] / ALPHA - 1))
vth = {}
for th in THETAS:
    rw = RadialWalk(-1, th, 30 * A0)
    w, vv = spl.eigs(rw.matrix(), k=1, sigma=np.exp(-1j * (1 - E["1s1/2"][th]) * th), which="LM", tol=0)
    gg, ff = rw.gf(rw.palindromic(vv[:, 0].reshape(rw.N, 2).T))
    rr = np.abs(gg)**2 + np.abs(ff)**2
    s2 = (rw.r > 0.01 * A0) & (rw.r < 8 * A0)
    vth[th] = float(np.sum((2 * np.abs(gg.real * ff.real) / rr)[s2] * rr[s2]) / np.sum(rr[s2]) / ALPHA)
v_ex, v_err = extrapolate(vth)
res["surfer_1s"] = dict(v_over_alpha_c_weighted=vw, max_dev_over_r=vmaxdev, by_theta={str(k): x for k, x in vth.items()},
                        extrapolated=v_ex, f_over_g_median=float(np.median((f.real / g.real)[sel])),
                        f_over_g_dirac=-np.sqrt((1 - np.sqrt(1 - ALPHA**2)) / (1 + np.sqrt(1 - ALPHA**2))))
print(f"  1s surfer: |v|/(alpha c) at the equator, density-weighted = {vw:.6f} (theta = 0.05), max deviation over r in "
      f"[0.01, 8] a0 = {vmaxdev:.1e}; by theta {', '.join(f'{k}: {x:.6f}' for k, x in vth.items())} -> a->0 {v_ex:.7f}; "
      f"f/g median {res['surfer_1s']['f_over_g_median']:.8f} (Dirac {res['surfer_1s']['f_over_g_dirac']:.8f})", flush=True)
check("H5a 1s surfer speed alpha c sin(theta) (theta = 0.05)", abs(vw - 1) < TOL["surf"], f"{vw:.6f} alpha c = {vw * ALPHA * C:.5e} m/s (alpha c = {ALPHA * C:.5e})")
check("H5b 1s surfer speed, a -> 0", abs(v_ex - 1) < TOL["surf_ex"], f"{v_ex:.7f} alpha c")

# save radial functions for the viewer / figures (theta = 0.05, every 20th site up to 30 a0)
sel = states["1s1/2"][0] < 30 * A0
np.savez_compressed(os.path.join(OUT, "zz_radial.npz"), **{nm.replace("/", "_"): np.array([s[0][sel][::20], s[1].real[sel][::20], s[2].real[sel][::20]])
                                                           for nm, s in states.items()})
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for nm in ("1s1/2", "2s1/2", "2p3/2"):
        th = np.array(sorted(E[nm]))
        y = np.array([E[nm][t] for t in th])
        ax[0].plot(th**2, (y / Ex[nm] - 1), "o-", label=nm)
    ax[0].set_xlabel("theta^2 = (a / lambda_C)^2"); ax[0].set_ylabel("E_bind(walk)/E_bind(a->0) - 1"); ax[0].legend(fontsize=7)
    ax[0].set_title("hydrogen levels: lattice convergence")
    for nm, (n, l) in (("1s1/2", (1, 0)), ("2s1/2", (2, 0)), ("2p3/2", (2, 1))):
        r, g, f = states[nm]
        ax[1].plot(r / A0, (g.real**2 + f.real**2) * A0, label=f"walk {nm}")
        ax[1].plot(r / A0, (r * R_nl(n, l, r))**2 * A0, "k:", lw=0.8)
    ax[1].set_xlim(0, 15); ax[1].set_xlabel("r / a0"); ax[1].set_ylabel("r^2 |psi|^2"); ax[1].legend(fontsize=7)
    ax[1].set_title("orbitals (dotted: Schroedinger R_nl)")
    r, g, f = states["1s1/2"]
    s = r < 8 * A0
    ax[2].plot(r[s] / A0, (2 * np.abs(g.real * f.real) / (g.real**2 + f.real**2))[s] / ALPHA)
    ax[2].set_ylim(0.99, 1.01); ax[2].set_xlabel("r / a0"); ax[2].set_ylabel("|v| / (alpha c) at the equator")
    ax[2].set_title("1s surfer speed")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "zz_hydrogen.png"), dpi=90)
except Exception as ex:
    print("figure skipped:", ex)
res["verdict"] = verdict
res["seconds"] = time.time() - t_start
json.dump(res, open(os.path.join(OUT, "zz_hydrogen.json"), "w"), indent=1, default=float)
npass = sum(v for _, v in verdict)
print(f"=== {npass} PASS / {len(verdict) - npass} FAIL  ({time.time() - t_start:.0f} s)", flush=True)
