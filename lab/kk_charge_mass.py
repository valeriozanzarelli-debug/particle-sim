"""kk_charge_mass.py - tests 1 (charge quantized + conserved) and 2 (rest-energy cost of n != 0, KK numbers).
Seeded; writes lab/out/kk_charge_mass.json."""
import json, time
import numpy as np
import kk_core as K
from zz_core import bloch_batch

t0 = time.time()
rng = np.random.default_rng(K.SEED)
out = {}

# ---- TEST 1a: a RANDOM state on the hidden circle -> its spectrum in time is discrete ----
# spatially uniform (nx = ny = 1), N = 32 hidden sites, theta = m a = 0.1 (lambda_C = 10 cells)
N, th, T = 32, 0.1, 4096
W = K.KKWalk(1, 1, N, th)
psi = rng.normal(size=(4, 1, 1, N)) + 1j * rng.normal(size=(4, 1, 1, N))
psi /= np.linalg.norm(psi)
sig = np.empty((T, 4, N), complex)
for t in range(T):
    sig[t] = psi[:, 0, 0, :]
    psi = W.tick(psi)
# frequencies present: FFT in time of every component, keep peaks
# the rule alternates the axis order (x,y,w then w,y,x): the true period is TWO ticks, sample every 2
sig = sig[::2]; T = sig.shape[0]
spec = np.sum(np.abs(np.fft.fft(sig * np.hanning(T)[:, None, None], axis=0))**2, axis=(1, 2))
om = -2 * np.pi * np.fft.fftfreq(T) / 2                   # angular frequency per tick
# predicted lattice levels for integer n (two-tick Bloch matrix, k_x = k_y = 0)
R = N / (2 * np.pi)
ns = np.arange(-N // 2 + 1, N // 2 + 1)
KA = np.stack([0 * ns, 0 * ns, 2 * np.pi * ns / N], 1).astype(float)
lev = np.sort(np.abs(np.angle(np.linalg.eigvals(bloch_batch(KA, th))) / 2).ravel())
pk = [i for i in range(1, T - 1) if spec[i] > spec[i - 1] and spec[i] > spec[i + 1] and spec[i] > 1e-4 * spec.max()]
pk_om = np.sort(np.abs(om[pk]))
dist = [float(np.min(np.abs(lev - w))) for w in pk_om]
# 'half-integer' candidates: where would n = 0.5, 1.5 sit? none of the peaks there
half = np.sqrt(th**2 + ((np.array([0.5, 1.5, 2.5])) / R)**2)
out["test1a"] = dict(N=N, theta=th, ticks=T, n_peaks=len(pk), freq_resolution=np.pi / T,
                     max_dist_peak_to_integer_level=max(dist),
                     min_dist_peak_to_half_integer=float(np.min(np.abs(pk_om[:, None] - half[None]))),
                     lowest_peaks=pk_om[:6].tolist(), integer_levels_n0to3=[float(np.sqrt(th**2 + (n / R)**2)) for n in range(4)])
print("TEST 1a  random hidden profile: %d spectral peaks; max |peak - integer-n level| = %.2e (resolution %.2e); "
      "closest peak to a half-integer level = %.2e" % (len(pk), max(dist), np.pi / T, out["test1a"]["min_dist_peak_to_half_integer"]))

# ---- TEST 1b: charge of a packet is integer and conserved (random potential + random twist) ----
nx = ny = 32; N = 6; th = 0.2
res = {}
for label, lumpy in (("smooth hidden circle", False), ("LUMPY hidden circle", True)):
    for n in (1, -1, 2):
        psi = K.packet(nx, ny, N, th, n, k0=(0.4, 0.2), sigma=4)
        V2 = 0.3 * rng.random((nx, ny, 1))
        V = V2 + (0.3 * rng.random((1, 1, N)) * (1 + 0 * V2) if lumpy else 0)
        D = {0: 0.5 * rng.normal(size=(nx, ny)), 1: 0.5 * rng.normal(size=(nx, ny))}
        Wk = K.KKWalk(nx, ny, N, th, V=V, D=D)
        q0 = K.charge(psi)[0]
        qs = []
        for t in range(300):
            psi = Wk.tick(psi)
            qs.append(K.charge(psi)[0])
        drift = float(np.max(np.abs(np.array(qs) - q0)))
        res[f"{label} n={n}"] = dict(q0=q0, max_drift=drift, norm_err=float(abs(np.sum(np.abs(psi)**2) - 1)))
        print(f"TEST 1b  {label:22s} n={n:+d}: Q(0) = {q0:+.12f}, max |Q(t)-Q(0)| over 300 ticks = {drift:.2e}")
out["test1b"] = res

# ---- TEST 2: rest energy of mode n vs hidden radius, measured from the walk ----
th = 0.1
rows = []
for N in (16, 32, 63, 126, 252):
    R = N / (2 * np.pi)
    for n in (0, 1, 2):
        psi = np.zeros((4, 1, 1, N), complex)
        # positive-energy rest state of mode n
        psi = K.packet(1, 1, N, th, n, sigma=0.0)
        Wk = K.KKWalk(1, 1, N, th)
        p0 = psi.copy(); Tm = 2000
        ov = np.empty(Tm, complex)
        for t in range(Tm):
            ov[t] = np.vdot(p0, psi)
            psi = Wk.tick(psi)
        # two-tick phase -> energy per tick (exact eigenstate so the phase is clean)
        E = float(-np.angle(np.mean(ov[2::2] / ov[:-2:2])) / 2)
        cont = float(np.sqrt(th**2 + (n / R)**2))
        rows.append(dict(N=N, R_over_lambdaC=R * th, n=n, E_walk=E, E_continuum=cont,
                         E_over_m=E / th, rel_err=(E - cont) / cont))
        print(f"TEST 2  N={N:3d} R={R*th:5.2f} lambda_C  n={n}: E_walk = {E:.6f}, sqrt(m^2+(n/R)^2) = {cont:.6f}"
              f"  (E/m = {E/th:6.3f}, rel.err {(E-cont)/cont:+.1e})")
out["test2_walk"] = rows

# ---- KK numbers (SI, CODATA) ----
hbar, c, G, me, eV = 1.054571817e-34, 299792458.0, 6.67430e-11, 9.1093837015e-31, 1.602176634e-19
alpha = 7.2973525693e-3
lP = np.sqrt(hbar * G / c**3)
R_min = hbar / (me * c)                                   # n=1 costs <= m_e c^2  <=> hbar/(R c) <= m_e c
R_KK = 2 * lP / np.sqrt(alpha)                             # KK: alpha = 4 G hbar / (c^3 R^2)  (e = sqrt(16 pi G) hbar/(R c), Heaviside-Lorentz)
M_KK_GeV = hbar * c / R_KK / eV / 1e9
alpha_at_Rmin = 4 * lP**2 / R_min**2
out["kk_numbers"] = dict(l_Planck_m=lP, R_min_for_cost_le_me_m=R_min, R_KK_for_alpha_m=R_KK,
                         first_mode_mass_GeV=M_KK_GeV, ratio_to_me=M_KK_GeV * 1e9 / (me * c**2 / eV),
                         coupling_if_R_eq_Rmin=alpha_at_Rmin)
print("KK  R >= hbar/(m_e c) = %.4e m to keep the n=1 cost below m_e c^2" % R_min)
print("KK  R needed for alpha = 1/137: 2 l_P / sqrt(alpha) = %.3e m  -> first charged mode %.2e GeV (%.1e x m_e)"
      % (R_KK, M_KK_GeV, out["kk_numbers"]["ratio_to_me"]))
print("KK  with R = hbar/(m_e c) the coupling would be 4 l_P^2/R^2 = %.2e (instead of 7.3e-3)" % alpha_at_Rmin)
out["time_s"] = time.time() - t0
json.dump(out, open("out/kk_charge_mass.json", "w"), indent=1)
print("time %.1f s" % out["time_s"])
