"""zz_free.py - T2: the FREE ELECTRON in the zig-zag rule (spin-1/2 Weyl walkers + mass flip theta_m).

Units hbar = c = m_e = 1; lattice a = theta_m lambda_C (the flip angle IS the electron mass: input).
Tests (tolerances in TOL, fixed before the first run):
  F1 dispersion: two-tick Bloch quasi-energies vs E = sqrt(p^2 c^2 + m^2 c^4).
     rest energy exact (1e-12); relative error <= 1 % for p <= 1 mc at theta = 0.1; error ~ theta^2
     (ratio theta 0.2 -> 0.1 in [3, 5]). Kinetic mass at small p reported (expected (1 + theta^2/3) m).
  F2 Zitterbewegung 1-D (same packet as lab/dirac_free.py: sigma = 5, spinor (1,1)/sqrt2 in the Dirac
     representation = a purely RIGHT-handed packet here): frequency and amplitude vs the continuum
     numbers (omega = 2.0095, A = 0.4905 from dirac_free.json): |d omega|/omega < 1 %, |dA|/A < 2 % at
     theta = 0.1. Positive-energy (walk eigenvectors) packet: no ZBW, residual < 1e-6.
  F3 Zitterbewegung 3-D: spinor (L, R) = ((0, i), (1, 0))/sqrt2 at rest -> <r> runs on a circle of
     radius ~ hbar/2mc at 2mc^2/hbar. Walk vs exact continuum Dirac (same packet, spectral): amplitude
     and frequency within 2 %.
  F4 spin 1/2: (a) Berry phase of the R walker band around loops of directions: great circle -> pi
     (exp = -1: a full turn gives a minus sign), cone 60 deg -> pi/2; at lambda = 64a, tolerance 1e-2 rad.
     (b) the lattice rotation by 90 deg about z commutes with 40 ticks of the walk only if the spinor is
     turned by exp(-i chi Sigma_z / 2) with chi = 90 deg (half-angle law -> four quarter turns = -1):
     argmin over chi (15 deg grid) = 90 deg and d(90)/d(180) < 0.1.
  F5 magnetic moment from the Dirac (Gordon) current of a positive-energy spin-up packet at rest:
     g = 2 mu / S, mu = (1/2) Int (x j_y - y j_x), S = <Sigma_z>/2. Walk vs continuum |dg| < 2e-3;
     |g - 2| < 1 % at sigma = 8.
  F6 surfer (1-D, mixed packet, 2000 surfers, RK4 over two ticks): |v| <= c (1e-12) at every stage;
     equivariance: each surfer keeps its quantile of |psi|^2 within 1e-2.
Run: timeout 590 python lab/zz_free.py
"""
import json, os, time
import numpy as np
from scipy.optimize import curve_fit
from zz_core import (Walk, bloch_batch, quasi_energies, gaussian3, dirac_velocity, PAULI, I2, SEED,
                     LAMBDA_C, T_UNIT)

TOL = dict(rest=1e-12, disp=0.01, disp_ratio=(3.0, 5.0), zbw_w=0.01, zbw_A=0.02, zbw_pos=1e-6,
           zbw3=0.02, berry=1e-2, rot_ratio=0.1, g_walk=2e-3, g2=0.01, vmax=1 + 1e-12, quant=1e-2)
OUT = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUT, exist_ok=True)
res, verdict = {}, []
t_start = time.time()
rng = np.random.default_rng(SEED)
DIRS = {"axis": np.array([1.0, 0, 0]), "face": np.array([1.0, 1, 0]) / np.sqrt(2), "body": np.array([1.0, 1, 1]) / np.sqrt(3)}
BETA = np.block([[0 * I2, I2], [I2, 0 * I2]])
SIGZ4 = np.diag([1, -1, 1, -1]).astype(complex)


def check(name, ok, msg):
    verdict.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def walk_E(p_vec, theta):
    """positive quasi-energies (units m c^2) of the two-tick walk at momentum p (units m c)"""
    om, v = quasi_energies(np.asarray(p_vec) * theta, theta, 0.5)
    return om / theta, v


# ---------------- F1 dispersion ----------------
print("=== F1 Dirac dispersion of the walk ===", flush=True)
res["F1"] = {}
ps = [0.0, 0.1, 0.5, 1.0, 2.0]
table = {}
for th in (0.4, 0.2, 0.1):
    for dn, d in DIRS.items():
        row = []
        for p in ps:
            E, _ = walk_E(p * d, th)
            Ep = np.sort(E)[-2:]                     # two positive (spin up/down)
            row.append(float(np.max(np.abs(Ep / np.sqrt(1 + p * p) - 1))))
        table[f"{th} {dn}"] = row
        print(f"  theta={th:4.2f} {dn:5s} " + "  ".join(f"p={p}: {e:.2e}" for p, e in zip(ps, row)), flush=True)
res["F1"]["p"] = ps
res["F1"]["rel_err"] = table
rest = max(table[f"{th} {dn}"][0] for th in (0.4, 0.2, 0.1) for dn in DIRS)
check("F1 rest energy = m c^2 exactly", rest < TOL["rest"], f"max |E(0)/mc^2 - 1| = {rest:.1e}")
e01 = max(max(table[f"0.1 {dn}"][:4]) for dn in DIRS)
e02 = max(max(table[f"0.2 {dn}"][:4]) for dn in DIRS)
check("F1 dispersion error, theta = 0.1, p <= mc", e01 <= TOL["disp"], f"{100 * e01:.3f} %")
check("F1 error ~ theta^2", TOL["disp_ratio"][0] <= e02 / e01 <= TOL["disp_ratio"][1], f"ratio 0.2/0.1 = {e02 / e01:.2f}")
meff = {}
for th in (0.4, 0.2, 0.1):
    p = 1e-3
    E, _ = walk_E(p * DIRS["axis"], th)
    m_eff = p * p / (2 * (np.max(E) - 1))
    meff[th] = m_eff
    print(f"  kinetic mass at theta={th}: {m_eff:.6f} m  (1 + theta^2/3 = {1 + th * th / 3:.6f})", flush=True)
res["F1"]["kinetic_mass"] = meff
# largest group velocity anywhere in the Brillouin zone (information)
K = rng.uniform(-np.pi, np.pi, (20000, 3))
h = 1e-5
vmax_g = 0.0
for th in (0.2,):
    w0 = np.sort(-np.angle(np.linalg.eigvals(bloch_batch(K, th))) / 2, axis=1)
    gs = []
    for q in range(3):
        dK = np.zeros(3); dK[q] = h
        wp = np.sort(-np.angle(np.linalg.eigvals(bloch_batch(K + dK, th))) / 2, axis=1)
        gs.append((wp - w0) / h)
    g = np.sqrt(sum(x**2 for x in gs))
    ok = np.abs(g) < 10       # drop crossings / branch jumps
    vmax_g = float(np.max(g[ok]))
    p95 = float(np.percentile(g[ok], 99))
print(f"  largest group velocity of the walk bands over the whole zone (theta=0.2): {vmax_g:.4f} c (99th percentile {p95:.4f})", flush=True)
res["F1"]["max_group_velocity_zone"] = vmax_g


# ---------------- F2 1-D Zitterbewegung ----------------
def zbw_fit(t, x, tmax):
    sel = t <= tmax
    f = lambda t, c0, c1, A, w, ph: c0 + c1 * t + A * np.sin(w * t + ph)
    best = None
    for ph0 in (0, 1.5, 3, 4.5):
        try:
            p, _ = curve_fit(f, t[sel], x[sel], p0=[x[0], 0.0, 0.5, 2.0, ph0], maxfev=20000)
        except RuntimeError:
            continue
        r = np.sum((f(t[sel], *p) - x[sel])**2)
        if best is None or r < best[1]:
            best = (p, r)
    p = best[0]
    return abs(p[2]), abs(p[3]), p


print("=== F2 Zitterbewegung, 1-D (sigma = 5, purely right-handed packet at rest) ===", flush=True)
ref = dict(omega=2.0095, amp=0.4905)
try:
    dj = json.load(open(os.path.join(OUT, "dirac_free.json")))
    ref = dict(omega=dj["1d_zbw"]["omega"], amp=dj["1d_zbw"]["amp"])
except Exception:
    pass
res["F2"] = dict(continuum=ref)


def walk1d_packet(th, sigma, positive=False):
    a = th
    n = int(round(400.0 / a))
    W = Walk("dirac1", n, a, dim=1, theta=th)
    x = W.x1
    g = np.exp(-x**2 / (4 * sigma**2)).astype(complex)
    psi = np.array([g, 0 * g])
    if positive:
        k = 2 * np.pi * np.fft.fftfreq(n)                     # k a
        pk = np.fft.fft(psi, axis=1)
        out = np.zeros_like(pk)
        for i in range(n):
            M = np.array([[np.cos(th), -1j * np.sin(th)], [-1j * np.sin(th), np.cos(th)]]) @ np.diag([np.exp(-1j * k[i]), np.exp(1j * k[i])])
            w, v = np.linalg.eig(M)
            j = int(np.argmax(-np.angle(w)))                   # positive frequency
            u = v[:, j] / np.linalg.norm(v[:, j])
            out[:, i] = u * np.vdot(u, pk[:, i])
        psi = np.fft.ifft(out, axis=1)
    psi /= np.sqrt(np.sum(np.abs(psi)**2))
    return W, psi


for th in (0.2, 0.1, 0.05):
    rows = {}
    for positive in (False, True):
        W, psi = walk1d_packet(th, 5.0, positive)
        nt = int(round(40.0 / th))
        ts, xs = [], []
        for t in range(nt + 1):
            rho = np.sum(np.abs(psi)**2, axis=0)
            ts.append(t * th); xs.append(np.sum(rho * W.x1) / rho.sum())
            psi = W.tick(psi)
        ts, xs = np.array(ts), np.array(xs)
        if not positive:
            A, w, p = zbw_fit(ts, xs, 4 * np.pi)
            rows.update(amp=A, omega=w)
        else:
            resid = float(np.max(np.abs(xs - np.polyval(np.polyfit(ts, xs, 2), ts))))
            rows.update(positive_residual=resid)
    res["F2"][f"theta {th}"] = rows
    print(f"  theta={th}: omega = {rows['omega']:.5f} mc^2/hbar, A = {rows['amp']:.5f} lambda_C, positive-only residual {rows['positive_residual']:.1e}", flush=True)
r1 = res["F2"]["theta 0.1"]
check("F2 ZBW frequency (theta = 0.1)", abs(r1["omega"] / ref["omega"] - 1) < TOL["zbw_w"],
      f"{r1['omega']:.5f} vs continuum {ref['omega']:.5f} ({r1['omega'] / T_UNIT:.4e} rad/s; 2mc^2/hbar = {2 / T_UNIT:.4e})")
check("F2 ZBW amplitude (theta = 0.1)", abs(r1["amp"] / ref["amp"] - 1) < TOL["zbw_A"],
      f"{r1['amp']:.5f} vs continuum {ref['amp']:.5f} ({r1['amp'] * LAMBDA_C:.4e} m; hbar/2mc = {0.5 * LAMBDA_C:.4e} m)")
check("F2 no ZBW for positive energy (walk eigenvectors)", r1["positive_residual"] < TOL["zbw_pos"], f"{r1['positive_residual']:.1e}")


# ---------------- F3 3-D Zitterbewegung circle ----------------
print("=== F3 Zitterbewegung, 3-D: circle of radius hbar/2mc ===", flush=True)
th3, n3, sig3 = 0.25, 96, 3.0
W3 = Walk("dirac", n3, th3, theta=th3)
spin3 = np.array([0, 1j, 1, 0]) / np.sqrt(2)
psi0 = gaussian3(W3, sig3, spin3)
x = W3.x1
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")


def means(psi):
    rho = np.sum(np.abs(psi)**2, axis=0)
    s = rho.sum()
    return np.array([np.sum(rho * X), np.sum(rho * Y), np.sum(rho * Z)]) / s


# continuum: exact spectral Dirac evolution (chiral representation) on the same grid
kk = 2 * np.pi * np.fft.fftfreq(n3, d=th3)
KX, KY, KZ = np.meshgrid(kk, kk, kk, indexing="ij")
Ek = np.sqrt(KX**2 + KY**2 + KZ**2 + 1)


def H_apply(phi):
    L, R = phi[0:2], phi[2:4]
    def sdotk(s):
        a, b = s
        return np.array([KZ * a + (KX - 1j * KY) * b, (KX + 1j * KY) * a - KZ * b])
    return np.concatenate([-sdotk(L) + R, sdotk(R) + L])


phi0 = np.fft.fftn(psi0, axes=(1, 2, 3))
nt3 = int(round(4 * np.pi / th3))
tw, mw, mc = [], [], []
psi = psi0.copy()
for t in range(nt3 + 1):
    tt = t * th3
    tw.append(tt); mw.append(means(psi))
    if t % 2 == 0:
        ph = np.cos(Ek * tt) * phi0 - 1j * np.sin(Ek * tt) * H_apply(phi0) / Ek
        mc.append(means(np.fft.ifftn(ph, axes=(1, 2, 3))))
    else:
        mc.append(mc[-1])
    psi = W3.tick(psi)
tw, mw, mc = np.array(tw), np.array(mw), np.array(mc)
r3 = {}
for lab_, arr, step in (("walk", mw, 1), ("continuum", mc, 2)):
    out = {}
    for ax, nm in ((0, "x"), (2, "z")):
        A, w, p = zbw_fit(tw[::step], arr[::step, ax], 4 * np.pi)
        out[nm] = dict(amp=A, omega=w)
    r3[lab_] = out
    print(f"  {lab_:9s}: x: A = {out['x']['amp']:.4f}, omega = {out['x']['omega']:.4f};  z: A = {out['z']['amp']:.4f}, omega = {out['z']['omega']:.4f}"
          f"   (|<y>| max {np.max(np.abs(arr[:, 1])):.1e})", flush=True)
res["F3"] = r3
okA = all(abs(r3["walk"][c]["amp"] / r3["continuum"][c]["amp"] - 1) < TOL["zbw3"] for c in "xz")
okW = all(abs(r3["walk"][c]["omega"] / r3["continuum"][c]["omega"] - 1) < TOL["zbw3"] for c in "xz")
check("F3 3-D ZBW circle: walk vs continuum", okA and okW,
      f"radius {r3['walk']['x']['amp']:.4f}/{r3['walk']['z']['amp']:.4f} vs {r3['continuum']['x']['amp']:.4f}/{r3['continuum']['z']['amp']:.4f} (hbar/2mc = 0.5 for a wide packet)")
np.save(os.path.join(OUT, "zz_zbw3d_track.npy"), np.column_stack([tw, mw, mc]))


# ---------------- F4 spin 1/2 ----------------
print("=== F4 spin 1/2 ===", flush=True)


def berry_weyl(kmag, cone_deg, npts=400):
    vs = []
    th = np.radians(cone_deg)
    for p in np.linspace(0, 2 * np.pi, npts, endpoint=False):
        kv = kmag * np.array([np.sin(th) * np.cos(p), np.sin(th) * np.sin(p), np.cos(th)])
        om, v = quasi_energies(kv, 0.0, 0.5)
        # R block, positive frequency (moves along +k)
        R = v[2:, :]
        j = int(np.argmax(om * (np.linalg.norm(R, axis=0) > 0.5)))
        u = v[:, j] / np.linalg.norm(v[:, j])
        vs.append(u)
    prod = 1.0 + 0j
    for i in range(npts):
        prod *= np.vdot(vs[i], vs[(i + 1) % npts])
    return -np.angle(prod)


res["F4"] = {"berry": {}}
for cone, expect in ((90, np.pi), (60, np.pi / 2)):
    for lam in (16, 32, 64):
        g = berry_weyl(2 * np.pi / lam, cone)
        # the sign of gamma depends on the orientation of the loop and on the helicity: compare |gamma|
        dev = float(abs(np.angle(np.exp(1j * g))) - expect)
        res["F4"]["berry"][f"cone {cone} lambda {lam}a"] = dict(gamma=g, dev=dev)
        print(f"  Weyl R walker, cone {cone} deg, lambda = {lam}a: gamma = {g:+.5f} (expected {expect:.5f}), exp(i gamma) = {np.exp(1j * g):.4f}", flush=True)
        if lam == 64:
            check(f"F4a Berry phase spin 1/2, cone {cone}", abs(dev) < TOL["berry"], f"gamma - expected = {dev:+.2e} rad")

nR, thR = 40, 0.2
WR = Walk("dirac", nR, 1.0, theta=thR)
spinR = rng.normal(size=4) + 1j * rng.normal(size=4)
psiR = gaussian3(WR, 4.0, spinR, k0=(0.25, -0.15, 0.1))


def rot90(psi, chi):
    """lattice rotation by +90 deg about z: psi'(x, y) = U psi(y, -x), U = exp(-i chi Sigma_z / 2)"""
    U = np.diag(np.exp(-0.5j * chi * np.array([1, -1, 1, -1])))
    out = np.rot90(psi, k=1, axes=(1, 2))
    return np.tensordot(U, out, axes=(1, 0))


def evolve(psi, nt=40):
    W = Walk("dirac", nR, 1.0, theta=thR)
    for _ in range(nt):
        psi = W.tick(psi)
    return psi


base = evolve(psiR)
chis = np.radians(np.arange(0, 361, 15))
ds = []
for chi in chis:
    d = np.linalg.norm(evolve(rot90(psiR, chi)) - rot90(base, chi)) / np.linalg.norm(psiR)
    ds.append(d)
ds = np.array(ds)
res["F4"]["rotation_scan_deg"] = np.degrees(chis).tolist()
res["F4"]["rotation_mismatch"] = ds.tolist()
print("  mismatch after 40 ticks vs spinor angle chi: " + ", ".join(f"{int(round(np.degrees(c)))}: {d:.3f}" for c, d in zip(chis, ds)), flush=True)
i90, i180 = 6, 12
check("F4b quarter turn needs the half-angle spinor turn", int(np.argmin(ds[:13])) == i90 and ds[i90] / ds[i180] < TOL["rot_ratio"],
      f"argmin chi = {np.degrees(chis[np.argmin(ds[:13])]):.0f} deg, d(90)/d(180) = {ds[i90] / ds[i180]:.3f}; "
      f"=> four quarter turns: exp(-i pi Sigma_z) = -1, eight = +1 (4 pi periodicity)")


# ---------------- F5 g-factor ----------------
print("=== F5 magnetic moment from the Dirac current (g) ===", flush=True)


def g_factor(sigma, n, a, walk=True):
    W = Walk("dirac", n, a, theta=a)
    psi = gaussian3(W, sigma, [1, 0, 1, 0])
    pk = np.fft.fftn(psi, axes=(1, 2, 3)).reshape(4, -1)
    kk = 2 * np.pi * np.fft.fftfreq(n, d=a)
    KX, KY, KZ = [q.ravel() for q in np.meshgrid(kk, kk, kk, indexing="ij")]
    wgt = np.sum(np.abs(pk)**2, axis=0)
    sel = np.nonzero(wgt > 1e-14 * wgt.max())[0]
    out = np.zeros_like(pk)
    if walk:
        Wb = bloch_batch(np.column_stack([KX[sel], KY[sel], KZ[sel]]) * a, a, 0.5)
        w, v = np.linalg.eig(Wb)
        om = -np.angle(w)
        # positive-frequency subspace (2-fold degenerate: spin): eig does not return orthonormal vectors
        # inside a degenerate pair, so orthonormalize them (QR) before projecting (bug of the first two runs)
        o = np.argsort(-om, axis=1)[:, :2]
        Vp = np.take_along_axis(v, o[:, None, :], axis=2)
        Q, _ = np.linalg.qr(Vp)
        c = np.einsum("nab,na->nb", Q.conj(), pk[:, sel].T)
        out[:, sel] = np.einsum("nab,nb->na", Q, c).T
    else:
        E = np.sqrt(KX[sel]**2 + KY[sel]**2 + KZ[sel]**2 + 1)
        ph = pk[:, sel]
        L, R = ph[0:2], ph[2:4]
        def sdotk(s):
            a_, b_ = s
            return np.array([KZ[sel] * a_ + (KX[sel] - 1j * KY[sel]) * b_, (KX[sel] + 1j * KY[sel]) * a_ - KZ[sel] * b_])
        Hph = np.concatenate([-sdotk(L) + R, sdotk(R) + L])
        out[:, sel] = 0.5 * (ph + Hph / E)
    psi = np.fft.ifftn(out.reshape(4, n, n, n), axes=(1, 2, 3))
    psi /= np.sqrt(np.sum(np.abs(psi)**2))
    v, rho = dirac_velocity(psi)                 # v = j / rho, so j = v rho
    x = W.x1
    X, Y, _ = np.meshgrid(x, x, x, indexing="ij")
    mu = 0.5 * np.sum(X * v[1] * rho - Y * v[0] * rho)
    S = 0.5 * np.real(np.sum(np.conj(psi) * np.tensordot(SIGZ4, psi, axes=(1, 0))))
    return 2 * mu / S, S, mu


res["F5"] = {}
for sigma, n, a in ((4.0, 64, 0.5), (8.0, 128, 0.5)):
    gw, Sw, muw = g_factor(sigma, n, a, True)
    gc, Sc, muc = g_factor(sigma, n, a, False)
    res["F5"][f"sigma {sigma}"] = dict(g_walk=gw, g_continuum=gc, S=Sw, mu=muw)
    print(f"  sigma = {sigma}: walk g = {gw:.5f} (S = {Sw:.5f}, mu = {muw:.5f}); continuum g = {gc:.5f}", flush=True)
r8 = res["F5"]["sigma 8.0"]
check("F5 g walk vs continuum", abs(r8["g_walk"] - r8["g_continuum"]) < TOL["g_walk"], f"{r8['g_walk']:.5f} vs {r8['g_continuum']:.5f}")
check("F5 g = 2 (sigma = 8)", abs(r8["g_walk"] / 2 - 1) < TOL["g2"], f"g = {r8['g_walk']:.5f}")


# ---------------- F6 surfers 1-D ----------------
print("=== F6 surfers, 1-D mixed packet ===", flush=True)
th6 = 0.1
W6, psi6 = walk1d_packet(th6, 5.0, False)
x6 = W6.x1
N6 = 2000
rho0 = np.sum(np.abs(psi6)**2, axis=0)
cdf0 = np.cumsum(rho0) - 0.5 * rho0
cdf0 /= rho0.sum()
pos = np.interp(rng.random(N6), cdf0, x6)


def vel(psi, p):
    R = np.interp(p, x6, psi[0].real) + 1j * np.interp(p, x6, psi[0].imag)
    L = np.interp(p, x6, psi[1].real) + 1j * np.interp(p, x6, psi[1].imag)
    return (np.abs(R)**2 - np.abs(L)**2) / (np.abs(R)**2 + np.abs(L)**2)


def quant(psi, p):
    rho = np.sum(np.abs(psi)**2, axis=0)
    c = np.cumsum(rho) - 0.5 * rho
    return np.interp(p, x6, c / rho.sum())


q0 = quant(psi6, pos)
vmax, qerr = 0.0, 0.0
nsteps = int(round(40.0 / (2 * th6)))
fields = [psi6]
traj = [pos.copy()]
for s in range(nsteps):
    p1 = W6.tick(fields[-1]); p2 = W6.tick(p1)
    dt = 2 * th6
    k1 = vel(fields[-1], pos); k2 = vel(p1, pos + 0.5 * dt * k1); k3 = vel(p1, pos + 0.5 * dt * k2); k4 = vel(p2, pos + dt * k3)
    vmax = max(vmax, *(float(np.max(np.abs(k))) for k in (k1, k2, k3, k4)))
    pos = pos + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    fields = [p2]
    qerr = max(qerr, float(np.max(np.abs(quant(p2, pos) - q0))))
    if s % 5 == 0:
        traj.append(pos.copy())
res["F6"] = dict(vmax=vmax, quantile_err=qerr)
check("F6 surfer |v| <= c", vmax <= TOL["vmax"], f"max |v| = {vmax:.15f}")
check("F6 equivariance (quantile kept)", qerr < TOL["quant"], f"max |F(x_i(t)) - F(x_i(0))| = {qerr:.2e}")
np.save(os.path.join(OUT, "zz_surf1d.npy"), np.array(traj)[:, :40])

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for th in (0.4, 0.2, 0.1):
        ax[0].semilogy(ps[1:], table[f"{th} body"][1:], "o-", label=f"theta = {th}")
    ax[0].set_xlabel("p / mc (body diagonal)"); ax[0].set_ylabel("|E_walk/E - 1|"); ax[0].legend(fontsize=7)
    ax[0].set_title("electron walk: dispersion error")
    ax[1].plot(mw[:, 0], mw[:, 2], "r-", label="walk")
    ax[1].plot(mc[::2, 0], mc[::2, 2], "k:", label="Dirac (continuum)")
    ax[1].set_aspect("equal"); ax[1].set_xlabel("<x> [hbar/mc]"); ax[1].set_ylabel("<z>"); ax[1].legend(fontsize=7)
    ax[1].set_title("3-D Zitterbewegung circle")
    ax[2].plot(np.degrees(chis), ds, "o-")
    ax[2].axvline(90, color="k", ls=":")
    ax[2].set_xlabel("spinor turn chi [deg] with a 90 deg lattice turn"); ax[2].set_ylabel("mismatch after 40 ticks")
    ax[2].set_title("spin 1/2: half-angle law")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "zz_free.png"), dpi=90)
except Exception as ex:
    print("figure skipped:", ex)

res["verdict"] = verdict
res["seconds"] = time.time() - t_start
with open(os.path.join(OUT, "zz_free.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
npass = sum(v for _, v in verdict)
print(f"=== {npass} PASS / {len(verdict) - npass} FAIL  ({time.time() - t_start:.0f} s)", flush=True)
