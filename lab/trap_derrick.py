"""trap_derrick.py - Parte B2c: forme localizzate STABILI richiedono topologia + termini extra?

(A) Teorema di Derrick in 3D per un campo scalare non lineare (senza topologia):
      phi_tt = lap phi - phi + phi^3 ;  soluzione statica Q(r) (stato fondamentale)
    - energia lungo la famiglia di dilatazioni Q(mu x): massimo a mu = 1 (instabile)
    - evoluzione nel tempo da (1 +/- 1e-3) Q: collasso o dispersione
(B) Skyrmione 3D (pioni, Skyrme 1961): profilo a riccio (hedgehog) minimizzato, numero barionico B,
    viriale E2 = E4, stabilita' in scala; poi il calcolo di Adkins-Nappi-Witten 1983: si FITTANO
    F_pi ed e su M_N e M_Delta e si PREDICONO raggio, momenti magnetici, g_A (confronto con misure).
    Controllo: senza termine di Skyrme il riccio collassa (Derrick).
(C) Baby-skyrmione 2D (campo unitario n in S^2, Piette-Schroers-Zakrzewski 1995): minimizzazione
    dell'energia sul campo 2D completo (nessuna simmetria imposta) partendo da una forma deformata;
    carica topologica Q, energia, confronto con il profilo radiale; senza termine di Skyrme collassa.

uso: timeout 590 python lab/trap_derrick.py
"""
import json, os, sys, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

OUT = os.path.join(os.path.dirname(__file__), "out")
T0 = time.time()
rng = np.random.default_rng(11)
res = {}


def log(*a):
    print(f"[{time.time()-T0:6.1f}s]", *a, flush=True)


# ====================================================================== (A) Derrick 3D
log("(A) Derrick: solitone scalare 3D senza topologia")


def shoot_Q(Q0, rmax=25):
    f = lambda r, y: [y[1], -2 * y[1] / r + y[0] - y[0] ** 3]
    ev = lambda r, y: y[0]
    ev.terminal = True
    ev2 = lambda r, y: y[1]
    ev2.terminal = True
    return solve_ivp(f, [1e-4, rmax], [Q0, 1e-4 / 3 * (Q0 - Q0**3)], rtol=1e-12, atol=1e-14, events=[ev, ev2], dense_output=True)


lo, hi = 3.0, 6.0
for _ in range(90):
    mid = 0.5 * (lo + hi)
    s = shoot_Q(mid)
    if s.t_events[0].size:
        hi = mid
    else:
        lo = mid
Q0 = lo
r = np.linspace(1e-4, 30, 30001)
sol = shoot_Q(Q0, 30)
Qr = np.zeros_like(r)
mask = r <= sol.t[-1]
Qr[mask] = sol.sol(r[mask])[0]
cut = np.argmin(np.where(r < sol.t[-1], np.abs(Qr), 9))  # punto dove la bisezione si stacca
Qr[r > min(sol.t[-1], 14.0)] = 0
Qp = np.gradient(Qr, r)
T = 4 * np.pi * np.trapezoid(0.5 * Qp**2 * r**2, r)
V2 = 4 * np.pi * np.trapezoid(0.5 * Qr**2 * r**2, r)
V4 = 4 * np.pi * np.trapezoid(0.25 * Qr**4 * r**2, r)
mus = np.array([0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2])
Emu = T / mus + (V2 - V4) / mus**3
log(f"  Q(0) = {Q0:.5f} (letteratura 4.3374);  T = {T:.3f}, V2 = {V2:.3f}, V4 = {V4:.3f};  "
    f"Derrick/Pohozaev: T + 3(V2-V4) = {T + 3*(V2-V4):.2e} (deve essere 0)")
log(f"  E(mu) per mu = {mus.tolist()}: {np.round(Emu, 4).tolist()}  -> massimo a mu=1: instabile")
dE2 = 2 * T + 12 * (V2 - V4)
log(f"  d2E/dmu2 a mu=1 = {dE2:.3f} (< 0)")
res["derrick_3d"] = dict(Q0=Q0, T=T, V2=V2, V4=V4, mu=mus.tolist(), E=Emu.tolist(), d2E=dE2)

# evoluzione radiale: u = r phi,  u_tt = u_rr - u + u^3 / r^2 ; spugna oltre r = 40
dr, Rm = 0.02, 60.0
rg = np.arange(1, int(Rm / dr) + 1) * dr
Qg = np.interp(rg, r, Qr, right=0.0)
sig = np.where(rg > 40, 0.5 * ((rg - 40) / 20) ** 2, 0.0)
dyn = {}
for fac in (0.999, 1.0, 1.001):
    u = fac * rg * Qg
    v = np.zeros_like(u)
    dt = 0.01
    tr, pk = [], []
    status = "resta (entro il tempo simulato)"
    for i in range(int(60 / dt) + 1):
        if i % 200 == 0:
            phi0 = u[0] / rg[0]
            tr.append(i * dt)
            pk.append(float(phi0))
        lap = np.empty_like(u)
        lap[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / dr**2
        lap[0] = (u[1] - 2 * u[0]) / dr**2          # u(0) = 0
        lap[-1] = (-2 * u[-1] + u[-2]) / dr**2
        a = lap - u + u**3 / rg**2 - sig * v
        v += dt * a
        u += dt * v
        if abs(u[0] / rg[0]) > 50:
            status = f"COLLASSO (phi(0) > 50) a t = {i*dt:.2f}"
            tr.append(i * dt); pk.append(float(u[0] / rg[0]))
            break
        if i * dt > 5 and abs(u[0] / rg[0]) < 0.05 * Q0:
            status = f"SI DISPERDE (phi(0) < 5% di Q0) a t = {i*dt:.2f}"
            tr.append(i * dt); pk.append(float(u[0] / rg[0]))
            break
    dyn[str(fac)] = dict(t=tr, phi0=pk, status=status)
    log(f"  partenza {fac} x Q: {status}")
res["derrick_3d"]["dynamics"] = dyn

# ====================================================================== (B) Skyrmione 3D
log("(B) Skyrmione 3D a riccio (unita' x = e F_pi r)")
# Energia (ANW 1983): M = 4 pi (F_pi/e) Int dx [x^2 F'^2/8 + sin^2F/4 + sin^4F/(2x^2) + sin^2F F'^2]
# Eulero-Lagrange: (x^2/4 + 2 s^2) F'' + x F'/2 + sin2F F'^2 - sin2F/4 - s^2 sin2F/x^2 = 0
# Soluzione per TIRO: F(0)=pi, F ~ pi - a x; si cerca a tale che F -> 0 (coda F ~ A/x^2, pioni senza massa).


def sk_rhs(x, y):
    F, Fp = y
    s, s2F = np.sin(F), np.sin(2 * F)
    return [Fp, -(x * Fp / 2 + s2F * Fp**2 - s2F / 4 - s * s * s2F / x**2) / (x**2 / 4 + 2 * s * s)]


def sk_shoot(a, xmax=30.0):
    x0 = 1e-3
    evz = lambda x, y: y[0]
    evz.terminal = True
    evu = lambda x, y: y[1]
    evu.terminal = True
    return solve_ivp(sk_rhs, [x0, xmax], [np.pi - a * x0, -a], rtol=1e-12, atol=1e-14,
                     events=[evz, evu], dense_output=True)


alo, ahi = 1.0, 3.0
for _ in range(60):
    am = 0.5 * (alo + ahi)
    s_ = sk_shoot(am)
    if s_.t_events[0].size:     # F attraversa 0: pendenza troppo grande
        ahi = am
    else:
        alo = am
a_sk = 0.5 * (alo + ahi)
s_ = sk_shoot(a_sk)
xend = s_.t[-1]
# raccordo con la coda A/x^2 dove F e' ancora affidabile (F ~ 0.05)
xx = np.linspace(1e-3, xend, 200001)
FF = s_.sol(xx)[0]
im = int(np.argmax(FF < 0.05)) if (FF < 0.05).any() else len(xx) // 2
xmatch = xx[im]
A_tail = FF[im] * xmatch**2
xt = np.geomspace(xmatch, 5e3, 20001)[1:]
x_all = np.concatenate([xx[:im + 1], xt])
F_all = np.concatenate([FF[:im + 1], A_tail / xt**2])
Fp_all = np.gradient(F_all, x_all)
s2 = np.sin(F_all) ** 2
I = lambda f: np.trapezoid(f, x_all)
E2 = 4 * np.pi * I(x_all**2 * Fp_all**2 / 8 + s2 / 4)
E4 = 4 * np.pi * I(s2**2 / (2 * x_all**2) + s2 * Fp_all**2)
M_coef = E2 + E4
B = -(2 / np.pi) * I(s2 * Fp_all)
Lam = I(x_all**2 * s2 * (1 + 4 * (Fp_all**2 + s2 / x_all**2)))
mus = np.array([0.8, 0.9, 1.0, 1.1, 1.25])
Esc = E2 / mus + E4 * mus
log(f"  pendenza al centro a = {a_sk:.6f} (in unita' e f_pi r: {2*a_sk:.4f}; letteratura 2.007), coda F = {A_tail:.3f}/x^2 da x = {xmatch:.2f}")
log(f"  M = {M_coef:.3f} F_pi/e (ANW 36.5);  M / limite di Bogomolny (12 pi^2 F_pi/(4e) = {3*np.pi**2:.3f} F_pi/e)"
    f" = {M_coef/(3*np.pi**2):.4f} (letteratura 1.232)")
log(f"  E2 = {E2:.3f}, E4 = {E4:.3f} (viriale di Derrick: uguali), B = {B:.6f}, Lambda = {Lam:.3f} (ANW 50.9)")
log(f"  energia sotto dilatazione mu = {mus.tolist()}: {np.round(Esc, 3).tolist()} -> MINIMO a mu=1: stabile")

# --- controllo: SENZA termine di Skyrme il riccio collassa (minimizzazione diretta su griglia)
N = 1200
xs = np.linspace(0, 12.0, N + 1)
x = xs[1:-1]
h = xs[1] - xs[0]
xm = 0.5 * (xs[1:] + xs[:-1])


def E_sigma(Fi):
    F = np.concatenate([[np.pi], Fi, [0.0]])
    Fp = np.diff(F) / h
    Fm = 0.5 * (F[1:] + F[:-1])
    s = np.sin(Fm); c = np.cos(Fm)
    e = (1 / 8) * (Fp**2 * xm**2 + 2 * s * s)
    dp = (1 / 4) * Fp * xm**2
    dm = (1 / 8) * 4 * s * c
    g = np.zeros_like(F)
    g[1:] += dp / h + 0.5 * dm
    g[:-1] += -dp / h + 0.5 * dm
    return 4 * np.pi * h * e.sum(), 4 * np.pi * h * g[1:-1]


F0 = 4 * np.arctan(np.exp(-x))
opt0 = minimize(E_sigma, F0, jac=True, method="L-BFGS-B", options=dict(maxiter=4000))
F00 = np.concatenate([[np.pi], opt0.x, [0.0]])
half = xs[np.argmin(np.abs(F00 - np.pi / 2))]
log(f"  SENZA termine di Skyrme: energia {opt0.fun:.4f} (partenza {E_sigma(F0)[0]:.3f}), raggio dove F = pi/2: "
    f"{half:.4f} (partenza {np.log(1/np.tan(np.pi/8)):.3f}) -> il riccio collassa fino al passo di griglia ({h:.3f})")

# --- fit ANW: M_N, M_Delta -> F_pi, e ; poi predizioni a parametri fissati
MN, MD = 938.92, 1232.0
inv_lam = (MD - MN) * 2 / 3          # 1/lambda,  M_D - M_N = 3/(2 lambda)
M = MN - 3 / 8 * inv_lam
e = (M_coef * (2 * np.pi / 3) * Lam / (M / inv_lam)) ** 0.25
Fpi = M * e / M_coef
hbarc = 197.3269804
r2_I0 = -(2 / np.pi) * I(x_all**2 * s2 * Fp_all) / (e * Fpi) ** 2   # MeV^-2
r_I0_fm = np.sqrt(r2_I0) * hbarc
mu_p = r2_I0 * MN * (MD - MN) / 9 + MN / (2 * (MD - MN))
mu_n = r2_I0 * MN * (MD - MN) / 9 - MN / (2 * (MD - MN))
s2F = np.sin(2 * F_all)
Dint = I(x_all**2 * Fp_all + x_all * s2F + 4 * (x_all * s2F * Fp_all**2 + 2 * s2 * Fp_all + s2 * s2F / x_all))
gA = -np.pi / (3 * e**2) * Dint
pred = dict(F_pi_MeV=Fpi, e=e, M_classical=M, r_I0_fm=r_I0_fm, mu_p=mu_p, mu_n=mu_n,
            mu_ratio=abs(mu_p / mu_n), g_A=gA)
expt = dict(F_pi_MeV=184.2, r_I0_fm=0.77, mu_p=2.7928, mu_n=-1.9130, mu_ratio=1.46, g_A=1.2754)
log(f"  FIT (2 numeri usati: M_N = {MN}, M_Delta = {MD}): F_pi = {Fpi:.1f} MeV (misura 184 con la convenzione"
    f" F_pi = 2 f_pi; ANW 129), e = {e:.3f} (ANW 5.45)")
log(f"  PREDIZIONI a F_pi, e fissati: raggio isoscalare {r_I0_fm:.3f} fm (ANW 0.59; misura 0.77 = sqrt(r_p^2 + r_n^2)), mu_p = {mu_p:.3f} "
    f"(2.793), mu_n = {mu_n:.3f} (-1.913), |mu_p/mu_n| = {abs(mu_p/mu_n):.3f} (1.46), g_A = {gA:.3f} (ANW 0.61; 1.275)")
res["skyrmion"] = dict(a=a_sk, A_tail=A_tail, M_coef=M_coef, E2=E2, E4=E4, B=B, Lambda=Lam,
                       scale_mu=mus.tolist(), scale_E=Esc.tolist(),
                       no_skyrme_final_E=float(opt0.fun), no_skyrme_halfradius=float(half),
                       fit_inputs=dict(M_N=MN, M_Delta=MD), predictions=pred, experiment=expt,
                       profile_x=xx[:im:400].tolist(), profile_F=FF[:im:400].tolist())

# ====================================================================== (C) baby-skyrmione 2D
log("(C) baby-skyrmione 2D sul campo completo (nessuna simmetria imposta)")
n2, L2 = 96, 16.0
hh = L2 / n2
g1 = (np.arange(n2) - n2 / 2 + 0.5) * hh
XX, YY = np.meshgrid(g1, g1, indexing="ij")


# differenze IN AVANTI: le differenze centrali hanno un modo a scacchiera "invisibile" (energia zero)
# attraverso il quale il campo si srotola sul reticolo; le differenze in avanti non lo hanno.
def D1(f):
    return (np.roll(f, -1, axis=1) - f) / hh


def D2(f):
    return (np.roll(f, -1, axis=2) - f) / hh


def D1T(g):      # trasposto di D1
    return (np.roll(g, 1, axis=1) - g) / hh


def D2T(g):
    return (np.roll(g, 1, axis=2) - g) / hh


def cross(a, b):
    return np.array([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]])


def baby(n, kap2, mu2):
    a, b = D1(n), D2(n)
    axb = cross(a, b)
    q = np.sum(n * axb, 0)
    E2 = 0.5 * np.sum(a**2 + b**2) * hh**2
    E4 = 0.5 * kap2 * np.sum(q**2) * hh**2
    E0 = mu2 * np.sum(1 - n[2]) * hh**2
    g = D1T(a) + D2T(b)
    g += kap2 * (q * axb + D1T(q * cross(b, n)) + D2T(q * cross(n, a)))
    g[2] += -mu2
    return E2, E4, E0, g * hh**2, q


def charge(n):
    return np.sum(baby(n, 0, 0)[4]) * hh**2 / (4 * np.pi)


def hedgehog(fprof, deform=0.0, noise=0.0):
    Xd, Yd = XX * (1 + deform), YY / (1 + deform)
    rr = np.hypot(Xd, Yd)
    th = np.arctan2(Yd, Xd)
    f = fprof(rr)
    n = np.array([np.sin(f) * np.cos(th), np.sin(f) * np.sin(th), np.cos(f)])
    n += noise * rng.normal(size=n.shape)
    return n / np.sqrt(np.sum(n**2, 0))


def arrested_flow(n, kap2, mu2, steps, dt, every=500):
    """flusso di Newton 'arrestato': n'' = -grad E (tangente), velocita' azzerata se E sale."""
    v = np.zeros_like(n)
    E2, E4, E0, g, _ = baby(n, kap2, mu2)
    Eold = E2 + E4 + E0
    hist = []
    for it in range(steps):
        g = g - n * np.sum(n * g, 0)
        v = v - dt * g / hh**2
        v = v - n * np.sum(n * v, 0)
        n = n + dt * v
        n = n / np.sqrt(np.sum(n**2, 0))
        E2, E4, E0, g, _ = baby(n, kap2, mu2)
        E = E2 + E4 + E0
        if E > Eold:
            v[:] = 0
        Eold = E
        if it % every == 0 or it == steps - 1:
            rho = np.sqrt(np.sum(n[2] < 0) * hh**2 / np.pi)   # raggio della zona n3 < 0 (dove f > pi/2)
            hist.append(dict(it=it, E=float(E), Q=float(charge(n)), size=float(rho)))
    return n, hist


kap2, mu2 = 1.0, 0.1
n0 = hedgehog(lambda rr: np.pi * np.exp(-rr / 2.0), deform=0.35, noise=0.02)
dn = rng.normal(size=n0.shape)
dn -= n0 * np.sum(n0 * dn, 0)
ee = 1e-6
Ef = lambda n: sum(baby(n, kap2, mu2)[:3])
fd = (Ef(n0 + ee * dn) - Ef(n0 - ee * dn)) / (2 * ee)
log(f"  controllo gradiente: analitico {np.sum(baby(n0, kap2, mu2)[3] * dn):.7f}  differenze finite {fd:.7f}")
Q_start, E_start = charge(n0), Ef(n0)
nf, hist = arrested_flow(n0, kap2, mu2, 12000, 0.02 * hh)
E2b, E4b, E0b, _, q = baby(nf, kap2, mu2)
Qf = charge(nf)
w = q**2
Ixx, Iyy, Ixy = np.sum(w * XX**2), np.sum(w * YY**2), np.sum(w * XX * YY)
ev = np.linalg.eigvalsh(np.array([[Ixx, Ixy], [Ixy, Iyy]]))
aniso = ev[1] / ev[0] - 1
w0 = baby(n0, 0, 0)[4] ** 2
ev0 = np.linalg.eigvalsh(np.array([[np.sum(w0 * XX**2), np.sum(w0 * XX * YY)], [np.sum(w0 * XX * YY), np.sum(w0 * YY**2)]]))
log(f"  partenza deformata (ellisse 35% + rumore 2%): Q = {Q_start:.4f}, E = {E_start:.4f}, anisotropia {ev0[1]/ev0[0]-1:.3f}")
for hst in hist[::4] + [hist[-1]]:
    log(f"    iter {hst['it']:6d}: E = {hst['E']:.5f}  Q = {hst['Q']:.4f}  raggio(n3<0) = {hst['size']:.4f}")
log(f"  finale: E = {E2b+E4b+E0b:.4f} = E2 {E2b:.4f} (>= 4pi|Q| = {4*np.pi*abs(Qf):.4f}) + E4 {E4b:.4f} + E0 {E0b:.4f}"
    f" (viriale di Derrick 2D: E4 = E0), Q = {Qf:.4f}, anisotropia residua {aniso:.2e}")
rr1 = np.linspace(0, 16, 3201)
hr = rr1[1] - rr1[0]
rm = 0.5 * (rr1[1:] + rr1[:-1])


def Erad(fi):
    f = np.concatenate([[np.pi], fi, [0.0]])
    fp = np.diff(f) / hr
    fm = 0.5 * (f[1:] + f[:-1])
    s, c = np.sin(fm), np.cos(fm)
    dens = 0.5 * (fp**2 + s * s / rm**2) + 0.5 * kap2 * s * s * fp**2 / rm**2 + mu2 * (1 - c)
    dp = fp + kap2 * s * s * fp / rm**2
    dm = s * c / rm**2 + kap2 * s * c * fp**2 / rm**2 + mu2 * s
    g = np.zeros_like(f)
    g[1:] += (dp / hr + 0.5 * dm) * rm
    g[:-1] += (-dp / hr + 0.5 * dm) * rm
    return 2 * np.pi * hr * np.sum(dens * rm), 2 * np.pi * hr * g[1:-1]


opr = minimize(Erad, np.pi * np.exp(-rr1[1:-1] / 2.0), jac=True, method="L-BFGS-B",
               options=dict(maxiter=50000, ftol=1e-15, gtol=1e-12))
log(f"  profilo radiale (1D, quasi continuo): E = {opr.fun:.4f}  -> differenza col 2D completo "
    f"{(E2b+E4b+E0b)/opr.fun-1:+.2e} (errore della griglia 2D, h = {hh:.3f})")
n0b = hedgehog(lambda rr: np.pi * np.exp(-rr / 2.0))
nb, histb = arrested_flow(n0b, 0.0, mu2, 2400, 0.02 * hh, every=150)
for hst in histb:
    log(f"    kappa=0  iter {hst['it']:6d}: E = {hst['E']:.5f}  Q = {hst['Q']:.4f}  raggio(n3<0) = {hst['size']:.4f}")
log("  SENZA termine di Skyrme: vedi raggio e Q qui sopra (Derrick 2D: E2 invariante di scala, E0 ~ raggio^2)")
res["baby_skyrmion"] = dict(kappa2=kap2, mu2=mu2, n=n2, L=L2, Q_start=Q_start, E_start=E_start, Q_end=Qf,
                            E=E2b + E4b + E0b, E2=E2b, E4=E4b, E0=E0b, anisotropy=aniso, E_radial=opr.fun,
                            history=hist, no_skyrme_history=histb)
np.savez_compressed(os.path.join(OUT, "trap_baby.npz"), q=q, n=nf, L=L2)

with open(os.path.join(OUT, "trap_derrick.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
log("scritto lab/out/trap_derrick.json")
