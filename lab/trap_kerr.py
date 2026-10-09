"""trap_kerr.py - Parte B2b: una non linearita' tipo Kerr (l'indice cresce con l'intensita') lega la luce?

(i)  Fascio parassiale (la luce corre lungo z a velocita' c, si guarda solo la sezione trasversa):
       i psi_z + lap_perp psi + N(|psi|^2) psi = 0      (equazione di Schroedinger non lineare)
     - 1 dimensione trasversa, N = |psi|^2: solitone sech, STABILE (Zakharov-Shabat 1972)
     - 2 dimensioni trasverse, N = |psi|^2: solitone di Townes, INSTABILE (collasso o diffrazione)
     - 2 dimensioni trasverse, N saturabile |psi|^2/(1+|psi|^2): STABILE
(ii) Luce FERMA (impulso totale zero) in 2D, equazioni di Maxwell complete (modo TM) in un mezzo
     di Kerr istantaneo  D = (1 + De(E^2)) E :
       cubico     De = chi E^2
       saturabile De = delta (E^2/S)/(1+E^2/S)
     un anello che gira (m = 2 pi R / lambda), spugna assorbente ai bordi, RK4 pseudo-spettrale.
     Misura: frazione dell'energia iniziale ancora entro r < 2R a vari tempi, picco del campo.

uso: timeout 590 python lab/trap_kerr.py beam      (fasci parassiali)
     timeout 590 python lab/trap_kerr.py ring K    (anello fermo, caso K = 0..6; i casi girano in parallelo)
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from trap_core import grid2d, tm2d_positive_frequency, TWO_PI

OUT = os.path.join(os.path.dirname(__file__), "out")
T0 = time.time()
rng = np.random.default_rng(7)


def log(*a):
    print(f"[{time.time()-T0:6.1f}s]", *a, flush=True)


PART = sys.argv[1] if len(sys.argv) > 1 else "beam"
CASE = int(sys.argv[2]) if len(sys.argv) > 2 else -1
res = {}

# ============================================================ (i) fasci parassiali
def nls(psi0, KX, KY, dz, nz, nonlin, every):
    k2 = KX**2 + (0 if KY is None else KY**2)
    half = np.exp(-1j * k2 * dz / 2)
    psi = psi0.copy()
    out = []
    for i in range(nz + 1):
        if i % every == 0:
            out.append(psi.copy())
        if i == nz:
            break
        F = np.fft.fftn
        psi = np.fft.ifftn(half * F(psi))
        psi = psi * np.exp(1j * nonlin(np.abs(psi) ** 2) * dz)
        psi = np.fft.ifftn(half * F(psi))
        if not np.isfinite(psi).all() or np.abs(psi).max() > 8 * np.abs(psi0).max():
            out.append(psi.copy())
            break
    return out


if PART == "beam":
    log("(i-a) NLS 1+1: solitone sech")
    n, L = 2048, 160.0
    x = (np.arange(n) - n // 2) * L / n
    k = TWO_PI * np.fft.fftfreq(n, L / n)
    r1 = {}
    for A in (1.0, 1.3):
        psi0 = A * np.sqrt(2) / np.cosh(x) + 0j      # A = 1: solitone esatto sqrt(2) sech(x)
        outs = nls(psi0, k, None, 0.005, 20000, lambda I: I, 2000)   # z = 100 = 100 lunghezze di dispersione
        peaks = [float(np.abs(p).max()) for p in outs]
        fwhm = []
        for p in outs:
            a = np.abs(p) ** 2
            fwhm.append(float((a > a.max() / 2).sum() * L / n))
        r1[str(A)] = dict(z=list(np.arange(len(peaks)) * 10.0), peak=peaks, fwhm=fwhm)
        log(f"  ampiezza {A}: picco lungo z = {np.round(peaks, 3).tolist()}")
        log(f"               FWHM(|psi|^2)  = {np.round(fwhm, 3).tolist()}")
    res["nls_1d"] = r1

    log("(i-b) NLS 2+1: solitone di Townes (critico) e saturabile")
    # profilo di Townes: R'' + R'/r - R + R^3 = 0, R(inf)=0, per tiro
    from scipy.integrate import solve_ivp


    def shoot(R0, sat=False):
        def f(r, y):
            R, dR = y
            nl = R**3 / (1 + R**2) if sat else R**3
            return [dR, -dR / max(r, 1e-9) + R - nl - (0 if not sat else 0)]
        sol = solve_ivp(f, [1e-6, 30], [R0, 0], rtol=1e-11, atol=1e-13, dense_output=True,
                        events=[lambda r, y: y[0], lambda r, y: y[1]])
        return sol


    def ground(sat=False, beta=1.0):
        # stato fondamentale: bisezione sull'ampiezza centrale (sotto: R ha min >0 e risale; sopra: attraversa 0)
        lo, hi = (1.0, 5.0) if not sat else (0.1, 40.0)
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            def f(r, y):
                R, dR = y
                nl = R**3 / (1 + R**2) if sat else R**3
                return [dR, -dR / max(r, 1e-9) + beta * R - nl]
            ev = lambda r, y: y[0]
            ev.terminal = True
            ev2 = lambda r, y: y[1]
            ev2.terminal = True
            nl0 = mid**3 / (1 + mid**2) if sat else mid**3
            sol = solve_ivp(f, [1e-4, 40], [mid, 0.5e-4 * (beta * mid - nl0)], rtol=1e-11, atol=1e-13, events=[ev, ev2])
            if sol.t_events[0].size:     # attraversa zero -> troppo alto
                hi = mid
            else:
                lo = mid
        def f(r, y):
            R, dR = y
            nl = R**3 / (1 + R**2) if sat else R**3
            return [dR, -dR / max(r, 1e-9) + beta * R - nl]
        nl0 = lo**3 / (1 + lo**2) if sat else lo**3
        sol = solve_ivp(f, [1e-4, 40], [lo, 0.5e-4 * (beta * lo - nl0)], rtol=1e-11, atol=1e-13, dense_output=True, events=[ev, ev2])
        return lo, sol


    R0T, solT = ground(False)
    rr = np.linspace(1e-4, 12, 20000)
    RT = np.where(rr <= solT.t[-1], np.clip(solT.sol(np.minimum(rr, solT.t[-1]))[0], 0, None), 0.0)
    RT[rr >= solT.t[-1]] = 0  # taglia la coda della bisezione (dove il tiro si stacca)
    Pcr = TWO_PI * np.trapezoid(RT**2 * rr, rr)
    log(f"  Townes: R(0) = {R0T:.5f}, potenza critica P_cr = {Pcr:.4f} (letteratura 11.70)")
    n2, L2 = 128, 24.0
    X, Y, KX, KY = grid2d(n2, L2)
    rgrid = np.hypot(X, Y)
    Rprof = np.interp(rgrid, rr, RT, right=0.0)
    r2 = {"P_cr": Pcr, "R0": R0T}
    for fac in (0.95, 1.0, 1.05):
        outs = nls(fac * Rprof + 0j, KX, KY, 0.004, 2500, lambda I: I, 125)
        peaks = [float(np.abs(p).max()) for p in outs]
        r2[str(fac)] = dict(peak=peaks, dz=0.5, z_end=0.5 * (len(peaks) - 1))
        status = "COLLASSO (picco > 4x iniziale; poi limitato dalla griglia)" if (max(peaks) > 4 * peaks[0] or not np.isfinite(peaks[-1])) else ("DIFFRAZIONE (picco < 1/2 iniziale)" if peaks[-1] < 0.5 * peaks[0] else "resta circa uguale")
        log(f"  Townes x {fac}: picco ogni 0.5 in z: {np.round(peaks, 2).tolist()}  -> {status}")
    res["nls_2d_kerr"] = r2

    beta = 0.5
    R0S, solS = ground(True, beta)
    RS = np.where(rr <= solS.t[-1], np.clip(solS.sol(np.minimum(rr, solS.t[-1]))[0], 0, None), 0.0)
    PS = TWO_PI * np.trapezoid(RS**2 * rr, rr)
    X, Y, KX, KY = grid2d(128, 40.0)
    rgrid = np.hypot(X, Y)
    Rprof = np.interp(rgrid, rr, RS, right=0.0)
    log(f"  raggio dove R < 1e-3 R(0): {rr[np.argmax(RS < 1e-3 * R0S)]:.2f} (box +/- 20)")
    log(f"  saturabile (beta = {beta}): R(0) = {R0S:.4f}, potenza {PS:.3f}")
    r3 = {"R0": R0S, "P": PS}
    pert = 1.10 * Rprof * (1 + 0.05 * np.cos(2 * np.arctan2(Y, X)))   # +10% e deformazione ellittica
    outs = nls(pert + 0j, KX, KY, 0.02, 5000, lambda I: I / (1 + I), 500)
    peaks = [float(np.abs(p).max()) for p in outs]
    mom2 = [float(np.sum(np.abs(p) ** 2 * rgrid**2) / np.sum(np.abs(p) ** 2)) for p in outs]
    r3.update(peak=peaks, r2=mom2, dz=10.0)
    log(f"  saturabile perturbato (+10%, ellittico 5%): picco ogni 10 in z: {np.round(peaks, 3).tolist()}")
    log(f"     <r^2>: {np.round(mom2, 3).tolist()}   (profilo esatto: picco {R0S:.3f})")
    res["nls_2d_sat"] = r3


# ============================================================ (ii) luce ferma, Maxwell completo
log("(ii) anello di luce FERMO in un mezzo di Kerr istantaneo, Maxwell 2D completo")
n, L = 160, 24.0
X, Y, KX, KY = grid2d(n, L)
r = np.hypot(X, Y)
dx = L / n
edge = np.maximum(np.abs(X), np.abs(Y))
sig = np.where(edge > L / 2 - 4, 2.0 * ((edge - (L / 2 - 4)) / 4) ** 2, 0.0)   # spugna


KXr, KYr = KX[:, : n // 2 + 1], KY[:, : n // 2 + 1]    # griglia per rfft2


def grad2(f):
    F = np.fft.rfft2(f)
    return np.fft.irfft2(1j * KXr * F, s=f.shape), np.fft.irfft2(1j * KYr * F, s=f.shape)


def curlz(Bx, By):     # dBy/dx - dBx/dy
    return np.fft.irfft2(1j * KXr * np.fft.rfft2(By) - 1j * KYr * np.fft.rfft2(Bx), s=Bx.shape)


def make_model(kind, p):
    if kind == "lin":
        De = lambda s: 0 * s
        dDe = lambda s: 0 * s
        intE = lambda s: 0 * s
    elif kind == "cubic":
        De = lambda s: p * s
        dDe = lambda s: p + 0 * s
        intE = lambda s: p * s**2 / 4
    else:
        dl, S = p
        De = lambda s: dl * s / (S + s)
        dDe = lambda s: dl * S / (S + s) ** 2
        intE = lambda s: 0.5 * dl * (s - S * np.log1p(s / S))

    def E_of_D(D, E):
        for _ in range(30):
            s = E * E
            g = E * (1 + De(s)) - D
            gp = 1 + De(s) + 2 * s * dDe(s)
            dE = g / gp
            E = E - dE
            if np.max(np.abs(dE)) < 1e-12:
                break
        return E

    def W(E):
        s = E * E
        return s * (1 + De(s)) - s / 2 - intE(s)
    return E_of_D, W


def run_ring(kind, p, A, R=3.0, w=1.0, tmax=24.0, dt=0.03):
    E_of_D, W = make_model(kind, p)
    m = int(round(TWO_PI * R))
    psi = A * np.exp(-((r - R) / w) ** 2) * np.exp(1j * m * np.arctan2(Y, X))
    Ek, Bxk, Byk, K = tm2d_positive_frequency(psi, KX, KY)
    E = np.fft.ifft2(Ek).real
    Bx = np.fft.ifft2(Bxk).real
    By = np.fft.ifft2(Byk).real
    D = E * (1 + {"lin": 0}.get(kind, 0)) if kind == "lin" else None
    s = E * E
    # D dal campo iniziale
    E_of_D_, W_ = E_of_D, W
    if kind == "lin":
        D = E.copy()
    elif kind == "cubic":
        D = E * (1 + p * s)
    else:
        D = E * (1 + p[0] * s / (p[1] + s))
    Ecur = [E.copy()]

    def rhs(D, Bx, By):
        Ee = E_of_D_(D, Ecur[0])
        Ecur[0] = Ee
        gx, gy = grad2(Ee)
        dD = curlz(Bx, By) - sig * D
        dBx = -gy - sig * Bx
        dBy = gx - sig * By
        return dD, dBx, dBy

    def energy(D, Bx, By):
        Ee = E_of_D_(D, Ecur[0])
        return W_(Ee) + 0.5 * (Bx**2 + By**2), Ee

    e0, _ = energy(D, Bx, By)
    Etot0 = e0.sum()
    inside = r < 2 * R
    rec = dict(t=[], frac_in=[], peakE=[], total=[])
    nsteps = int(round(tmax / dt))
    tsave = set(int(round(tt / dt)) for tt in np.arange(0, tmax + 1e-9, 2.0))
    for i in range(nsteps + 1):
        if i in tsave:
            e, Ee = energy(D, Bx, By)
            rec["t"].append(i * dt)
            rec["frac_in"].append(float(e[inside].sum() / Etot0))
            rec["total"].append(float(e.sum() / Etot0))
            rec["peakE"].append(float(np.abs(Ee).max()))
            log(f"     t={i*dt:5.1f}  frazione entro 2R {rec['frac_in'][-1]:.4f}  energia totale {rec['total'][-1]:.4f}  picco |E| {rec['peakE'][-1]:.3f}")
            if not np.isfinite(e).all() or np.abs(Ee).max() > 50 * A + 5:
                rec["blowup_t"] = i * dt
                break
        if i == nsteps:
            break
        k1 = rhs(D, Bx, By)
        k2 = rhs(D + 0.5 * dt * k1[0], Bx + 0.5 * dt * k1[1], By + 0.5 * dt * k1[2])
        k3 = rhs(D + 0.5 * dt * k2[0], Bx + 0.5 * dt * k2[1], By + 0.5 * dt * k2[2])
        k4 = rhs(D + dt * k3[0], Bx + dt * k3[1], By + dt * k3[2])
        D = D + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        Bx = Bx + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        By = By + dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
    return rec


def t_half(rec):
    t, f = rec["t"], rec["frac_in"]
    for i in range(1, len(t)):
        if f[i] < 0.5 <= f[i - 1]:
            return t[i - 1] + (f[i - 1] - 0.5) * (t[i] - t[i - 1]) / (f[i - 1] - f[i])
    return float("nan")


if PART == "ring":
    ring = []
    for kind, p, A in [("lin", None, 1.0), ("cubic", 1.0, 0.2), ("cubic", 1.0, 0.4), ("cubic", 1.0, 0.6),
                       ("sat", (0.5, 0.04), 0.3), ("sat", (0.5, 0.04), 1.0), ("sat", (1.0, 0.04), 1.0)][CASE:CASE + 1]:
        rec = run_ring(kind, p, A)
        th = t_half(rec)
        idx = {tt: i for i, tt in enumerate(rec["t"])}
        f10 = rec["frac_in"][idx[10.0]] if 10.0 in idx else float("nan")
        f40 = rec["frac_in"][-1]
        dmax = (p if kind == "cubic" else (p[0] if kind == "sat" else 0))
        ring.append(dict(kind=kind, p=p, A=A, t_half=th, frac_t10=f10, frac_end=f40, t_end=rec["t"][-1],
                         blowup_t=rec.get("blowup_t"), peakE=rec["peakE"], frac=rec["frac_in"], t=rec["t"]))
        extra = f" ESPLOSO a t={rec['blowup_t']}" if "blowup_t" in rec else ""
        log(f"  {kind:5s} p={p} A={A}: t1/2(r<2R)={th:6.2f}  frazione a t=10: {f10:.3f}  a t={rec['t'][-1]:.0f}: {f40:.4f}"
            f"  picco |E| max={max(rec['peakE']):.2f}{extra}")
    res["ring_fullwave"] = ring


with open(os.path.join(OUT, f"trap_kerr_{PART}{CASE if CASE >= 0 else ''}.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
log(f"scritto lab/out/trap_kerr_{PART}{CASE if CASE >= 0 else ''}.json")
