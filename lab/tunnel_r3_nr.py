"""tunnel_r3_nr.py - corsa R3: limite non relativistico (Buttiker 1983), forma di Ramos 2020 (barriera gaussiana),
e controllo (POST-HOC, non pre-registrato) della riflessione del rivelatore assorbente ideale della regola.
Esami E9, E10. Uscita lab/out/tunnel_r3.json"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tunnel_core import times_lattice, times_schrod, schrod_steps, lattice_scatter, _Msite, _modes, packet, k_of_omega, step

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
t0 = time.time()
res = {}

# ------------------------------------------------ E9: limite NR, rettangolo kappa0 L = 3
VNR = [0.04, 0.02, 0.01]
NCELLS = [200, 400, 800]
XS = [0.2, 0.4, 0.6, 0.8]
e9 = {}
worst = 0.0
for x in XS:
    rows = {}
    for V0 in VNR:
        L = 3 / np.sqrt(2 * V0)
        vals = []
        for Nc in NCELLS:
            a = L / Nc
            w = 1.0 + x * V0
            r = times_lattice(w, a, a, np.full(Nc, a), np.full(Nc, V0))
            vals.append([r["TV"].real * V0, -r["TV"].imag * V0, r["tphi"] * V0])
        vals = np.array(vals)
        aa = np.array([L / n for n in NCELLS])
        M = np.stack([np.ones(3), aa, aa ** 2], 1)
        ex = np.linalg.solve(M, vals)[0]
        rows[str(V0)] = dict(per_a=vals.tolist(), extrap=ex.tolist())
    # V0/m -> 0 (lineare, minimi quadrati sui 3 punti)
    vv = np.array(VNR)
    Y = np.array([rows[str(V)]["extrap"] for V in VNR])
    M = np.stack([np.ones(3), vv], 1)
    c, *_ = np.linalg.lstsq(M, Y, rcond=None)
    s = times_schrod(x, 1.0, 3 / np.sqrt(2))
    ref = [s["TV"].real, -s["TV"].imag, s["tphi"]]
    err = [abs(c[0][i] / ref[i] - 1) for i in range(3)]
    worst = max(worst, err[0], err[1])
    e9[str(x)] = dict(rows=rows, NR_limit=c[0].tolist(), slope_in_V0=c[1].tolist(), buttiker=ref, rel_err=err)
res["E9"] = dict(data=e9, max_rel_err_tau_y_tau_z=worst, soglia=0.01, PASS=bool(worst < 0.01))
print("E9", worst, time.time() - t0, flush=True)


# ------------------------------------------------ E10: barriera gaussiana (forma di Ramos)
def gauss_lattice(x, V0, kw, Nw):
    w = kw / np.sqrt(2 * V0)
    a = w / Nw
    xs = np.arange(-int(3 * Nw), int(3 * Nw) + 1) * a
    prof = V0 * np.exp(-2 * xs ** 2 / w ** 2)
    # orologio con campo proporzionale al profilo (come in Ramos: la Raman la fa il fascio della barriera),
    # normalizzato al picco: T_V = i d ln t / dV0.  (BUG della prima esecuzione: campo uniforme su +-3w; corretto)
    h = 1e-6
    w0 = 1.0 + x * V0
    th_ = np.full(len(xs), a)
    tp = lattice_scatter(w0, a, a, th_, prof * (1 + h))[0]
    tm = lattice_scatter(w0, a, a, th_, prof * (1 - h))[0]
    TV = 1j * np.log(tp / tm) / (2 * h * V0)
    T = abs(lattice_scatter(w0, a, a, th_, prof)[0]) ** 2
    return TV.real * V0, -TV.imag * V0, T


def gauss_schrod(x, kw, n=1200, h=1e-6):
    w = kw / np.sqrt(2.0)
    xs = np.linspace(-3 * w, 3 * w, n + 1)
    xc = 0.5 * (xs[1:] + xs[:-1])
    dx = xs[1] - xs[0]
    prof = np.exp(-2 * xc ** 2 / w ** 2)
    tp = schrod_steps(x, (1 + h) * prof, dx)
    tm = schrod_steps(x, (1 - h) * prof, dx)
    TV = 1j * np.log(tp / tm) / (2 * h)
    return TV.real, -TV.imag


e10 = {}
worst = 0.0
mono = {}
for kw in (1.0, 2.0):
    ys = []
    for x in (0.3, 0.5, 0.7, 0.9):
        lat = [gauss_lattice(x, 0.01, kw, Nw) for Nw in (25, 50, 100)]
        aa = np.array([1 / 25, 1 / 50, 1 / 100])
        M = np.stack([np.ones(3), aa, aa ** 2], 1)
        ex = np.linalg.solve(M, np.array([[l[0], l[1]] for l in lat]))[0]
        lat5 = gauss_lattice(x, 0.005, kw, 100)
        sch = gauss_schrod(x, kw)
        err = abs(ex[0] / sch[0] - 1)
        worst = max(worst, err)
        ys.append(ex[0])
        e10[f"{kw}_{x}"] = dict(rule_V0_001=ex.tolist(), rule_V0_0005_Nw100=[lat5[0], lat5[1]], schrod=list(sch),
                                rel_err_tau_y=err, T=lat[-1][2])
    mono[str(kw)] = ys
pred_shape = bool(np.all(np.diff(mono["2.0"]) > 0)) and not bool(np.all(np.diff(mono["1.0"]) > 0))
res["E10"] = dict(data=e10, tau_y_V0_vs_E=mono, max_rel_err=worst, soglia=0.01, forma_come_prevista=pred_shape,
                  PASS=bool(worst < 0.01 and pred_shape))
print("E10", worst, mono, time.time() - t0, flush=True)

# Ramos in ms (ESPLORATIVO): Rb87, w = 1.3 um (preso come raggio 1/e^2, ipotesi), E/kB = 1 nK (ordine dichiarato)
HBAR, KB, U = 1.054571817e-34, 1.380649e-23, 1.66053906660e-27
M_RB = 86.909180527 * U
w_m = 1.3e-6
E_nK = 1.0
ram = []
for V_nK in (1.25, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
    V0 = V_nK * 1e-9 * KB
    kw = np.sqrt(2 * M_RB * V0) * w_m / HBAR
    x = E_nK / V_nK
    ty, tz, T = gauss_lattice(x, 0.01, kw, 100)
    ram.append(dict(V0_nK=V_nK, E_over_V0=x, kappa0_w=kw, tau_y_ms=ty * HBAR / V0 * 1e3, tau_z_ms=tz * HBAR / V0 * 1e3, T=T))
res["ramos_ms"] = dict(ipotesi="w = 1.3 um raggio 1/e^2; E/kB = 1 nK; altezza ignota -> scansione", rows=ram,
                       misura_ms=0.62, errore_ms=0.07)
print("Ramos", [(r["V0_nK"], round(r["tau_y_ms"], 3), round(r["T"], 3)) for r in ram], flush=True)

# ------------------------------------------------ POST-HOC: riflessione del rivelatore assorbente ideale
post = {}
for E0 in (2.0, 5.0):
    a = 0.025
    th = a
    Mf = _Msite(E0, a, th, 0.0)
    lam, v, f = _modes(Mf)
    iin = int(np.argmax(f))
    vin = v[:, iin] / np.sqrt(abs(f[iin]))
    vref = v[:, 1 - iin] / np.sqrt(abs(f[1 - iin]))
    r_stat = abs(vin[1] / vref[1]) ** 2               # condizione L(X_D) = 0
    k = np.sqrt(E0 ** 2 - 1)
    vel = k / E0
    sig = 20.0
    x0 = -7 * sig - 10
    nx = int((14 * sig + 400) / a)
    i0 = int((7 * sig + 10 + 7 * sig + 150) / a)
    psi, x, _ = packet(nx, a, th, x0, sig, k_of_omega(E0, a, th), cut=7, i0=i0)
    iD = i0
    nt = int((7 * sig + 10 + 7 * sig + 20) / vel / a)
    for n in range(nt):
        psi = step(psi, np.cos(th), np.sin(th), 1.0)
        psi[:, iD:] = 0
    Rpk = float(np.sum(np.abs(psi) ** 2))
    post[str(E0)] = dict(v=vel, R_packet=Rpk, R_stationary_lattice=float(r_stat), R_continuum=(1 - vel) / (1 + vel))
res["posthoc_absorber"] = post
print("absorber", post, flush=True)
res["time_s"] = time.time() - t0
json.dump(res, open(os.path.join(OUT, "tunnel_r3.json"), "w"), indent=1, default=float)
print("tempo %.1f s" % (time.time() - t0))
