"""tunnel_r2_pack.py - corsa R2: pacchetti d'onda della regola attraverso la barriera (tempo reale, tick per tick).
Uso: python3 tunnel_r2_pack.py PK1 0.025     (un pezzo = una configurazione e un passo; checkpoint = un JSON per pezzo)

Copie evolute insieme (stessa regola, stessa partenza):
  0 base (barriera V0)        1,2 orologio di Larmor V0 +- d/2 dentro la barriera
  3,4 orologio di massa (= contatore di rotazioni) th +- d a/2 dentro la barriera
  5 libero (senza barriera)   6 barriera + rivelatore assorbente ideale in X_D
Misure: cono di luce (esatto), P_barriera(t) (dwell), picco trasmesso vs libero, angoli degli orologi,
traiettorie dei surfisti (quantili 1-D = Bohm), flusso e assorbimento al rivelatore (tempo di arrivo).
Esami E1, E5, E6 (parte pacchetto), E7, E12, E13."""
import json, os, sys, time
import numpy as np
import numba as nb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tunnel_core import packet, k_of_omega, dispersion, lattice_scatter, times_lattice

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
PRE = json.load(open(os.path.join(OUT, "tunnel_prereg.json")))
DELTA = 1e-4
NQ = 200


@nb.njit(cache=True)
def tick_all(psi, out, ct, st, ph):
    C, _, n = psi.shape
    for c in range(C):
        for x in range(n):
            R = psi[c, 0, x - 1] if x > 0 else 0j
            L = psi[c, 1, x + 1] if x < n - 1 else 0j
            out[c, 0, x] = ph[c, x] * (ct[c, x] * R - 1j * st[c, x] * L)
            out[c, 1, x] = ph[c, x] * (ct[c, x] * L - 1j * st[c, x] * R)


@nb.njit(cache=True)
def support(p):
    """primo e ultimo indice con ampiezza non nulla (esattamente != 0)"""
    n = p.shape[1]
    lo, hi = -1, -1
    for x in range(n):
        if p[0, x] != 0 or p[1, x] != 0:
            lo = x
            break
    for x in range(n - 1, -1, -1):
        if p[0, x] != 0 or p[1, x] != 0:
            hi = x
            break
    return lo, hi


@nb.njit(cache=True)
def quantiles_right(p, qs, a, x0idx):
    """posizioni X_q dove la probabilita' alla destra vale q (q crescenti, piccoli): scansione da destra"""
    n = p.shape[1]
    X = np.full(qs.size, np.nan)
    Cum = 0.0
    j = 0
    for x in range(n - 1, -1, -1):
        rho = abs(p[0, x]) ** 2 + abs(p[1, x]) ** 2
        Cnew = Cum + rho
        while j < qs.size and qs[j] <= Cnew:
            frac = (qs[j] - Cum) / rho if rho > 0 else 0.0
            X[j] = (x - x0idx) * a + 0.5 * a - frac * a
            j += 1
        Cum = Cnew
        if j >= qs.size:
            break
    return X


def run(cfg, a, sig_override=None):
    P = dict(PRE["pacchetti"][cfg])
    if sig_override is not None:          # ESPLORATIVO (non pre-registrato): dipendenza dalla larghezza
        s_ = sig_override
        P["t_final"] = P["t_final"] - (abs(P["x0"]) + 4 * P["sigma"]) / P["v"] + (7 * s_ + 10 + 4 * s_) / P["v"]
        P["x0"], P["sigma"] = -7 * s_ - 10, s_
    E0, V0, Lb, sig = P["E"], P["V0"], P["L"], P["sigma"]
    x0 = P["x0"]
    th = a
    k0 = k_of_omega(E0, a, th)
    _, v0 = dispersion(k0, a, th)
    N = int(round(Lb / a))
    tf = P["t_final"]
    nt = int(round(tf / a))
    cut = 7.0
    xmin = x0 - cut * sig - tf - 20
    xmax = x0 + cut * sig + tf + 20
    i0 = int(round(-xmin / a))
    nx = int(round((xmax - xmin) / a)) + 1
    psi0, x, kg = packet(nx, a, th, x0, sig, k0, cut=cut, i0=i0)
    inB = (x > -0.5 * a) & (x < Lb - 0.5 * a)
    XD = Lb + 20.0
    iD = int(round(XD / a)) + i0
    C = 7
    psi = np.ascontiguousarray(np.broadcast_to(psi0, (C, 2, nx))).copy()
    out = np.empty_like(psi)
    ct = np.full((C, nx), np.cos(th))
    st = np.full((C, nx), np.sin(th))
    V = np.zeros((C, nx))
    for c in (0, 1, 2, 3, 4, 6):
        V[c, inB] = V0
    V[1, inB] += DELTA / 2
    V[2, inB] -= DELTA / 2
    ct[3, inB], st[3, inB] = np.cos(th + DELTA * a / 2), np.sin(th + DELTA * a / 2)
    ct[4, inB], st[4, inB] = np.cos(th - DELTA * a / 2), np.sin(th - DELTA * a / 2)
    ph = np.exp(-1j * V * a)
    # previsioni stazionarie della regola allo stesso passo (media sullo spettro del pacchetto)
    ks = k0 + np.linspace(-4, 4, 161) / (2 * sig)
    wgt = np.exp(-2 * sig ** 2 * (ks - k0) ** 2)
    Ts, TVs, Tms, Ds, tph = [], [], [], [], []
    for kk in ks:
        w, _ = dispersion(kk, a, th)
        r = times_lattice(w, a, th, np.full(N, th), np.full(N, V0))
        Ts.append(r["T"]); TVs.append(r["TV"]); Tms.append(r["Tm"]); Ds.append(r["dwell"]); tph.append(r["tphi"])
    Ts, TVs, Tms, Ds = map(np.array, (Ts, TVs, Tms, Ds))
    PT_pred = float(np.sum(wgt * Ts) / wgt.sum())
    stat = dict(PT=PT_pred,
                tau_y=float(np.sum(wgt * Ts * TVs.real) / np.sum(wgt * Ts)),
                tau_z=float(-np.sum(wgt * Ts * TVs.imag) / np.sum(wgt * Ts)),
                tau_my=float(np.sum(wgt * Ts * Tms.real) / np.sum(wgt * Ts)),
                tau_mz=float(-np.sum(wgt * Ts * Tms.imag) / np.sum(wgt * Ts)),
                tphi_k0=float(tph[80]), v0=float(v0))
    qs = PT_pred * (np.arange(NQ) + 0.5) / NQ
    t_inB = np.zeros(NQ)
    Xq_hist = []
    lo0, hi0 = support(psi[0])
    cone_ok = True
    cone_front = []
    PB = np.zeros(nt)
    PBT = np.zeros(nt)
    J = np.zeros(nt)
    Jfree = np.zeros(nt)
    F = np.zeros(nt)
    rho0 = np.sum(np.abs(psi0) ** 2, 0)
    ipk0 = int(np.argmax(rho0))
    t0 = time.time()
    iB0, iB1 = i0, i0 + N           # celle della barriera: [iB0, iB1)
    for n in range(nt):
        # flusso attraverso il legame (iD-1 -> iD) durante questo tick: |R(iD-1)|^2 - |L(iD)|^2
        J[n] = abs(psi[0, 0, iD - 1]) ** 2 - abs(psi[0, 1, iD]) ** 2
        Jfree[n] = abs(psi[5, 0, iD - 1]) ** 2 - abs(psi[5, 1, iD]) ** 2
        tick_all(psi, out, ct, st, ph)
        psi, out = out, psi
        # rivelatore assorbente ideale: toglie tutto quello che e' entrato in x >= X_D
        F[n] = np.sum(np.abs(psi[6, :, iD:]) ** 2)
        psi[6, :, iD:] = 0
        # cono di luce, copia base, ogni tick (esatto)
        lo, hi = support(psi[0])
        if lo != lo0 - (n + 1) or hi != hi0 + (n + 1):
            cone_ok = False
        if n % 200 == 0:
            for c in (1, 2, 3, 4, 5):
                l2, h2 = support(psi[c])
                if l2 != lo0 - (n + 1) or h2 != hi0 + (n + 1):
                    cone_ok = False
            cone_front.append([n + 1, hi - hi0])
        rhoB = np.abs(psi[0, :, iB0:iB1]) ** 2
        PB[n] = rhoB.sum()
        # surfisti (Bohm 1-D): quantili trasmessi
        Xq = quantiles_right(psi[0], qs, a, i0)
        inq = (Xq > -0.5 * a) & (Xq < Lb - 0.5 * a)
        t_inB += inq * a
        if n % max(1, nt // 400) == 0:
            Xq_hist.append([n * a] + list(Xq[::20]))
        # parte della barriera "davanti" alla separatrice q = P_T (appartiene ai trasmessi)
        CB1 = np.sum(np.abs(psi[0, :, iB1:]) ** 2)
        PBT[n] = max(0.0, min(CB1 + PB[n], PT_pred) - min(CB1, PT_pred))
    elapsed = time.time() - t0
    trR = x >= Lb + 5.0
    rfR = x < -5.0
    rho = np.sum(np.abs(psi) ** 2, 1)
    PT = float(rho[0, trR].sum())
    PR = float(rho[0, rfR].sum())
    ov = np.vdot(psi[1][:, trR].ravel(), psi[2][:, trR].ravel())
    tau_y = float(np.angle(ov) / DELTA)
    tau_z = float(np.log(rho[2, trR].sum() / rho[1, trR].sum()) / (2 * DELTA))
    ovR = np.vdot(psi[1][:, rfR].ravel(), psi[2][:, rfR].ravel())
    tau_yR = float(np.angle(ovR) / DELTA)
    ovm = np.vdot(psi[3][:, trR].ravel(), psi[4][:, trR].ravel())
    tau_my = float(np.angle(ovm) / DELTA)
    tau_mz = float(np.log(rho[4, trR].sum() / rho[3, trR].sum()) / (2 * DELTA))
    dwell = float(PB.sum() * a)

    def peak(r, mask):
        idx = np.where(mask)[0]
        rr = r[idx]
        j = int(np.argmax(rr))
        w = int(round(0.5 * sig / a))
        sl = slice(max(j - w, 0), min(j + w + 1, len(idx)))
        xx = x[idx][sl]
        yy = np.log(rr[sl])
        c2, c1, _ = np.polyfit(xx - xx.mean(), yy, 2)
        return float(xx.mean() - c1 / (2 * c2))
    xpT = peak(rho[0], trR)
    xpF = peak(rho[5], x > -1e9)
    xp0 = peak(rho0, x > -1e9)
    tfin = nt * a
    adv = xpT - xpF
    margin = xpT - (xp0 + tfin)
    # Bohm
    tauB_T = float(t_inB.mean())
    tauB_T_min = float(t_inB.min())
    tauB_T_dens = float(PBT.sum() * a / PT_pred)
    tauB_R = float((PB - PBT).sum() * a / max(PR, 1e-300))
    # arrivo
    tt = (np.arange(nt) + 1) * a
    def mean_t(f):
        return float(np.sum(tt * f) / np.sum(f))
    passed = float(rho[0, x >= XD - 0.5 * a].sum())
    res = dict(cfg=cfg, a=a, nx=nx, nt=nt, elapsed_s=elapsed, k0=float(k0), v0=float(v0),
               stationary=stat, PT=PT, PR=PR, PT_pred=PT_pred, dwell=dwell,
               tau_y=tau_y, tau_z=tau_z, tau_yR=tau_yR, tau_my=tau_my, tau_mz=tau_mz,
               dwell_identity=PT * tau_y + PR * tau_yR,
               cone_exact=bool(cone_ok), cone_front=cone_front[:5] + cone_front[-3:],
               front_speed=float((cone_front[-1][1]) / cone_front[-1][0]),
               peak=dict(x0=xp0, xT=xpT, xfree=xpF, t=tfin, advance=adv, margin_vs_light=margin,
                         advance_pred_lattice=float(Lb - v0 * stat["tphi_k0"])),
               bohm=dict(tauB_T_traj=tauB_T, tauB_T_min=tauB_T_min, tauB_T_density=tauB_T_dens, tauB_R=tauB_R,
                         identity=PT * tauB_T + PR * tauB_R, Xq_hist=Xq_hist),
               arrival=dict(XD=XD, int_F=float(F.sum()), int_J=float(J.sum()), passed=passed,
                            mean_t_absorber=mean_t(F), mean_t_flux=mean_t(J), mean_t_free=mean_t(Jfree),
                            min_J_over_max=float(J.min() / J.max()),
                            F_series=F[::max(1, nt // 300)].tolist(), J_series=J[::max(1, nt // 300)].tolist(),
                            Jfree_series=Jfree[::max(1, nt // 300)].tolist()),
               PB_series=PB[::max(1, nt // 300)].tolist(),
               snap=dict(x=x[::max(1, int(1 / a) // 4)].tolist(), rhoT=rho[0, ::max(1, int(1 / a) // 4)].tolist(),
                         rhoF=rho[5, ::max(1, int(1 / a) // 4)].tolist(), rho0=rho0[::max(1, int(1 / a) // 4)].tolist()))
    fn = os.path.join(OUT, f"tunnel_r2_{cfg}_{a}.json" if sig_override is None else f"tunnel_r2x_{cfg}_{a}_s{sig}.json")
    json.dump(res, open(fn, "w"), default=float)
    print(cfg, a, f"{elapsed:.1f}s", "PT", PT, PT_pred, "tau_y", tau_y, stat["tau_y"], "tau_z", tau_z, stat["tau_z"],
          "tau_my", tau_my, stat["tau_my"], "tau_mz", tau_mz, stat["tau_mz"],
          "dwell", dwell, PT * tau_y + PR * tau_yR, "cone", cone_ok, "adv", adv, res["peak"]["advance_pred_lattice"],
          "margin", margin, "Bohm", tauB_T, tauB_T_min, tauB_T_dens, "arr", res["arrival"]["int_F"], passed,
          res["arrival"]["mean_t_absorber"], res["arrival"]["mean_t_flux"], res["arrival"]["mean_t_free"], flush=True)


if __name__ == "__main__":
    run(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]) if len(sys.argv) > 3 else None)
