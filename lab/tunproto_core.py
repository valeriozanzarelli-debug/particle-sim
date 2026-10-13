"""tunproto_core.py - PROTOCOLLO per misurare l'orologio di massa nel tunnel (studio tunproto_).

Un'ANCELLA (qubit, stati 0/1) cambia la massa dentro la barriera di +-delta/2 (orologio di massa) oppure il
potenziale di +-delta/2 (orologio di Larmor/Buttiker, il riferimento). Si prepara l'ancella in (|0>+|1>)/sqrt2 e si
legge il suo vettore di Bloch sui soli trasmessi:  s_z = -tanh(delta tau_z),  s_x + i s_y = e^{-i delta tau_y}/cosh.
Stimatori:  tau_y^est = -arg(s_x + i s_y)/delta,   tau_z^est = -atanh(s_z)/delta.

Tre modi di calcolare la stessa misura:
  rule_run()      regola zig-zag (tunnel_core.step) tick per tick, due copie = i due rami dell'ancella
  spectral_rule() previsione dallo stazionario esatto della regola mediato sullo spettro del pacchetto
  dirac_run()     Dirac continuo 1+1 D (H = sigma_z p + m sigma_x + V), split-step di Fourier - SOLO ESAME
  dirac_spectral() lo stesso, dallo stazionario di Dirac (tunnel_core.dirac_scatter) - SOLO ESAME
  wg_*            modello delle guide d'onda binarie (equazioni dei modi accoppiati) = modello dell'APPARATO
Unita' hbar = c = m = 1 (per la regola: a = passo, th = m a).
"""
import os, sys
import numpy as np
import numba as nb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tunnel_core import packet, k_of_omega, dispersion, lattice_scatter, dirac_scatter, bloch_plus  # noqa

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


# ============================================================ lettura dell'ancella
def bloch_from_amps(p0, p1, c01):
    """p0, p1 = probabilita' trasmesse nei due rami, c01 = sum psi0 psi1^* (coerenza). -> (s_x, s_y, s_z, P_T)"""
    P = p0 + p1
    s = 2 * c01 / P
    return np.real(s), np.imag(s), (p0 - p1) / P, P / 2


def estimators(sx, sy, sz, delta):
    ty = -np.angle(sx + 1j * sy) / delta
    tz = -np.arctanh(np.clip(sz, -1 + 1e-15, 1 - 1e-15)) / delta
    return ty, tz


# ============================================================ profili di barriera
def profile(xc, L, w):
    """barriera alta 1 su [0, L]; bordi lisci di larghezza w (tanh); w = 0: rettangolo. xc = centri delle celle"""
    if w <= 0:
        return ((xc >= 0) & (xc < L)).astype(float)
    return 0.5 * (np.tanh(xc / w) - np.tanh((xc - L) / w))


# ============================================================ regola: tick per tick
@nb.njit(cache=True)
def _tick(psi, out, ct, st, ph):
    C, _, n = psi.shape
    for c in range(C):
        for x in range(n):
            R = psi[c, 0, x - 1] if x > 0 else 0j
            L = psi[c, 1, x + 1] if x < n - 1 else 0j
            out[c, 0, x] = ph[c, x] * (ct[c, x] * R - 1j * st[c, x] * L)
            out[c, 1, x] = ph[c, x] * (ct[c, x] * L - 1j * st[c, x] * R)


def rule_setup(a, E0, sigma, x0, xmax, xmin=None, cut=8.0):
    th = a
    k0 = k_of_omega(E0, a, th)
    _, v0 = dispersion(k0, a, th)
    if xmin is None:
        xmin = x0 - 9 * sigma
    i0 = int(round(-xmin / a))
    nx = i0 + int(round(xmax / a))
    psi, x, kg = packet(nx, a, th, x0, sigma, k0, cut=cut, i0=i0)
    return psi, x, kg, k0, v0


def rule_run(a, E0, V0, L, sigma, delta, clock="m", w=0.0, x0=None, tfinal=None, margin=None,
             noise=None, seed=0, nreal=1, pdl=0.0, clock_profile=False, cut=8.0, nt=None):
    """Misura con l'ancella sulla regola. clock: 'm' (massa) o 'V' (potenziale) o 'mV' (entrambi, 4 copie).
    noise: None oppure dict(phase=s_phi, coin=s_th, common=True): rumore per sito e per tick
      phase: fase casuale e^{i phi}, phi ~ N(0, s_phi^2) (potenziale rumoroso del cammino);
      coin:  angolo della moneta th -> th + N(0, s_th^2) (massa rumorosa);
      common: lo stesso rumore nei due rami dell'ancella (True) o indipendente (False).
    pdl: perdita dipendente dalla polarizzazione: la componente R perde una frazione pdl di ampiezza^2 per tick
    clock_profile: l'orologio segue il profilo liscio (True) o e' un rettangolo su [0, L] (False).
    Ritorna dict con s_x, s_y, s_z, P_T, stime tau_y, tau_z (per ogni orologio richiesto)."""
    th = a
    if x0 is None:
        x0 = -6 * sigma - 4 * max(w, 0)
    if margin is None:
        margin = 4 * max(w, 0) + 2.0
    k0 = k_of_omega(E0, a, th)
    _, v0 = dispersion(k0, a, th)
    if tfinal is None:
        tfinal = (L + margin - x0 + 6 * sigma) / v0 + 8 * sigma
    xmax = max(L + margin, x0 + 8 * sigma) + tfinal + 10 * a
    if nt is not None:
        tfinal = nt * a
        xmax = L + margin + tfinal + 10 * a + max(0.0, x0 + cut * sigma)
    psi0, x, kg, k0, v0 = rule_setup(a, E0, sigma, x0, xmax, xmin=x0 - cut * sigma - 2 * a, cut=cut)
    nx = x.size
    xc = x
    prof = profile(xc, L, w)
    clk = prof.copy() if clock_profile else ((xc >= 0) & (xc < L)).astype(float)
    clocks = {"m": ["m"], "V": ["V"], "mV": ["m", "V"]}[clock]
    C = 2 * len(clocks)
    thx = np.full((C, nx), th)
    Vx = np.tile(V0 * prof, (C, 1))
    for j, c in enumerate(clocks):
        for b, sgn in ((0, +1), (1, -1)):
            if c == "m":
                thx[2 * j + b] += sgn * 0.5 * delta * a * clk
            else:
                Vx[2 * j + b] += sgn * 0.5 * delta * clk
    nt = int(round(tfinal / a))
    sel = xc >= L + margin
    rng = np.random.default_rng(seed)
    acc = {c: np.zeros(3, complex) for c in clocks}     # p0, p1, c01 (mediati sulle realizzazioni)
    for r in range(nreal):
        psi = np.repeat(psi0[None], C, axis=0).copy()
        out = np.empty_like(psi)
        ct0, st0 = np.cos(thx), np.sin(thx)
        ph0 = np.exp(-1j * Vx * a)
        for n in range(nt):
            if noise is None:
                ct, st, ph = ct0, st0, ph0
            else:
                com = noise.get("common", True)
                shp = (1, nx) if com else (C, nx)
                ph = ph0
                ct, st = ct0, st0
                if noise.get("phase", 0) > 0:
                    ph = ph0 * np.exp(1j * noise["phase"] * rng.standard_normal(shp))
                if noise.get("coin", 0) > 0:
                    dth = noise["coin"] * rng.standard_normal(shp)
                    ct, st = np.cos(thx + dth), np.sin(thx + dth)
            _tick(psi, out, ct, st, ph)
            psi, out = out, psi
            if pdl > 0:
                psi[:, 0, :] *= np.sqrt(1 - pdl)
        for j, c in enumerate(clocks):
            p0 = np.sum(np.abs(psi[2 * j][:, sel]) ** 2)
            p1 = np.sum(np.abs(psi[2 * j + 1][:, sel]) ** 2)
            c01 = np.sum(psi[2 * j][:, sel] * np.conj(psi[2 * j + 1][:, sel]))
            acc[c] += np.array([p0, p1, c01]) / nreal
    res = dict(v0=float(v0), k0=float(k0), tfinal=float(tfinal), nt=nt, nx=nx,
               light_edge_ok=bool(np.all(psi[:, :, -1] == 0)))
    for c in clocks:
        p0, p1, c01 = acc[c]
        sx, sy, sz, PT = bloch_from_amps(np.real(p0), np.real(p1), c01)
        ty, tz = estimators(sx, sy, sz, delta)
        res[c] = dict(sx=float(sx), sy=float(sy), sz=float(sz), PT=float(PT), ty=float(ty), tz=float(tz),
                      r=float(np.hypot(sx, sy)))
    return res


# ============================================================ regola: previsione spettrale (stazionario esatto)
def _zone(a, L, w):
    ext = 4 * w if w > 0 else 0.0
    i_lo = int(np.floor(-ext / a))
    i_hi = int(np.ceil((L + ext) / a))
    xc = np.arange(i_lo, i_hi) * a
    return xc


def t_rule(om, a, V0, L, w=0.0, dm=0.0, dV=0.0, clock_profile=False):
    th = a
    xc = _zone(a, L, w)
    prof = profile(xc, L, w)
    clk = prof if clock_profile else ((xc >= 0) & (xc < L)).astype(float)
    t, r = lattice_scatter(om, a, th, th + dm * a * clk, V0 * prof + dV * clk)
    # la zona comincia in i_lo: lattice_scatter misura t rispetto alla propagazione libera sulla zona -> ok
    return t


def spectral_weights(a, E0, sigma, nk=401, width=6.0):
    """spettro del pacchetto (stesso di tunnel_core.packet): |g(k)|^2 ~ exp(-2 sigma^2 (k-k0)^2)"""
    th = a
    k0 = k_of_omega(E0, a, th)
    sk = 1 / (2 * sigma)
    k = k0 + np.linspace(-width, width, nk) * sk
    k = k[(k > 0) & (k < np.pi / a)]
    wgt = np.exp(-2 * sigma ** 2 * (k - k0) ** 2)
    om = np.arccos(np.cos(th) * np.cos(k * a)) / a
    return k, om, wgt / wgt.sum()


def spectral_rule(a, E0, V0, L, sigma, delta, clock="m", w=0.0, clock_profile=False, nk=401):
    """previsione della lettura dell'ancella per il pacchetto, dallo stazionario della regola"""
    k, om, wgt = spectral_weights(a, E0, sigma, nk)
    t0 = np.empty(len(k), complex)
    t1 = np.empty(len(k), complex)
    for j, o in enumerate(om):
        if clock == "m":
            t0[j] = t_rule(o, a, V0, L, w, dm=+delta / 2, clock_profile=clock_profile)
            t1[j] = t_rule(o, a, V0, L, w, dm=-delta / 2, clock_profile=clock_profile)
        else:
            t0[j] = t_rule(o, a, V0, L, w, dV=+delta / 2, clock_profile=clock_profile)
            t1[j] = t_rule(o, a, V0, L, w, dV=-delta / 2, clock_profile=clock_profile)
    p0 = np.sum(wgt * abs(t0) ** 2)
    p1 = np.sum(wgt * abs(t1) ** 2)
    c01 = np.sum(wgt * t0 * np.conj(t1))
    sx, sy, sz, PT = bloch_from_amps(p0, p1, c01)
    ty, tz = estimators(sx, sy, sz, delta)
    return dict(sx=float(sx), sy=float(sy), sz=float(sz), PT=float(PT), ty=float(ty), tz=float(tz),
                r=float(np.hypot(sx, sy)))


def stationary_rule(a, E, V0, L, w=0.0, h=1e-6, clock_profile=False):
    """tempi stazionari della regola a quasi-energia E: T, tau_my, tau_mz, tau_y, tau_z (derivate esatte)"""
    t = t_rule(E, a, V0, L, w, clock_profile=clock_profile)
    Tm = 1j * np.log(t_rule(E, a, V0, L, w, dm=h, clock_profile=clock_profile) /
                     t_rule(E, a, V0, L, w, dm=-h, clock_profile=clock_profile)) / (2 * h)
    TV = 1j * np.log(t_rule(E, a, V0, L, w, dV=h, clock_profile=clock_profile) /
                     t_rule(E, a, V0, L, w, dV=-h, clock_profile=clock_profile)) / (2 * h)
    return dict(T=float(abs(t) ** 2), tau_my=float(Tm.real), tau_mz=float(-Tm.imag),
                tau_y=float(TV.real), tau_z=float(-TV.imag))


# ============================================================ Dirac continuo (solo esame)
def dirac_t(E, V0, L, dm=0.0, dV=0.0):
    t, r = dirac_scatter(E, 1.0, 1.0 + dm, V0 + dV, L)
    return t


def stationary_dirac(E, V0, L, h=1e-6):
    t = dirac_t(E, V0, L)
    Tm = 1j * np.log(dirac_t(E, V0, L, dm=h) / dirac_t(E, V0, L, dm=-h)) / (2 * h)
    TV = 1j * np.log(dirac_t(E, V0, L, dV=h) / dirac_t(E, V0, L, dV=-h)) / (2 * h)
    return dict(T=float(abs(t) ** 2), tau_my=float(Tm.real), tau_mz=float(-Tm.imag),
                tau_y=float(TV.real), tau_z=float(-TV.imag))


def dirac_spectral(E0, V0, L, sigma, delta, clock="m", nk=401, width=6.0):
    k0 = np.sqrt(E0 ** 2 - 1)
    sk = 1 / (2 * sigma)
    k = k0 + np.linspace(-width, width, nk) * sk
    k = k[k > 0]
    wgt = np.exp(-2 * sigma ** 2 * (k - k0) ** 2)
    wgt /= wgt.sum()
    E = np.sqrt(k ** 2 + 1)
    kw = dict(dm=delta / 2) if clock == "m" else dict(dV=delta / 2)
    kw1 = {kk: -vv for kk, vv in kw.items()}
    t0 = np.array([dirac_t(e, V0, L, **kw) for e in E])
    t1 = np.array([dirac_t(e, V0, L, **kw1) for e in E])
    p0, p1 = np.sum(wgt * abs(t0) ** 2), np.sum(wgt * abs(t1) ** 2)
    c01 = np.sum(wgt * t0 * np.conj(t1))
    sx, sy, sz, PT = bloch_from_amps(p0, p1, c01)
    ty, tz = estimators(sx, sy, sz, delta)
    return dict(sx=float(sx), sy=float(sy), sz=float(sz), PT=float(PT), ty=float(ty), tz=float(tz),
                r=float(np.hypot(sx, sy)))


def dirac_run(E0, V0, L, sigma, delta, clock="m", dx=0.025, dt=0.0125, x0=None, margin=2.0, tfinal=None):
    """Dirac continuo, H = sigma_z p + m(x) sigma_x + V(x), split-step di Fourier (Strang). Due copie = rami
    dell'ancella. Pacchetto gaussiano di energia positiva (proiettato in k). SOLO ESAME."""
    if x0 is None:
        x0 = -6 * sigma
    k0 = np.sqrt(E0 ** 2 - 1)
    v0 = k0 / E0
    if tfinal is None:
        tfinal = (L + margin - x0 + 6 * sigma) / v0 + 8 * sigma
    xmin = x0 - 9 * sigma - tfinal          # spazio per il riflesso (bordo periodico)
    xmax = max(L + margin, x0 + 9 * sigma) + tfinal + 10
    n = int(2 ** np.ceil(np.log2((xmax - xmin) / dx)))
    x = xmin + dx * np.arange(n)
    kk = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    # stato di energia positiva per H(k) = sigma_z k + sigma_x: (E + k, 1)/norm
    Ek = np.sqrt(kk ** 2 + 1)
    u = np.stack([Ek + kk, np.ones_like(kk)])
    u /= np.linalg.norm(u, axis=0)
    g = np.exp(-(kk - k0) ** 2 * sigma ** 2) * np.exp(-1j * kk * (x0 - xmin))
    spec = g * u
    psi = np.fft.ifft(spec, axis=1)
    psi /= np.sqrt(np.sum(abs(psi) ** 2))
    inB = ((x >= 0) & (x < L)).astype(float)
    C = 2
    psi = np.repeat(psi[None], C, 0)
    m = np.ones((C, n))
    V = np.tile(V0 * inB, (C, 1))
    if clock == "m":
        m[0] += 0.5 * delta * inB
        m[1] -= 0.5 * delta * inB
    else:
        V[0] += 0.5 * delta * inB
        V[1] -= 0.5 * delta * inB
    h = dt / 2
    cm, sm, pv = np.cos(m * h), np.sin(m * h), np.exp(-1j * V * h)
    kin = np.stack([np.exp(-1j * kk * dt), np.exp(1j * kk * dt)])
    nt = int(round(tfinal / dt))

    def loc(p):
        R, Lc = p[:, 0], p[:, 1]
        return np.stack([pv * (cm * R - 1j * sm * Lc), pv * (cm * Lc - 1j * sm * R)], 1)

    for _ in range(nt):
        psi = loc(psi)
        psi = np.fft.ifft(np.fft.fft(psi, axis=2) * kin[None], axis=2)
        psi = loc(psi)
    sel = x >= L + margin
    p0 = np.sum(abs(psi[0][:, sel]) ** 2)
    p1 = np.sum(abs(psi[1][:, sel]) ** 2)
    c01 = np.sum(psi[0][:, sel] * np.conj(psi[1][:, sel]))
    sx, sy, sz, PT = bloch_from_amps(p0, p1, c01)
    ty, tz = estimators(sx, sy, sz, delta)
    return dict(sx=float(sx), sy=float(sy), sz=float(sz), PT=float(PT), ty=float(ty), tz=float(tz),
                r=float(np.hypot(sx, sy)), tfinal=float(tfinal), n=int(n))


# ============================================================ guide d'onda binarie (modello dell'apparato)
def wg_hamiltonian(N, kap, sig, Phi):
    """i dc_n/dz = -kap (c_{n+1} + c_{n-1}) + (-1)^n sig_n c_n + Phi_n c_n   (Longhi 2010).
    sig, Phi: array per guida. Ritorna la matrice tridiagonale."""
    H = np.diag((-1.0) ** np.arange(N) * sig + Phi)
    H += np.diag(-kap * np.ones(N - 1), 1) + np.diag(-kap * np.ones(N - 1), -1)
    return H


def wg_band(q, kap, sig):
    """bande del reticolo binario (passo = 1 guida, cella di 2 guide): beta = +-sqrt(sig^2 + 4 kap^2 cos^2 q)"""
    return np.sqrt(sig ** 2 + 4 * kap ** 2 * np.cos(q) ** 2)


def wg_packet(N, kap, sig, x0, sx, p0):
    """pacchetto della banda alta del reticolo binario: celle j = guide (2j, 2j+1), H(K) = [[sig, -kap(1+e^{-iK})],
    [-kap(1+e^{iK}), -sig]]; vicino al punto di Dirac K = pi: momento per guida p = (K - pi)/2.
    x0, sx in guide (sx = larghezza in |psi|^2), p0 = momento centrale per guida."""
    M = N // 2
    K = 2 * np.pi * np.fft.fftfreq(M)              # per cella
    p = ((K - np.pi + np.pi) % (2 * np.pi) - np.pi) / 2   # (K - pi)/2 in (-pi/2, pi/2]
    g = np.exp(-(p - p0) ** 2 * sx ** 2) * np.exp(-1j * p * x0)  # posizione in guide: x = 2 j, e^{iKj} = (-1)^j e^{ipx}
    u = np.zeros((M, 2), complex)
    for i, k in enumerate(K):
        Hk = np.array([[sig, -kap * (1 + np.exp(-1j * k))], [-kap * (1 + np.exp(1j * k)), -sig]])
        w, v = np.linalg.eigh(Hk)
        uu = v[:, 1]
        u[i] = uu * np.exp(-1j * np.angle(uu[0]))
    j = np.arange(M)
    cells = np.zeros((M, 2), complex)
    for c in range(2):
        cells[:, c] = np.fft.ifft(g * u[:, c]) * M
    # fase di Bloch: e^{iKj} con K = pi + 2p -> il fattore e^{i pi j} e' gia' nel K della FFT
    psi = cells.reshape(-1)
    psi /= np.linalg.norm(psi)
    return psi


def wg_run(kap, sig, V0, Lg, sx, E0, delta, Z, clock="m", x0=None, margin=6, eta_kappa=0.0, disorder=0.0,
           seed=0, N=None, nreal=1):
    """guide binarie: barriera di Lg guide (Phi = V0) a partire dalla guida i_b; orologio = polarizzazione:
    clock 'm': sig_n -> sig +- delta/2 dentro la barriera (birifrangenza alternata); 'V': Phi_n +- delta/2.
    eta_kappa: accoppiamento dipendente dalla polarizzazione kap(1 +- eta/2) OVUNQUE (contaminazione).
    disorder: rumore di fabbricazione sulle costanti di propagazione (std, cm^-1), uguale per le due polarizzazioni.
    E0 = energia (cm^-1) del pacchetto incidente; p0 da beta = sqrt(sig^2 + (2 kap sin p)^2)."""
    c = 2 * kap
    p0 = np.arcsin(np.sqrt(E0 ** 2 - sig ** 2) / (2 * kap))
    if x0 is None:
        x0 = -4 * sx
    if N is None:
        N = int(2 * np.ceil((Lg + margin + c * Z + 8 * sx - x0 + 40) / 2))
    off = int(-x0 + 4 * sx + 20)
    off += off % 2
    n = np.arange(N) - off                       # guida 0 = inizio barriera
    inB = ((n >= 0) & (n < Lg)).astype(float)
    psi0 = wg_packet(N, kap, sig, x0 + off, sx, p0)
    sel = n >= Lg + margin
    rng = np.random.default_rng(seed)
    acc = np.zeros(3, complex)
    for r in range(nreal):
        dis = disorder * rng.standard_normal(N) if disorder > 0 else 0.0
        out = []
        for b, sgn in ((0, 1), (1, -1)):
            sg = sig + (sgn * delta / 2 * inB if clock == "m" else 0)
            Ph = V0 * inB + dis + (sgn * delta / 2 * inB if clock == "V" else 0)
            H = wg_hamiltonian(N, kap * (1 + sgn * eta_kappa / 2), sg, Ph)
            w, U = np.linalg.eigh(H)
            out.append(U @ (np.exp(-1j * w * Z) * (U.conj().T @ psi0)))
        p0_ = np.sum(abs(out[0][sel]) ** 2)
        p1_ = np.sum(abs(out[1][sel]) ** 2)
        c01 = np.sum(out[0][sel] * np.conj(out[1][sel]))
        acc += np.array([p0_, p1_, c01]) / nreal
    sx_, sy_, sz_, PT = bloch_from_amps(np.real(acc[0]), np.real(acc[1]), acc[2])
    ty, tz = estimators(sx_, sy_, sz_, delta)
    # tempo in unita' di 1/sig (massa), lunghezze in Compton = c/sig guide
    return dict(sx=float(sx_), sy=float(sy_), sz=float(sz_), PT=float(PT), ty=float(ty), tz=float(tz),
                ty_m=float(ty * sig), tz_m=float(tz * sig), N=int(N), L_compton=float(Lg * sig / c),
                front_ok=bool(abs(out[0][-1]) < 1e-8))
