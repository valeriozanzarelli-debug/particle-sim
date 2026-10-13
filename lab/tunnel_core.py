"""tunnel_core.py - TEMPO DI TUNNEL e TEMPO DI ARRIVO sulla regola zig-zag (1+1 D).

La regola e' quella congelata di zz_core ('dirac1'): psi = (R, L); R corre di una cella verso +x, L verso -x
(a velocita' c), poi il flip  R' = e^{-iVa}(cos th R - i sin th L),  L' = e^{-iVa}(cos th L - i sin th R),
th = m c a / hbar. Unita' hbar = c = m_e = 1, tick dt = a.

Unica variante (permessa dalle regole comuni): V e th possono dipendere dal sito (barriera di potenziale V0 su N
celle; "orologio" = piccola variazione di V o di th solo dentro la barriera). Niente altro.

Contenuto:
  step()             un tick della regola con th(x), V(x) per sito
  lattice_scatter()  soluzione stazionaria esatta (Floquet) della regola: matrice di trasferimento sito per sito
  dirac_scatter()    lo stesso nel continuo (equazione di Dirac 1+1 D, rappresentazione chirale) - SOLO ESAME
  schrod_t()         barriera rettangolare di Schroedinger (formula chiusa) - SOLO ESAME (limite non relativistico)
  times_*            tempi: fase (Wigner/Hartman), Larmor (Buttiker tau_y, tau_z), dwell, orologio di massa
                     (= conteggio delle rotazioni), tutti come derivate esatte della regola
  packet()           pacchetto d'onda di energia positiva della regola (autovettori di Bloch, fase liscia)
"""
import numpy as np
from scipy.linalg import expm

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2, dtype=complex)


# ===================================================================== la regola (un tick)
def step(psi, ct, st, ph):
    """un tick: R corre a +x, L a -x, poi flip e fase, con ct = cos th(x), st = sin th(x), ph = e^{-i V(x) a}
    (array per sito oppure scalari). psi shape (2, n) oppure (k, 2, n) per k copie (orologi)."""
    R = np.roll(psi[..., 0, :], 1, axis=-1)
    L = np.roll(psi[..., 1, :], -1, axis=-1)
    out = np.empty_like(psi)
    out[..., 0, :] = ph * (ct * R - 1j * st * L)
    out[..., 1, :] = ph * (ct * L - 1j * st * R)
    return out


def dispersion(k, a, theta):
    """quasi-energia omega (ramo positivo) e velocita' di gruppo della regola libera: cos(w a) = cos th cos(k a)"""
    w = np.arccos(np.cos(theta) * np.cos(k * a)) / a
    v = np.cos(theta) * np.sin(k * a) / np.sin(w * a)
    return w, v


def k_of_omega(w, a, theta):
    return np.arccos(np.cos(w * a) / np.cos(theta)) / a


# ===================================================================== stazionario: la regola esatta
def _Msite(w, a, th, V):
    """matrice di trasferimento di un sito: w_x = (R(x-1), L(x)) -> w_{x+1} = (R(x), L(x+1))"""
    mu = np.exp(1j * (w - V) * a)
    c, s = np.cos(th), np.sin(th)
    return np.array([[mu, -1j * s], [1j * s, 1 / mu]]) / c


def _modes(M):
    lam, v = np.linalg.eig(M)
    f = np.abs(v[0]) ** 2 - np.abs(v[1]) ** 2           # flusso per tick del modo (|R|^2 - |L|^2 sul legame)
    return lam, v, f


def lattice_scatter(w, a, theta, thB, VB, want_density=False):
    """Soluzione stazionaria della regola con quasi-energia w (ramo positivo).
    Siti 0..N-1 della 'zona' hanno th = thB[x], V = VB[x]; fuori: th = theta, V = 0.
    Ritorna t (rispetto alla propagazione libera), r, e se richiesto la densita' rho(x) sui siti della zona
    (con flusso incidente 1 per tick)."""
    thB = np.atleast_1d(np.asarray(thB, float))
    VB = np.atleast_1d(np.asarray(VB, float))
    N = len(VB)
    if thB.size == 1:
        thB = np.full(N, thB[0])
    Mf = _Msite(w, a, theta, 0.0)
    lam, v, f = _modes(Mf)
    if not np.allclose(np.abs(lam), 1, atol=1e-9):
        raise ValueError("energia fuori dalla banda libera")
    iin = int(np.argmax(f))
    iref = 1 - iin
    vin = v[:, iin] / np.sqrt(abs(f[iin]))
    vin = vin * np.exp(-1j * np.angle(vin[0]))
    vref = v[:, iref] / np.sqrt(abs(f[iref]))
    T = I2.copy()
    Ms = []
    for x in range(N):
        Mx = _Msite(w, a, thB[x], VB[x])
        Ms.append(Mx)
        T = Mx @ T
    # base dei modi liberi (in, rif): T_mode (1, r) = (tau, 0);  det T = 1  ->  tau = 1 / T_mode[1,1]
    # (forma stabile anche per barriere spesse: niente cancellazioni tra esponenziali grandi)
    Vb = np.stack([vin, vref], 1)
    Tm = np.linalg.solve(Vb, T @ Vb)
    tau = 1.0 / Tm[1, 1]          # det T = 1 esattamente (det M_x = 1): niente cancellazioni
    r = -Tm[1, 0] / Tm[1, 1]
    t = tau / lam[iin] ** N
    if not want_density:
        return t, r
    # densita': propagazione all'indietro dal lato destro (stabile: la soluzione fisica cresce verso sinistra)
    wv = tau * vin
    W = [wv]
    for Mx in Ms[::-1]:
        wv = np.linalg.solve(Mx, wv)
        W.append(wv)
    W = np.array(W[::-1])                 # W[x] = (R(x-1), L(x)), x = 0..N
    rho = np.abs(W[1:, 0]) ** 2 + np.abs(W[:-1, 1]) ** 2     # sito x: |R(x)|^2 + |L(x)|^2
    return t, r, rho


def lattice_scatter_z(w, a, theta, thB, VB, z):
    """Come lattice_scatter, ma ogni ROTAZIONE (flip) fatta dentro la zona riceve un peso extra z sul cammino:
    moneta [[c, -i z s], [-i z s, c]]. Non e' una dinamica (non unitaria per z != 1): e' il CONTATORE di rotazioni
    dei cammini di Feynman a scacchiera, t(z) = sum_f A_f z^f. f_w = d ln t / d ln z a z = 1."""
    thB = np.asarray(thB, float)
    VB = np.asarray(VB, float)
    N = len(VB)
    Mf = _Msite(w, a, theta, 0.0)
    lam, v, f = _modes(Mf)
    iin = int(np.argmax(f))
    iref = 1 - iin
    vin = v[:, iin] / np.sqrt(abs(f[iin]))
    vin = vin * np.exp(-1j * np.angle(vin[0]))
    vref = v[:, iref] / np.sqrt(abs(f[iref]))
    T = I2.copy()
    for x in range(N):
        mu = np.exp(1j * (w - VB[x]) * a)
        c, s = np.cos(thB[x]), np.sin(thB[x])
        c11, c22, c12, c21 = c, c, -1j * z * s, -1j * z * s
        Mx = np.array([[mu * (c11 - c12 * c21 / c22), c12 / c22], [-c21 / c22, 1 / (mu * c22)]])
        T = Mx @ T
    Vb = np.stack([vin, vref], 1)
    Tm = np.linalg.solve(Vb, T @ Vb)
    tau = 1.0 / Tm[1, 1]          # det M_x = c11/c22 = 1 anche con z
    return tau / lam[iin] ** N


def flips_direct(w, a, theta, thB, VB, h=1e-6):
    """numero di rotazioni dentro la zona sui cammini trasmessi (valore debole), contato direttamente con z"""
    tp = lattice_scatter_z(w, a, theta, thB, VB, 1 + h)
    tm = lattice_scatter_z(w, a, theta, thB, VB, 1 - h)
    return np.log(tp / tm) / (2 * h)


# ===================================================================== stazionario: continuo (solo esame)
def dirac_scatter(E, m, mB, VB, L, nseg=1):
    """Dirac 1+1 D, H = sigma_z p + m(x) sigma_x + V(x); barriera rettangolare di lunghezza L (o profilo a
    gradini: mB, VB array di nseg pezzi uguali). Ritorna t (rispetto alla propagazione libera), r."""
    mB = np.atleast_1d(mB) * np.ones(nseg)
    VB = np.atleast_1d(VB) * np.ones(nseg)
    K0 = 1j * E * SZ + m * SY
    lam, v = np.linalg.eig(K0)
    f = np.real(np.einsum("ia,ij,ja->a", v.conj(), SZ, v))
    iin = int(np.argmax(f))
    iref = 1 - iin
    k = np.imag(lam[iin])
    vin = v[:, iin] / np.sqrt(abs(f[iin]))
    vin = vin * np.exp(-1j * np.angle(vin[0]))
    vref = v[:, iref] / np.sqrt(abs(f[iref]))
    T = I2.copy()
    for j in range(nseg):
        T = expm((1j * (E - VB[j]) * SZ + mB[j] * SY) * (L / nseg)) @ T
    Vb = np.stack([vin, vref], 1)
    Tm = np.linalg.solve(Vb, T @ Vb)
    tau = 1.0 / Tm[1, 1]          # det = 1 (K a traccia nulla)
    r = -Tm[1, 0] / Tm[1, 1]
    return tau * np.exp(-1j * k * L), r


def schrod_t(E, V0, L, M=1.0):
    """Schroedinger, barriera rettangolare (hbar = 1, massa M), t rispetto alla propagazione libera."""
    k = np.sqrt(2 * M * E + 0j)
    q = np.sqrt(2 * M * (E - V0) + 0j)
    # t_abs e^{ikL} convention: psi_T = t_abs e^{ikx}; t_rel = t_abs (onda libera = e^{ikx})
    num = 1.0
    den = np.cos(q * L) - 1j * (k ** 2 + q ** 2) / (2 * k * q) * np.sin(q * L)
    return num / den * np.exp(-1j * k * L)


def schrod_steps(E, Vlist, dx, M=1.0):
    """Schroedinger con potenziale a gradini (celle di larghezza dx), t rispetto alla propagazione libera.
    Matrice di trasferimento su (psi, psi'/ik) esatta per ogni gradino. SOLO ESAME."""
    k = np.sqrt(2 * M * E)
    T = np.eye(2, dtype=complex)
    for V in Vlist:
        q = np.sqrt(2 * M * (E - V) + 0j)
        if abs(q) < 1e-14:
            q = 1e-14
        T = np.array([[np.cos(q * dx), np.sin(q * dx) / q], [-q * np.sin(q * dx), np.cos(q * dx)]]) @ T
    # stato (psi, psi'): incidente e^{ikx} + r e^{-ikx} a x=0, trasmesso t e^{ik(x - Ltot)} a x = Ltot
    Ltot = dx * len(Vlist)
    A = np.array([[1, 1], [1j * k, -1j * k]])        # colonne: onda +k, onda -k (in x = 0)
    Bv = np.array([1, 1j * k])
    # T A (1, r) = tt * Bv
    lhs = np.stack([T @ A[:, 1], -Bv], 1)
    r, tt = np.linalg.solve(lhs, -T @ A[:, 0])
    return tt * np.exp(-1j * k * Ltot)


# ===================================================================== tempi (derivate esatte)
def _cdiff(f, x, h):
    return (f(x + h) - f(x - h)) / (2 * h)


def times_lattice(w, a, theta, thB, VB, clock=None, h=1e-6):
    """Tutti i tempi stazionari della regola alla quasi-energia w.
    clock: maschera dei siti dove agisce l'orologio (default: dove VB != 0 o thB != theta, cioe' la barriera).
    Ritorna dict con: T (trasmissione), tphi (fase/Wigner, attraversamento della zona), TV = tau_y - i tau_z
    (Larmor, orologio sul potenziale), Tm = tau_my - i tau_mz (orologio sulla massa), fw (numero di rotazioni
    dentro la zona, valore debole), nw (numero di tick dentro la zona, valore debole), dwell."""
    thB = np.atleast_1d(np.asarray(thB, float)) * np.ones(len(np.atleast_1d(VB)))
    VB = np.atleast_1d(np.asarray(VB, float)) * np.ones(len(thB))
    N = len(VB)
    if clock is None:
        clock = np.ones(N, bool)
    t0, r0, rho = lattice_scatter(w, a, theta, thB, VB, want_density=True)
    lnt = lambda tt: np.log(tt)
    # fase: d arg t / d w  + L / v
    k = k_of_omega(w, a, theta)
    _, vg = dispersion(k, a, theta)
    dphi = np.angle(lattice_scatter(w + h, a, theta, thB, VB)[0] / lattice_scatter(w - h, a, theta, thB, VB)[0]) / (2 * h)
    L = N * a
    tphi = dphi + L / vg
    # orologio sul potenziale (Larmor): T_V = i d ln t / dV
    tp = lattice_scatter(w, a, theta, thB, VB + h * clock)
    tm = lattice_scatter(w, a, theta, thB, VB - h * clock)
    TV = 1j * (np.log(tp[0] / tm[0])) / (2 * h)
    TVr = 1j * (np.log(tp[1] / tm[1])) / (2 * h)
    # orologio sulla massa: T_m = i d ln t / dm  (th = m a)
    sp = lattice_scatter(w, a, theta, thB + h * a * clock, VB)
    sm = lattice_scatter(w, a, theta, thB - h * a * clock, VB)
    Tm = 1j * (np.log(sp[0] / sm[0])) / (2 * h)
    # valori deboli sui cammini trasmessi: n_w = T_V / a ;  f_w = s c (d ln t/d th) + s^2 n_w
    thc = thB[clock][0] if clock.any() else theta
    nw = TV / a
    dlnt_dth = -1j * Tm / a
    fw = np.sin(thc) * np.cos(thc) * dlnt_dth + np.sin(thc) ** 2 * nw
    dwell = a * rho[clock].sum()            # flusso incidente 1 per tick
    return dict(T=abs(t0) ** 2, R=abs(r0) ** 2, t=t0, r=r0, tphi=tphi, TV=TV, TVr=TVr, Tm=Tm, nw=nw, fw=fw,
                dwell=dwell, L=L, vg=vg, k=k)


def times_dirac(E, m, V0, L, h=1e-6):
    """stessi tempi nel continuo di Dirac (barriera rettangolare, stessa massa dentro e fuori). SOLO ESAME."""
    t0, r0 = dirac_scatter(E, m, m, V0, L)
    k = np.sqrt(E ** 2 - m ** 2)
    vg = k / E
    dphi = np.angle(dirac_scatter(E + h, m, m, V0, L)[0] / dirac_scatter(E - h, m, m, V0, L)[0]) / (2 * h)
    tphi = dphi + L / vg
    TV = 1j * np.log(dirac_scatter(E, m, m, V0 + h, L)[0] / dirac_scatter(E, m, m, V0 - h, L)[0]) / (2 * h)
    TVr = 1j * np.log(dirac_scatter(E, m, m, V0 + h, L)[1] / dirac_scatter(E, m, m, V0 - h, L)[1]) / (2 * h)
    Tm = 1j * np.log(dirac_scatter(E, m, m + h, V0, L)[0] / dirac_scatter(E, m, m - h, V0, L)[0]) / (2 * h)
    return dict(T=abs(t0) ** 2, R=abs(r0) ** 2, tphi=tphi, TV=TV, TVr=TVr, Tm=Tm, L=L, vg=vg, k=k)


def times_schrod(E, V0, L, M=1.0, h=1e-7):
    """Buttiker 1983 (Schroedinger, rettangolo): tau_y = -d arg t/dV, tau_z = -d ln|t|/dV; fase. SOLO ESAME."""
    t0 = schrod_t(E, V0, L, M)
    k = np.sqrt(2 * M * E)
    vg = k / M
    TV = 1j * np.log(schrod_t(E, V0 + h * V0, L, M) / schrod_t(E, V0 - h * V0, L, M)) / (2 * h * V0)
    dphi = np.angle(schrod_t(E + h * E, V0, L, M) / schrod_t(E - h * E, V0, L, M)) / (2 * h * E)
    return dict(T=abs(t0) ** 2, TV=TV, tphi=dphi + L / vg)


# ===================================================================== pacchetti
def bloch_plus(k, a, theta):
    """autovettore di energia positiva del tick libero W(k) = C diag(e^{-ika}, e^{ika}), gauge: R reale > 0"""
    c, s = np.cos(theta), np.sin(theta)
    out = np.zeros((len(k), 2), complex)
    om = np.zeros(len(k))
    for j, kk in enumerate(k):
        W = np.array([[c, -1j * s], [-1j * s, c]]) @ np.diag([np.exp(-1j * kk * a), np.exp(1j * kk * a)])
        lam, v = np.linalg.eig(W)
        o = -np.angle(lam) / a
        i = int(np.argmax(o))
        u = v[:, i]
        u = u * np.exp(-1j * np.angle(u[0]))
        out[j] = u / np.linalg.norm(u)
        om[j] = o[i]
    return out, om


def packet(nx, a, theta, x0, sigma, k0, cut=10.0, i0=0):
    """pacchetto gaussiano di energia positiva: centro x0, larghezza sigma (in |psi|^2: sigma_x), impulso k0.
    Griglia x_j = (j - i0) a. Supporto compatto: azzerato fuori da |x - x0| > cut sigma (supporto esatto finito,
    per il test del cono di luce). Ritorna psi (2, nx), x, kgrid."""
    x = (np.arange(nx) - i0) * a
    kg = 2 * np.pi * np.fft.fftfreq(nx, d=a)
    g = np.exp(-(kg - k0) ** 2 * sigma ** 2) * np.exp(-1j * kg * x0)   # |g|^2 ~ exp(-2 s^2 dk^2): sigma_x = sigma
    sel = np.abs(g) > 1e-14
    u, _ = bloch_plus(kg[sel], a, theta)
    psi = np.zeros((2, nx), complex)
    shift = np.exp(-1j * kg * i0 * a)          # x_j = (j - i0) a
    for comp in (0, 1):
        spec = np.zeros(nx, complex)
        spec[sel] = g[sel] * u[:, comp] * shift[sel]
        psi[comp] = np.fft.ifft(spec) * nx
    if cut is not None:
        psi[:, np.abs(x - x0) > cut * sigma] = 0
    psi /= np.sqrt(np.sum(np.abs(psi) ** 2))
    return psi, x, kg


if __name__ == "__main__":
    import time
    t0 = time.time()
    a, th = 0.05, 0.05
    # 1. senza barriera t = 1
    t, r = lattice_scatter(2.0, a, th, np.full(40, th), np.zeros(40))
    print("libero: |t-1| =", abs(t - 1), " |r| =", abs(r))
    # 2. conservazione del flusso con barriera
    N = 200
    t, r = lattice_scatter(2.0, a, th, np.full(N, th), np.full(N, 2.0))
    print("barriera: T + R - 1 =", abs(t) ** 2 + abs(r) ** 2 - 1, " T =", abs(t) ** 2)
    td, rd = dirac_scatter(2.0, 1.0, 1.0, 2.0, N * a)
    print("continuo Dirac: T =", abs(td) ** 2, " T+R-1 =", abs(td) ** 2 + abs(rd) ** 2 - 1, " arg t latt/cont",
          np.angle(t), np.angle(td))
    for aa in (0.1, 0.05, 0.025, 0.0125):
        NN = int(round(10 / aa))
        tl, _ = lattice_scatter(2.0, aa, aa, np.full(NN, aa), np.full(NN, 1.5))
        print(aa, abs(tl) ** 2, np.angle(tl))
    td, _ = dirac_scatter(2.0, 1.0, 1.0, 1.5, 10.0)
    print("cont", abs(td) ** 2, np.angle(td))
    print("time %.1f s" % (time.time() - t0))
