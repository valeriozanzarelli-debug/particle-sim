"""dilat_core.py - DILATAZIONE QUANTISTICA DEL TEMPO con l'orologio interno della regola zig-zag (1+1 D).

La regola e' quella congelata di zz_core (Walk 'dirac1'): R corre a +x, L a -x di una cella per tick (a c),
poi il flip  R' = cos th R - i sin th L,  L' = cos th L - i sin th R,  th = m c a / hbar.  Unita' hbar = c = m = 1,
tick dt = a, quindi th = m a.  Niente potenziali: particella libera.

L'OROLOGIO.  Non aggiungiamo nessun ingranaggio alla regola: l'orologio e' la massa stessa.
  * "orologio a qubit" (sonda, DEBITO dichiarato): due copie dello stesso stato che evolvono con massa
    m - delta/2 e m + delta/2 (un livello interno che pesa un po' di piu' = l'orologio di Zych 2011,
    Smith-Ahmadi 2020, Paige-Plato-Kim 2020: H = sqrt(p^2 + M^2), M = m + H_clock).  Si legge la
    sovrapposizione O(t) = <psi_-(t)|psi_+(t)>:  tempo proprio tau = -arg O / delta, visibilita' V = |O|.
  * "contatore di rotazioni" (dentro la regola, nessun grado di liberta' in piu'): ogni flip riceve un peso
    z; f_w(n) = d ln <psi(n)|psi_z(n)> / d ln z a z = 1 (valore debole del numero di flip).  Identita' esatta
    (TEMPO_TUNNEL.md):  f_w = sin th cos th * d ln O/d th + sin^2 th * n  ->  Im f_w = -(sin 2th / 2a) tau.
    Quindi contare i flip = orologio di massa (delta -> 0).

Due motori, che devono coincidere:
  * nel tempo (zz_core.Walk('dirac1'), tick per tick, spazio reale) - la regola vera;
  * spettrale esatto: la regola libera e' diagonale in k, W(k) = C(th) diag(e^{-ika}, e^{ika}) in SU(2),
    W^n = cos(n w a) - i sin(n w a) N(k),  cos(w a) = cos th cos(k a).
Il continuo di Dirac (H = sigma_z p + m sigma_x, SOLO ESAME) si calcola con la stessa macchina.
"""
import numpy as np
import sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import Walk  # noqa: E402

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


# ------------------------------------------------------------------ propagatori in k (shape (nk, 2, 2))
def W_lattice(k, a, m, z=1.0):
    """un tick libero della regola in k; z = peso extra di ogni flip (z = 1: la regola)"""
    th = m * a
    c, s = np.cos(th), np.sin(th)
    e = np.exp(-1j * k * a)
    W = np.empty((len(k), 2, 2), complex)
    W[:, 0, 0] = c * e
    W[:, 0, 1] = -1j * z * s / e
    W[:, 1, 0] = -1j * z * s * e
    W[:, 1, 1] = c / e
    return W


def U_lattice(k, a, m, n):
    """W(k)^n esatto (W in SU(2))"""
    W = W_lattice(k, a, m)
    cw = 0.5 * np.real(W[:, 0, 0] + W[:, 1, 1])
    wa = np.arccos(np.clip(cw, -1, 1))
    sw = np.sin(wa)
    N = 1j * (W - cw[:, None, None] * I2) / sw[:, None, None]
    return np.cos(n * wa)[:, None, None] * I2 - 1j * np.sin(n * wa)[:, None, None] * N


def omega_lattice(k, a, m):
    return np.arccos(np.cos(m * a) * np.cos(k * a)) / a


def rate_lattice(k, a, m):
    """ritmo dell'orologio di massa della regola, d omega / d m (ramo positivo)"""
    wa = np.arccos(np.cos(m * a) * np.cos(k * a))
    return np.sin(m * a) * np.cos(k * a) / np.sin(wa)


def vgroup_lattice(k, a, m):
    wa = np.arccos(np.cos(m * a) * np.cos(k * a))
    return np.cos(m * a) * np.sin(k * a) / np.sin(wa)


def U_dirac(p, m, t):
    """exp(-i H t), H = sigma_z p + m sigma_x (continuo, SOLO ESAME)"""
    E = np.sqrt(p ** 2 + m ** 2)
    H = np.zeros((len(p), 2, 2), complex)
    H[:, 0, 0], H[:, 1, 1] = p, -p
    H[:, 0, 1] = H[:, 1, 0] = m
    return np.cos(E * t)[:, None, None] * I2 - 1j * (np.sin(E * t) / E)[:, None, None] * H


# ------------------------------------------------------------------ spinori di energia positiva
def u_plus_lattice(k, a, m, sign=+1):
    """autovettore di W(k) con W u = e^{-i w a} u, w > 0 (sign=+1) o w < 0 (sign=-1); gauge: R reale >= 0"""
    W = W_lattice(k, a, m)
    cw = 0.5 * np.real(W[:, 0, 0] + W[:, 1, 1])
    sw = np.sin(np.arccos(np.clip(cw, -1, 1)))
    N = 1j * (W - cw[:, None, None] * I2) / sw[:, None, None]       # W = cos - i sin N, autovalori di N = +-1
    return _eig_pm(N, sign)


def u_plus_dirac(p, m, sign=+1):
    E = np.sqrt(p ** 2 + m ** 2)
    N = np.zeros((len(p), 2, 2), complex)
    N[:, 0, 0], N[:, 1, 1] = p / E, -p / E
    N[:, 0, 1] = N[:, 1, 0] = m / E
    return _eig_pm(N, sign)


def _eig_pm(N, sign):
    # N hermitiana a traccia nulla, N^2 = 1: proiettore (1 + sign N)/2 applicato alla colonna piu' grande
    P = 0.5 * (I2 + sign * N)
    c0 = P[:, :, 0]
    c1 = P[:, :, 1]
    use1 = np.linalg.norm(c1, axis=1) > np.linalg.norm(c0, axis=1)
    u = np.where(use1[:, None], c1, c0)
    u = u / np.linalg.norm(u, axis=1)[:, None]
    ph = np.where(np.abs(u[:, 0]) > 1e-14, np.exp(-1j * np.angle(u[:, 0])), np.exp(-1j * np.angle(u[:, 1])))
    return u * ph[:, None]


# ------------------------------------------------------------------ stati
def two_packets(k, cos_al, sin_al, phi, p1, p2, Dl, x1=0.0, x2=0.0):
    """ampiezza scalare (in k) della sovrapposizione cos a |p1, x1> + e^{i phi} sin a |p2, x2>,
    pacchetti gaussiani alla Smith-Ahmadi: g_j(p) = exp(-(p - p_j)^2 / (2 Dl^2)) e^{-i p x_j} (non normalizzati
    singolarmente: lo stato totale viene normalizzato dopo).  Ritorna anche g1, g2 normalizzati."""
    g1 = np.exp(-(k - p1) ** 2 / (2 * Dl ** 2)) * np.exp(-1j * k * x1)
    g2 = np.exp(-(k - p2) ** 2 / (2 * Dl ** 2)) * np.exp(-1j * k * x2)
    g1 = g1 / np.sqrt(np.sum(np.abs(g1) ** 2))
    g2 = g2 / np.sqrt(np.sum(np.abs(g2) ** 2))
    g = cos_al * g1 + np.exp(1j * phi) * sin_al * g2
    return g / np.sqrt(np.sum(np.abs(g) ** 2)), g1, g2


def spinor_state(g, u):
    """psi(k) = g(k) u(k), shape (nk, 2)"""
    return g[:, None] * u


# ------------------------------------------------------------------ orologi (spettrale esatto)
def overlap_spec(psi0, k, n_or_t, m, delta, kind="lattice", a=None):
    """O = sum_k psi0^dag U_-^dag U_+ psi0 (U_+- con massa m +- delta/2), per una lista di tempi/tick."""
    out = []
    for nt in np.atleast_1d(n_or_t):
        if kind == "lattice":
            Up = U_lattice(k, a, m + delta / 2, nt)
            Um = U_lattice(k, a, m - delta / 2, nt)
        else:
            Up = U_dirac(k, m + delta / 2, nt)
            Um = U_dirac(k, m - delta / 2, nt)
        pp = np.einsum("kij,kj->ki", Up, psi0)
        pm = np.einsum("kij,kj->ki", Um, psi0)
        out.append(np.sum(np.conj(pm) * pp))
    return np.array(out)


def tau_spec(psi0, k, n_or_t, m, delta=1e-5, kind="lattice", a=None):
    """tempo proprio letto dall'orologio di massa (delta -> 0): tau = -arg O / delta"""
    O = overlap_spec(psi0, k, n_or_t, m, delta, kind, a)
    return -np.angle(O) / delta, np.abs(O)


def mean_rate(g, k, a=None, m=1.0, kind="lattice", sign=+1):
    """ritmo medio <d omega/dm> su |g(k)|^2 (stato di un solo ramo)"""
    if kind == "lattice":
        r = rate_lattice(k, a, m)
    else:
        r = m / np.sqrt(k ** 2 + m ** 2)
    w = np.abs(g) ** 2
    return sign * np.sum(w * r) / np.sum(w)


# ------------------------------------------------------------------ nel tempo: la regola vera (zz_core)
def to_real_space(psi_k, nx):
    """psi_k shape (nk = nx, 2) sulla griglia numpy fftfreq -> psi(x) shape (2, nx) (convenzione ifft)"""
    psi = np.fft.ifft(psi_k.T, axis=1) * np.sqrt(nx)
    return psi


def run_time_domain(psi_k, nx, a, m, delta, nticks, every, hz=1e-6):
    """evolve con la regola congelata (zz_core.Walk 'dirac1') tre copie: massa m-delta/2, m+delta/2 e
    il contatore di flip (copie con z = 1 +- hz, passo scritto qui: identico a Walk con i flip pesati).
    Ritorna lista di (n, O, f_w, norma, <x>)."""
    psi0 = to_real_space(psi_k, nx)
    wp = Walk("dirac1", nx, a, dim=1, theta=(m + delta / 2) * a)
    wm = Walk("dirac1", nx, a, dim=1, theta=(m - delta / 2) * a)
    w0 = Walk("dirac1", nx, a, dim=1, theta=m * a)
    pp, pm, p0 = psi0.copy(), psi0.copy(), psi0.copy()
    zp, zm = psi0.copy(), psi0.copy()
    th = m * a
    c, s = np.cos(th), np.sin(th)

    def ztick(ps, z):
        R = np.roll(ps[0], 1)
        L = np.roll(ps[1], -1)
        return np.array([c * R - 1j * z * s * L, c * L - 1j * z * s * R])

    x = np.arange(nx) * a
    rec = []
    for n in range(1, nticks + 1):
        pp = wp.tick(pp)
        pm = wm.tick(pm)
        p0 = w0.tick(p0)
        zp = ztick(zp, 1 + hz)
        zm = ztick(zm, 1 - hz)
        if n % every == 0:
            O = np.vdot(pm, pp)
            fw = (np.log(np.vdot(p0, zp)) - np.log(np.vdot(p0, zm))) / (2 * hz)
            rho = np.sum(np.abs(p0) ** 2, axis=0)
            rec.append((n, O, fw, rho.sum(), rho))
    return rec


# ------------------------------------------------------------------ formule della letteratura (SOLO ESAME)
def smith_ahmadi_gq(th, phi, p1, p2, Dl, m=1.0):
    """Smith & Ahmadi, Nat. Commun. 11, 5360 (2020), eq. 25: correzione quantistica gamma_Q^{-1}
    (approssimazione non relativistica, pacchetti nello stesso punto). Il loro theta = il nostro alpha."""
    num = np.cos(phi) * np.sin(2 * th) * ((p2 - p1) ** 2 - 2 * (p2 ** 2 - p1 ** 2) * np.cos(2 * th))
    den = 8 * m ** 2 * (np.cos(phi) * np.sin(2 * th) + np.exp((p2 - p1) ** 2 / (4 * Dl ** 2)))
    return num / den


def smith_ahmadi_gc(th, p1, p2, Dl, m=1.0):
    """eq. 24 con clock B a riposo, Delta_B -> 0 (solo il termine di A)"""
    return 1 - (p1 ** 2 * np.cos(th) ** 2 + p2 ** 2 * np.sin(th) ** 2) / (2 * m ** 2) - Dl ** 2 / (4 * m ** 2)


def gq_offset_nr(th, phi, p1, p2, Dl, d, m=1.0):
    """estensione (nostra, stessa algebra di S-A) a pacchetti separati di d = x2 - x1 nello spazio:
    Delta<1/gamma> = -(<p^2>_sup - <p^2>_mix)/(2 m^2), cross term <g1|p^2|g2> con centro complesso."""
    pb = 0.5 * (p1 + p2)
    O = np.exp(-(p2 - p1) ** 2 / (4 * Dl ** 2) - Dl ** 2 * d ** 2 / 4 - 1j * pb * d)
    pc = pb - 1j * Dl ** 2 * d / 2
    X = np.sin(2 * th) * np.exp(1j * phi) * O        # 2 Re[c s e^{i phi} <g1|g2>]  (con <g1|g2> = O)
    Nn = 1 + np.real(X)
    p2sup = (np.cos(th) ** 2 * (p1 ** 2 + Dl ** 2 / 2) + np.sin(th) ** 2 * (p2 ** 2 + Dl ** 2 / 2)
             + np.real(X * (pc ** 2 + Dl ** 2 / 2))) / Nn
    p2mix = np.cos(th) ** 2 * p1 ** 2 + np.sin(th) ** 2 * p2 ** 2 + Dl ** 2 / 2
    return -(p2sup - p2mix) / (2 * m ** 2)


def lattice_rate_a2(p, m=1.0):
    """previsione analitica (nostra) del primo errore reticolare: rate_lattice = (m/E)[1 + C a^2],
    C = -(p^2/3)(1 - m^2/(2 E^2))   [da omega^2 = p^2 + m^2 - m^2 p^2 a^2 / 3 + O(a^4)]"""
    E2 = p ** 2 + m ** 2
    return -(p ** 2 / 3) * (1 - m ** 2 / (2 * E2))
