"""trap_core.py - strumenti comuni per lo studio "luce intrappolata" (docs/LUCE_INTRAPPOLATA.md).

Unita': c = 1, lunghezza d'onda lambda = 1 (salvo dove detto), quindi k0 = 2*pi.

Perche' un propagatore spettrale "e'" la regola delle palline nel vuoto
-----------------------------------------------------------------------
Nel vuoto la regola delle palline-onda (docs/REGOLA_BASE.md) somma ondine esatte
(funzione di Green dell'equazione di Maxwell / delle onde).  La somma di ondine esatte
di una griglia di sorgenti che campiona il campo sopra Nyquist e' la stessa cosa del
propagatore esatto e^{-i|k|t} applicato nello spazio di Fourier: due modi di scrivere
la stessa soluzione lineare.  `kirchhoff_poisson_2d` fa il controllo esplicito: somma
"a palline" (formula di Poisson 2D, ondine di Green ritardate) contro il propagatore
spettrale.  Sul reticolo periodico il propagatore e' esatto finche' l'onda uscente
non rientra dal lato opposto: gli script si fermano prima (o usano una spugna).
"""
import numpy as np

TWO_PI = 2.0 * np.pi


# ----------------------------------------------------------------------------- griglie
def grid2d(n, L):
    x = (np.arange(n) - n // 2) * (L / n)
    X, Y = np.meshgrid(x, x, indexing="ij")
    k = TWO_PI * np.fft.fftfreq(n, d=L / n)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    return X, Y, KX, KY


def grid3d(n, L):
    x = (np.arange(n) - n // 2) * (L / n)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    k = TWO_PI * np.fft.fftfreq(n, d=L / n)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    return (X, Y, Z), (KX, KY, KZ)


# ------------------------------------------------------------- 2D: modo TM (E_z, Bx, By)
def tm2d_positive_frequency(psi, KX, KY):
    """Dato E_z complesso (a frequenza positiva, e^{-i|k|t}) restituisce gli spettri
    di E_z e (Bx, By) dello stesso modo.  -i w B = -curl(E_z z) => B = curl(E z)/(i w)."""
    K = np.sqrt(KX**2 + KY**2)
    K[0, 0] = 1.0
    Ek = np.fft.fft2(psi)
    Ek[0, 0] = 0.0
    # curl(E z) = (dE/dy, -dE/dx)
    Bxk = (1j * KY * Ek) / (1j * K)
    Byk = (-1j * KX * Ek) / (1j * K)
    return Ek, Bxk, Byk, K


def tm2d_evolve(Ek, Bxk, Byk, K, t):
    """Campo complesso a frequenza positiva al tempo t (propagatore esatto)."""
    ph = np.exp(-1j * K * t)
    return np.fft.ifft2(Ek * ph), np.fft.ifft2(Bxk * ph), np.fft.ifft2(Byk * ph)


def tm2d_energy_density(Ez, Bx, By):
    """Densita' di energia del campo reale (media sul ciclo non presa: campo istantaneo)."""
    return 0.5 * (Ez.real**2 + Bx.real**2 + By.real**2)


# ------------------------------------------------------------- 3D: Maxwell vettoriale
def transverse_project(Fk, K):
    KX, KY, KZ = K
    k2 = KX**2 + KY**2 + KZ**2
    k2[0, 0, 0] = 1.0
    dot = (KX * Fk[0] + KY * Fk[1] + KZ * Fk[2]) / k2
    return np.array([Fk[0] - KX * dot, Fk[1] - KY * dot, Fk[2] - KZ * dot])


def cross_khat(K, Fk):
    KX, KY, KZ = K
    k = np.sqrt(KX**2 + KY**2 + KZ**2)
    k[0, 0, 0] = 1.0
    kx, ky, kz = KX / k, KY / k, KZ / k
    return np.array([ky * Fk[2] - kz * Fk[1], kz * Fk[0] - kx * Fk[2], kx * Fk[1] - ky * Fk[0]])


def maxwell3d_evolve(Fk, K, t):
    """Vettore di Riemann-Silberstein F = E + iB (spettro, trasverso).
    dF/dt = -i curl F  ->  F_k(t) = cos(|k|t) F_k + sin(|k|t) (khat x F_k)  (esatto)."""
    KX, KY, KZ = K
    k = np.sqrt(KX**2 + KY**2 + KZ**2)
    return np.cos(k * t) * Fk + np.sin(k * t) * cross_khat(K, Fk)


def helicity3d(Ek, Bk, K, dV, N):
    """Elicita' magnetica + elettrica: (1/2) Int (A.B + C.E), A = i k x B / k^2 ... con
    curl A = B.  Nello spazio di Fourier (Parseval) con la normalizzazione di numpy."""
    KX, KY, KZ = K
    k2 = KX**2 + KY**2 + KZ**2
    k2[0, 0, 0] = 1.0

    def vecpot(Bk_):  # A_k = i k x B_k / k^2  (curl A = B per B trasverso)
        return 1j * np.array([KY * Bk_[2] - KZ * Bk_[1], KZ * Bk_[0] - KX * Bk_[2],
                              KX * Bk_[1] - KY * Bk_[0]]) / k2

    Ak = vecpot(Bk)
    Ck = vecpot(Ek)
    hm = np.sum((np.conj(Ak) * Bk).real) * dV / N
    he = np.sum((np.conj(Ck) * Ek).real) * dV / N
    return hm, he


# ------------------------------------------------- controllo "a palline" (Poisson 2D)
def kirchhoff_poisson_2d(u0, u1, X, Y, probes, t, nalpha=400, ntheta=720):
    """u(x,t) per l'equazione delle onde 2D con u(0)=u0, u_t(0)=u1, sommando le ondine
    ritardate di Green G = 1/(2 pi sqrt(t^2-r^2)) emesse da ogni punto (le 'palline' del
    campo iniziale).  Formula di Poisson: u = d/dt[t M_t(u0)] + t M_t(u1), con
    M_t(f)(x) = Int_0^{pi/2} Int_0^{2pi} f(x + t sin a (cos th, sin th)) sin a da dth / (2 pi).
    Interpolazione bilineare del campo iniziale."""
    from scipy.interpolate import RegularGridInterpolator
    xs = X[:, 0]
    i0 = RegularGridInterpolator((xs, xs), u0, bounds_error=False, fill_value=0.0)
    i1 = RegularGridInterpolator((xs, xs), u1, bounds_error=False, fill_value=0.0)
    # quadratura di Gauss-Legendre in alpha, rettangoli in theta (periodica)
    ga, gw = np.polynomial.legendre.leggauss(nalpha)
    a = 0.25 * np.pi * (ga + 1.0)
    wa = 0.25 * np.pi * gw
    th = np.linspace(0, TWO_PI, ntheta, endpoint=False)
    wth = TWO_PI / ntheta
    A, TH = np.meshgrid(a, th, indexing="ij")
    W = (wa[:, None] * np.sin(A)) * wth / TWO_PI

    def M(interp, px, py, tt):
        pts = np.stack([px + tt * np.sin(A) * np.cos(TH), py + tt * np.sin(A) * np.sin(TH)], -1)
        return np.sum(W * interp(pts.reshape(-1, 2)).reshape(A.shape))

    out = []
    h = 1e-3
    for (px, py) in probes:
        d = ((t + h) * M(i0, px, py, t + h) - (t - h) * M(i0, px, py, t - h)) / (2 * h)
        out.append(d + t * M(i1, px, py, t))
    return np.array(out)


def radial_fraction(e, R, rmax, dA):
    return np.sum(e[R < rmax]) * dA
