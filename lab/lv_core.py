"""lv_core.py - limiti di Lorentz per la regola zig-zag (studio LIMITI_LORENTZ, prefisso lv_).

Formule chiuse della dispersione della regola congelata (lab/zz_core.py), derivate a mano/sympy
(prima delle corse) e controllate numericamente in lv_r1.py:

  LUCE (spin 1, due tick palindromici x,y,z | z,y,x):
      cos(omega dt) = cos(kx a) cos(ky a) cos(kz a)
  ELETTRONE (spin 1/2, Dirac, theta = m c a / hbar):
      cos(2 omega dt) = cos^2(theta) cos(2kx a) cos(2ky a) cos(2kz a) - sin^2(theta)
  => massless:  omega_e(k) = omega_gamma(2k) / 2      (SU(2) -> SO(3): la pallina di spin 1/2 gira
                                                       di un angolo doppio)
  UN SOLO TICK (ordine fisso x,y,z, NON la regola: solo diagnosi del palindromo):
      cos(omega/2) = c_x c_y c_z + s s_x s_y s_z   (c_j = cos(k_j a/2), s_j = sin(k_j a/2), s = +-1)

Sviluppi (hbar = c = 1, E_a = hbar c / a, A(n) = nx^2 ny^2 + ny^2 nz^2 + nz^2 nx^2 in [0, 1/3]):
  luce:      omega = k [1 - (A/6) (k/E_a)^2 + O(k^4)]      v_g = 1 - (A/2)(k/E_a)^2
  elettrone: E^2 = m^2 + p^2 - (4/3) A p^4/E_a^2 - m^2 p^2/(3 E_a^2) + ...
"""
import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zz_core as zz

# ---- costanti (CODATA 2018) ----
HBARC_GEVM = 1.97326980459e-16      # hbar c in GeV m
E_PL = 1.220890e19                  # energia di Planck, GeV
L_P = 1.616255e-35                  # lunghezza di Planck, m
ME_GEV = 0.51099895e-3              # m_e c^2, GeV
LAMBDA_C = zz.LAMBDA_C              # 3.8616e-13 m


def A_of(n):
    n = np.asarray(n, float)
    n = n / np.linalg.norm(n, axis=-1, keepdims=True)
    x, y, z = n[..., 0], n[..., 1], n[..., 2]
    return x * x * y * y + y * y * z * z + z * z * x * x


# ---- formule chiuse (adimensionali: k = k a, theta = m a; omega in unita' 1/dt) ----
def omega_light_closed(ka):
    ka = np.atleast_2d(ka)
    return np.arccos(np.clip(np.prod(np.cos(ka), axis=1), -1, 1))


def omega_dirac_closed(ka, theta):
    ka = np.atleast_2d(ka)
    c = np.cos(theta) ** 2 * np.prod(np.cos(2 * ka), axis=1) - np.sin(theta) ** 2
    return 0.5 * np.arccos(np.clip(c, -1, 1))


# ---- numerica esatta dalla regola congelata ----
def omega_numeric(ka, theta=0.0, spin=1):
    """quasi-energie positive (per tick) dalla matrice di Bloch dei due tick di zz_core"""
    W = zz.bloch_batch(np.atleast_2d(ka), theta, 0.5 if spin == 0.5 else 1)
    ev = np.linalg.eigvals(W)
    om = -np.angle(ev) / 2.0
    return np.sort(om, axis=1)


def single_tick_omega(ka, spin=1):
    """quasi-energie per UN tick con ordine fisso x,y,z (diagnosi: senza palindromo)"""
    KA = np.atleast_2d(ka)
    N = KA.shape[0]
    out = []
    for i in range(N):
        T = np.eye(3, dtype=complex)
        for j in (0, 1, 2):
            T = zz.T_axis(KA[i, j], j, +1, 1) @ T
        out.append(np.sort(-np.angle(np.linalg.eigvals(T))))
    return np.array(out)


# ---- cinematica di soglia con le dispersioni ESATTE (unita' fisiche, GeV) ----
def E_photon(k, n, Ea):
    """energia del fotone di impulso k (GeV) lungo n, reticolo E_a = hbar c / a"""
    ka = (k / Ea) * np.asarray(n)[None, :] / np.linalg.norm(n)
    return omega_light_closed(ka) * Ea


def E_electron(p, n, Ea, m=ME_GEV):
    ka = (np.atleast_1d(p)[:, None] / Ea) * np.asarray(n)[None, :] / np.linalg.norm(n)
    return omega_dirac_closed(ka, m / Ea) * Ea


def photon_decay_excess(k, n, Ea, xs, m=ME_GEV):
    """Delta E = E_gamma(k) - E_e(x k) - E_e((1-x) k), collineare lungo n (max su x => soglia)"""
    k = float(k)
    Eg = E_photon(np.array([k]), n, Ea)[0]
    return Eg - E_electron(xs * k, n, Ea, m) - E_electron((1 - xs) * k, n, Ea, m)


def kth_photon_decay_analytic(A, Ea, m=ME_GEV):
    """soglia al primo ordine: A k^4 /(16 E_a^2) = m^2  ->  k_th = 2 sqrt(m E_a) A^(-1/4)"""
    return 2.0 * np.sqrt(m * Ea) * A ** (-0.25)


if __name__ == "__main__":
    rng = np.random.default_rng(1)
    k = rng.uniform(-1, 1, (5, 3))
    print(omega_numeric(k, 0, 1)[:, -1], omega_light_closed(k))
    print(omega_numeric(k, 0.3, 0.5)[:, -1], omega_dirac_closed(k, 0.3))
