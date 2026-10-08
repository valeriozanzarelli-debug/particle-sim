"""Image formation: shared pieces (camera, wave-level PSF, light sources, atom light rules).

No material tables.  What is in here:
  * physical constants (CODATA, exact or measured constants of nature only);
  * Camera: a DEVICE with its own parameters (pixel pitch, quantum efficiency, dark rate, read noise,
    full well, ADC gain).  These are the specification of a detector you buy, not properties of the scene;
  * airy_psf: level-3 (wave) image of a point emitter through a circular pupil (scalar, paraxial);
    its error against the exact vector Rayleigh–Sommerfeld integral is measured in img_pilot.py;
  * planck_photon_flux: black body (Planck's law, constants only);
  * hydrogen(): EXACT non-relativistic quantum mechanics of hydrogen solved numerically on a radial grid
    (energies, 1s–2p dipole, Einstein A, dynamic polarizability by Dalgarno–Lewis).
    >>> placeholder: replace with our balls' value (docs/MOLECOLE.md, when it exists) <<<
  * dipole collection fraction of a lens of numerical aperture NA.
"""
import numpy as np
from scipy.special import j1
from scipy.linalg import eigh_tridiagonal, solve_banded

# ---- constants of nature (SI, CODATA 2018) ----
H = 6.62607015e-34
HBAR = H / (2 * np.pi)
C = 299792458.0
KB = 1.380649e-23
QE_CHARGE = 1.602176634e-19
ME = 9.1093837015e-31
MP = 1.67262192369e-27
ALPHA = 7.2973525693e-3
A0 = 5.29177210903e-11          # Bohr radius (infinite nuclear mass)
EH = 4.3597447222071e-18        # Hartree energy
T_AU = HBAR / EH                # atomic unit of time
C_AU = 1 / ALPHA                # speed of light in atomic units


# ───────────────────────────── the camera (a device) ─────────────────────────────
class Camera:
    """Turns an expected photon image (photons per pixel arriving at the sensor) into a recorded image.

    Chain, per pixel, all independent:
      detected  = Binomial(arriving, QE)            (each photon clicks with probability QE)
                  -> for a Poisson number of arriving photons this is Poisson(QE * mean)
      dark      = Poisson(dark_rate * exposure)      (thermal electrons, look like photon clicks)
      electrons = detected + dark, clipped at full_well
      signal    = electrons + Normal(0, read_noise)  (amplifier noise, once per read)
      ADU       = round(signal / gain)
    Parameters are those of a camera datasheet; defaults are an idealised scientific camera.
    """

    def __init__(self, qe=1.0, dark_rate=0.0, read_noise=0.0, full_well=np.inf, gain=1.0, quantize=False):
        self.qe, self.dark_rate, self.read_noise = qe, dark_rate, read_noise
        self.full_well, self.gain, self.quantize = full_well, gain, quantize

    def expose(self, mean_photons, exposure, rng):
        """mean_photons: expected photons per pixel arriving in `exposure` seconds (array)."""
        lam = self.qe * np.asarray(mean_photons, float) + self.dark_rate * exposure
        e = rng.poisson(lam).astype(float)
        e = np.minimum(e, self.full_well)
        if self.read_noise > 0:
            e = e + rng.normal(0, self.read_noise, e.shape)
        s = e / self.gain
        return np.round(s) if self.quantize else s

    def clicks(self, arriving_counts, exposure, rng):
        """Same chain starting from a given integer number of arriving photons per pixel (e.g. photon balls
        that landed): QE as a binomial coin per photon."""
        d = rng.binomial(np.asarray(arriving_counts, np.int64), self.qe).astype(float)
        d = d + rng.poisson(self.dark_rate * exposure, d.shape)
        d = np.minimum(d, self.full_well)
        if self.read_noise > 0:
            d = d + rng.normal(0, self.read_noise, d.shape)
        return d / self.gain

    def snr_theory(self, mean_photons, exposure):
        s = self.qe * mean_photons
        return s / np.sqrt(s + self.dark_rate * exposure + self.read_noise ** 2)


# ───────────────────────────── level 3: PSF of a circular pupil ─────────────────────
def airy(r, na, lam):
    """normalised Airy intensity (peak 1) at radius r (same units as lam)."""
    v = 2 * np.pi * na * np.asarray(r, float) / lam
    out = np.ones_like(v)
    nz = v > 1e-9
    out[nz] = (2 * j1(v[nz]) / v[nz]) ** 2
    return out


def pixel_psf(n, pitch, na, lam, cx=0.0, cy=0.0, sub=8):
    """fraction of the photons of a point emitter (image at cx, cy; pixel centres at (i - n//2)·pitch)
    landing in each pixel of an n×n sensor.  Airy integrated over the pixel (sub×sub sub-samples),
    normalised by the analytic total power π/(...)·: the total over an infinite plane is 1."""
    c = (np.arange(n) - n // 2) * pitch
    s = ((np.arange(sub) + 0.5) / sub - 0.5) * pitch
    X = c[None, :, None, None] + s[None, None, None, :]
    Y = c[:, None, None, None] + s[None, None, :, None]
    I = airy(np.hypot(X - cx, Y - cy), na, lam).mean(axis=(2, 3))
    # ∫ airy dA over the plane = lam² / (π na²)
    return I * pitch ** 2 / (lam ** 2 / (np.pi * na ** 2))


def sample_airy_offsets(m, na, lam, rng, rmax_factor=60.0):
    """level-3 photons: m landing offsets drawn from the Airy pattern (radial CDF on a fine grid).
    Encircled energy of Airy: EE(v) = 1 - J0(v)² - J1(v)²  (v = 2π NA r / λ).  rmax: 60 Airy radii."""
    from scipy.special import j0
    vmax = rmax_factor * 3.8317
    v = np.linspace(0, vmax, 200001)
    ee = 1 - j0(v) ** 2 - j1(v) ** 2
    ee /= ee[-1]
    vv = np.interp(rng.random(m), ee, v)
    ph = rng.random(m) * 2 * np.pi
    r = vv * lam / (2 * np.pi * na)
    return r * np.cos(ph), r * np.sin(ph)


# ───────────────────────────── light sources ─────────────────────────────
def planck_photon_spectral_radiance(lam, T):
    """photons / (s · m² · sr · m) at wavelength lam (m), black body at temperature T (Planck)."""
    x = H * C / (lam * KB * T)
    return 2 * C / lam ** 4 / np.expm1(x)


def occupation(lam, T):
    """mean photons per mode n̄ = 1/(e^{hν/kT} − 1)."""
    return 1 / np.expm1(H * C / (lam * KB * T))


# ───────────────────────────── dipole emission collected by a lens ─────────────────
def collection_fraction(na, orient="perp", n=1.0):
    """fraction of the photons of an electric dipole emitter that enter a cone of half-angle asin(NA/n)
    around the optical axis z.  orient: 'perp' (dipole ⊥ axis), 'par' (dipole ∥ axis), 'iso' (average).
    Pattern: dP/dΩ = 3/(8π) sin²ψ, ψ = angle from the dipole axis (exact for an electric dipole)."""
    th = np.arcsin(na / n)
    t = np.linspace(0, th, 20001)
    phi = np.linspace(0, 2 * np.pi, 721)
    T, P = np.meshgrid(t, phi)
    if orient == "par":
        cospsi = np.cos(T)
    elif orient == "perp":
        cospsi = np.sin(T) * np.cos(P)
    else:
        return (1 - np.cos(th)) / 2
    f = 3 / (8 * np.pi) * (1 - cospsi ** 2) * np.sin(T)
    return np.trapezoid(np.trapezoid(f, phi, axis=0), t)


# ───────────────────────────── hydrogen, exact QM on a grid (PLACEHOLDER) ─────────────
def _radial_ops(l, h, rmax, mu):
    """finite-difference radial Hamiltonian for u = r R (atomic units, reduced mass mu in m_e)."""
    r = h * np.arange(1, int(rmax / h) + 1)
    d = 1 / (mu * h * h) + l * (l + 1) / (2 * mu * r * r) - 1 / r
    e = np.full(r.size - 1, -1 / (2 * mu * h * h))
    return r, d, e


def hydrogen(h=0.004, rmax=80.0, reduced_mass=True):
    """>>> placeholder: replace with our balls' value (docs/MOLECOLE.md) <<<
    Exact non-relativistic hydrogen (Schrödinger, Coulomb), solved here on a radial grid with
    Richardson extrapolation h -> h/2 (finite-difference error O(h²)).  Returns a dict in SI + au."""
    mu = 1 / (1 + ME / MP) if reduced_mass else 1.0

    def solve(hh):
        r, d0, e0 = _radial_ops(0, hh, rmax, mu)
        E0, U0 = eigh_tridiagonal(d0, e0, select="i", select_range=(0, 1))
        r1, d1, e1 = _radial_ops(1, hh, rmax, mu)
        E1, U1 = eigh_tridiagonal(d1, e1, select="i", select_range=(0, 0))
        u1s = U0[:, 0] / np.sqrt(hh)
        u2p = U1[:, 0] / np.sqrt(hh)
        u1s *= np.sign(u1s[np.argmax(np.abs(u1s))])
        u2p *= np.sign(u2p[np.argmax(np.abs(u2p))])
        I = np.sum(u1s * r * u2p) * hh               # ∫ R10 r R21 r² dr
        return dict(E1s=E0[0], E2s=E0[1], E2p=E1[0], I=I, r=r, u1s=u1s, d=(d1, e1))

    a, b = solve(h), solve(h / 2)
    ext = {k: (4 * b[k] - a[k]) / 3 for k in ("E1s", "E2s", "E2p", "I")}
    w = ext["E2p"] - ext["E1s"]                       # Hartree
    z2 = ext["I"] ** 2 / 3                            # |<1s|z|2p0>|²  (au, includes 1/μ length scaling)
    A = 4 * w ** 3 * z2 / (3 * C_AU ** 3) / T_AU      # s⁻¹
    lam = 2 * np.pi * C_AU / w * A0                   # m
    out = dict(mu=mu, E1s=ext["E1s"], E2p=ext["E2p"], E2s=ext["E2s"], omega_au=w, dip_z=np.sqrt(z2),
               A21=A, lam_lya=lam, fd_err_E1s=abs(a["E1s"] - b["E1s"]) / 3,
               exact_E1s=-mu / 2, exact_E2=-mu / 8, exact_z=128 * np.sqrt(2) / 243 / mu)
    out["_grid"] = b
    out["_h"] = h / 2
    return out


def polarizability(hyd, omega_au):
    """dynamic dipole polarizability α(ω) of hydrogen 1s (atomic units, a0³), Dalgarno–Lewis:
    α = (1/3)[⟨φ|(H₁−E₀−ω)⁻¹|φ⟩ + ⟨φ|(H₁−E₀+ω)⁻¹|φ⟩], φ = r·u_1s, H₁ = l=1 radial operator.
    Valid below the first resonance (ω < E2p − E1s).  >>> placeholder: replace with our balls' value <<<"""
    g = hyd["_grid"]
    hh = hyd["_h"]
    d1, e1 = g["d"]
    phi = g["r"] * g["u1s"]
    out = []
    for w in np.atleast_1d(omega_au):
        acc = 0
        for s in (+1, -1):
            ab = np.zeros((3, d1.size))
            ab[0, 1:] = e1
            ab[1] = d1 - g["E1s"] + s * w
            ab[2, :-1] = e1
            x = solve_banded((1, 1), ab, phi)
            acc += np.sum(phi * x) * hh
        out.append(acc / 3)
    return np.array(out)


def rayleigh_sigma(alpha_au, lam):
    """Rayleigh cross-section (m²) of a dipole of polarizability alpha (au) at wavelength lam (m):
    σ = (8π/3) k⁴ α_v², α_v = α·a0³ (volume polarizability)."""
    k = 2 * np.pi / lam
    return 8 * np.pi / 3 * k ** 4 * (alpha_au * A0 ** 3) ** 2


def save_png(path, fig):
    fig.savefig(path, dpi=80, bbox_inches="tight")
