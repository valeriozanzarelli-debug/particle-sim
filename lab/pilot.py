"""Rule R4 (pilot wave, all local):
  wave balls  — fly straight, turn phase 2π per λ, re-emitted at slit edges (Huygens)
  photon ball — one of the balls; its sideways velocity is read from the
                phases of the neighbours around it:
                  v = Im(ψ* ∂xψ) / (k |ψ|²),   ψ = Σ_neigh a_j e^{iφ_j}
Photon absorbed ⇔ the photon ball hits the wall. Click ⇔ where it lands.
"""
import numpy as np, time
from rules import fresnel_double_slit, K, LAM

rng = np.random.default_rng(4)
W0, D, A, L = 7.0, 8.0, 2.0, 160.0
TH = 0.45
BINS, HALF = 240, 60.0


def wave_after_slits(n_emit):
    """Balls of one photon that get through the slits, as Huygens rays."""
    x = rng.normal(0, W0 / 2, n_emit)
    inside = (np.abs(x - D / 2) <= A / 2) | (np.abs(x + D / 2) <= A / 2)
    s = x[inside]
    th = rng.uniform(-TH, TH, s.size)
    amp = 0.5 * (1 + np.cos(th)) * np.exp(s**2 / W0**2)
    return s, th, amp, x


def local_field(xm, z, s, th, amp, h):
    """ψ and ∂xψ at (xm, z) from the neighbour rays within ~h."""
    xr = s + z * np.tan(th)                    # where each ray is at this z
    d = xm[:, None] - xr[None, :]
    w = np.exp(-0.5 * (d / h) ** 2)
    ell = z / np.cos(th)
    r = np.maximum(ell, 1e-3)
    jac = z / np.cos(th) ** 2 if z > 0 else np.ones_like(th)
    c = amp * jac / np.sqrt(LAM * r)
    ph = K * ell[None, :] + K * np.sin(th)[None, :] * d
    e = c[None, :] * w * np.exp(1j * ph)
    psi = e.sum(1)
    dpsi = (e * (1j * K * np.sin(th))[None, :]).sum(1)
    return psi, dpsi


def run(n_photons=400, n_emit=30000, h=1.0, dz=0.5, share_wave=False):
    edges = np.linspace(-HALF, HALF, BINS + 1)
    xs = 0.5 * (edges[1:] + edges[:-1])
    clicks, absorbed = [], 0
    wave = wave_after_slits(n_emit) if share_wave else None
    t = time.time()
    for p in range(n_photons):
        s, th, amp, x_all = wave if share_wave else wave_after_slits(n_emit)
        # the photon is one ball, picked at random among all balls
        xm0 = x_all[rng.integers(x_all.size)]
        if not ((abs(xm0 - D / 2) <= A / 2) or (abs(xm0 + D / 2) <= A / 2)):
            absorbed += 1
            continue
        xm = np.array([xm0])
        z = 1e-3
        while z < L:
            hz = h * max(1.0, 0.0)
            psi, dpsi = local_field(xm, z, s, th, amp, hz)
            v = np.imag(np.conj(psi) * dpsi) / (K * np.abs(psi) ** 2 + 1e-30)
            xm = xm + dz * v
            z += dz
        clicks.append(xm[0])
    ref = fresnel_double_slit(xs, W0, D, A, L)
    hist, _ = np.histogram(clicks, edges)
    g = 8
    hc = hist.reshape(-1, g).sum(1); rc = ref.reshape(-1, g).sum(1)
    corr = np.corrcoef(hc, rc)[0, 1]
    return corr, len(clicks), absorbed, time.time() - t, hist, ref


if __name__ == "__main__":
    for h in (0.5, 1.0, 2.0):
        corr, n, ab, dt, *_ = run(n_photons=300, n_emit=20000, h=h, share_wave=True)
        print(f"R4 pilot, window h={h}: clicks {n}, absorbed {ab}, corr (coarse) with Fresnel {corr:.3f}  [{dt:.0f}s]", flush=True)
