import numpy as np, time
from pilot import *
import pilot

def batch(n_photons, n_emit, h, dz, waves=8):
    edges = np.linspace(-HALF, HALF, BINS + 1); xs = 0.5*(edges[1:]+edges[:-1])
    clicks=[]; absorbed=0; t=time.time()
    per = n_photons // waves
    for wv in range(waves):                      # a fresh wave every `per` photons
        s, th, amp, x_all = wave_after_slits(n_emit)
        xm = x_all[pilot.rng.integers(x_all.size, size=per)]
        ok = (np.abs(xm-D/2)<=A/2)|(np.abs(xm+D/2)<=A/2)
        absorbed += (~ok).sum(); xm = xm[ok]
        z = 1e-3
        while z < L:
            for i0 in range(0, xm.size, 200):
                sl = slice(i0, i0+200)
                psi, dpsi = local_field(xm[sl], z, s, th, amp, h)
                xm[sl] += dz * np.imag(np.conj(psi)*dpsi)/(K*np.abs(psi)**2+1e-30)
            z += dz
        clicks += list(xm)
    ref = fresnel_double_slit(xs, W0, D, A, L)
    hist,_ = np.histogram(clicks, edges)
    hc = hist.reshape(-1,8).sum(1); rc = ref.reshape(-1,8).sum(1)
    return np.corrcoef(hc, rc)[0,1], len(clicks), absorbed, time.time()-t, hist, ref, xs

if __name__ == "__main__":
    import sys
    for h in (0.5, 2.0):
        c, n, ab, dt, hist, ref, xs = batch(8000, 20000, h, 0.5)
        print(f"h={h}: {n} clicks, {ab} absorbed ({100*ab/(n+ab):.0f}%), corr with Fresnel {c:.3f} [{dt:.0f}s]", flush=True)
        np.save(f"hist_h{h}.npy", hist); np.save("ref.npy", ref)
