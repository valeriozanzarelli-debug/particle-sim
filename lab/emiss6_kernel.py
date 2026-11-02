"""emiss6_kernel.py - study "Emissione vera", phase 4 (docs/EMISSIONE.md §13): the RESPONSE of the zig-zag light rule
at the atom, computed once by running the light rule itself.

Why.  In the coupled runs the atom talks to the light only through two 3-vectors per light tick: it READS the
Gaussian-averaged field at its centre (Light3D.field_at_atom) and it WRITES its current there (transverse deposit,
Light3D.deposit_T).  The light rule (walk + sponge + deposit) is LINEAR, so the field the atom reads at light tick t is
exactly  E(t) = sum_{d < t} K_{p(d)}(t - d) . J_d  with K the response of the light rule to one unit deposit.  The rule
alternates the axis order every tick, so it repeats every TWO ticks: there are two kernels, one per parity p of the
deposit tick.  Field and deposit scale as 1/b^2 (b = light cell): K is computed at b = 1 and rescaled.
This is the same rule, rearranged (like the "two boxes" of ATOMO.md §3.8); its only approximation is the cut of the
kernel when it has died away (sponge), checked by exam G0 against the full box (lab/emiss6_g0.py).

Output: lab/out/emiss_ck/emiss6_kernel_<n>_<sigma>.npz  (K[p, k, i, c] = field component i at age k from unit J_c)
Run: python lab/emiss6_kernel.py [n=128] [sigma=1] [ages=700]
"""
import os, sys, time
import numpy as np
from emiss_core import Light3D, CK


def kernel(n=128, sigma=1.0, ages=700, sponge=12):
    f = os.path.join(CK, f"emiss6_kernel_{n}_{sigma:g}.npz")
    if os.path.exists(f):
        return np.load(f)["K"]
    L = Light3D(n, 1.0, sponge=sponge, sigma=sigma)
    L.make_transverse()
    K = np.zeros((2, ages + 1, 3, 3))
    t0 = time.time()
    for p in (0, 1):
        for c in range(3):
            L.F[...] = 0
            L.tick_count = p              # the deposit happens right after the tick with this count
            I = np.zeros(3); I[c] = 1.0
            L.deposit_T(I)
            L.sponge()
            for k in range(1, ages + 1):
                L.tick()
                K[p, k, :, c] = L.field_at_atom()
                L.sponge()
            print(f"  parity {p} column {c}: |K(1)| {np.abs(K[p, 1, :, c]).max():.3e}  |K({ages})| "
                  f"{np.abs(K[p, ages, :, c]).max():.3e}  ({time.time() - t0:.0f} s)", flush=True)
    np.savez(f, K=K, M=L.M, n=n, sigma=sigma, sponge=sponge)
    return K


def one_tick_rot(K, nt):
    """one-tick light matrices S_nt(k) (as zz_core.bloch_batch, spin 1, without the product of the two ticks)"""
    from zz_core import BASIS_ONE, order_of_tick
    N = K.shape[0]
    T = np.broadcast_to(np.eye(3, dtype=complex), (N, 3, 3)).copy()
    for j in order_of_tick(nt):
        lam, P = BASIS_ONE[j]
        ph = np.exp(-1j * lam[None, :] * K[:, j:j + 1])
        T = np.einsum("ab,nb,cb->nac", P, ph, P.conj()) @ T
    return T


def transverse_odd(n=128, sigma=1.0):
    """AMENDMENT (phase 4): deposit patterns for deposits made at ODD light-tick counts.  The light rule repeats every
    two ticks, and its static (longitudinal) mode at an odd count is the even one carried by one tick: n_odd = S0 n_even.
    The production projector uses n_even at every count (Light3D.make_transverse), so half of the deposits leave a
    piece of static field behind.  T_odd_c = S0 P_even S0^T (w e_c)."""
    f = os.path.join(CK, f"transverse_odd_{n}_{sigma:g}.npy")
    if os.path.exists(f):
        return np.load(f)
    from zz_core import bloch_batch
    L = Light3D(n, 1.0, sponge=0, sigma=sigma)
    kk = 2 * np.pi * np.fft.fftfreq(n)
    K = np.stack([q.ravel() for q in np.meshgrid(kk, kk, kk, indexing="ij")], -1)
    w = np.zeros((n, n, n)); w[L.idx] = L.w
    wk = np.fft.fftn(w).ravel()
    Pk = np.empty((K.shape[0], 3, 3))
    for b0 in range(0, K.shape[0], 100000):
        sl = slice(b0, b0 + 100000)
        R = bloch_batch(K[sl], spin=1).real
        ax = np.stack([R[:, 2, 1] - R[:, 1, 2], R[:, 0, 2] - R[:, 2, 0], R[:, 1, 0] - R[:, 0, 1]], -1)
        nr = np.linalg.norm(ax, axis=1)
        for i in np.where(nr < 1e-9)[0]:
            w_, v_ = np.linalg.eig(R[i]); j = np.argmin(np.abs(w_ - 1))
            ax[i] = np.real(v_[:, j]); nr[i] = np.linalg.norm(ax[i])
        ax /= nr[:, None]
        S0 = one_tick_rot(K[sl], 0)
        if np.max(np.abs(S0.imag)) > 1e-12:
            raise RuntimeError("one-tick matrix not real")
        S0 = S0.real
        ao = np.einsum("nab,nb->na", S0, ax)
        P = np.eye(3)[None] - ao[:, :, None] * ao[:, None, :]
        if b0 == 0:
            P[0] = 2.0 / 3.0 * np.eye(3)
        Pk[sl] = P
    T = np.empty((3, 3, n, n, n), np.float32)
    for c in range(3):
        Tc = np.fft.ifftn((Pk[:, :, c] * wk[:, None]).T.reshape(3, n, n, n), axes=(1, 2, 3))
        T[c] = Tc.real
    np.save(f, T)
    return T


def kernel_parity(n=128, sigma=1.0, ages=700, sponge=12):
    """the same response with the parity-correct transverse deposit (even counts: production T, odd counts: T_odd)"""
    f = os.path.join(CK, f"emiss6_kernelP_{n}_{sigma:g}.npz")
    if os.path.exists(f):
        return np.load(f)["K"]
    L = Light3D(n, 1.0, sponge=sponge, sigma=sigma)
    L.make_transverse()
    Te = L.T
    To = transverse_odd(n, sigma)
    Mo = np.array([[np.sum(L.w * To[c][cp][L.idx]) for c in range(3)] for cp in range(3)])
    K = np.zeros((2, ages + 1, 3, 3))
    t0 = time.time()
    for p in (0, 1):
        L.T = Te if p == 0 else To
        for c in range(3):
            L.F[...] = 0
            L.tick_count = p
            I = np.zeros(3); I[c] = 1.0
            L.deposit_T(I)
            L.sponge()
            for k in range(1, ages + 1):
                L.tick()
                K[p, k, :, c] = L.field_at_atom()
                L.sponge()
            print(f"  [P] parity {p} column {c}: |K(1)| {np.abs(K[p, 1, :, c]).max():.3e}  |K({ages})| "
                  f"{np.abs(K[p, ages, :, c]).max():.3e}  ({time.time() - t0:.0f} s)", flush=True)
    np.savez(f, K=K, M=L.M, Mo=Mo, n=n, sigma=sigma, sponge=sponge)
    return K


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 128
    s = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    ag = int(sys.argv[3]) if len(sys.argv) > 3 else 700
    K = kernel_parity(n, s, ag) if os.environ.get("PARITA") else kernel(n, s, ag)
    a = np.abs(K).max(axis=(0, 2, 3))
    for k in (1, 2, 5, 10, 20, 50, 100, 150, 200, 300, 400, 500, 600, ag):
        if k <= ag:
            print(f"age {k:4d}: max |K| = {a[k]:.3e}")
