"""zz_gfactor.py - follow-up of zz_free.py F5 (written AFTER its third run, information only).

F5 found g(walk) = 1.82 against g(continuum Dirac) = 1.989 for the same sigma = 8 packet on the theta = 0.5
lattice. Hypothesis (before this run): it is the lattice's kinetic-mass error m_kin = m (1 + theta^2/3)
found in F1 (the spin current is ~ 1/m_kin), so g(walk) -> g(continuum) like theta^2.
Here: sigma = 3 packet in a box of 24 lambda_C, theta = a / lambda_C = 0.5, 0.25, 0.125.
Run: timeout 590 python lab/zz_gfactor.py
"""
import json, os, time
import numpy as np
from zz_core import Walk, bloch_batch, gaussian3, dirac_velocity

OUT = os.path.join(os.path.dirname(__file__), "out")
SIGZ4 = np.diag([1, -1, 1, -1]).astype(complex)
t0 = time.time()


def g_factor(sigma, n, a, walk):
    W = Walk("dirac", n, a, theta=a)
    psi = gaussian3(W, sigma, [1, 0, 1, 0])
    pk = np.fft.fftn(psi, axes=(1, 2, 3)).reshape(4, -1)
    del psi
    kk = 2 * np.pi * np.fft.fftfreq(n, d=a)
    KX, KY, KZ = [q.ravel() for q in np.meshgrid(kk, kk, kk, indexing="ij")]
    wgt = np.sum(np.abs(pk)**2, axis=0)
    sel = np.nonzero(wgt > 1e-14 * wgt.max())[0]
    out = np.zeros_like(pk)
    for b0 in range(0, len(sel), 200000):
        sb = sel[b0:b0 + 200000]
        if walk:
            w, v = np.linalg.eig(bloch_batch(np.column_stack([KX[sb], KY[sb], KZ[sb]]) * a, a, 0.5))
            o = np.argsort(np.angle(w), axis=1)[:, :2]          # most negative angle = positive frequency
            Q, _ = np.linalg.qr(np.take_along_axis(v, o[:, None, :], axis=2))
            c = np.einsum("nab,na->nb", Q.conj(), pk[:, sb].T)
            out[:, sb] = np.einsum("nab,nb->na", Q, c).T
        else:
            E = np.sqrt(KX[sb]**2 + KY[sb]**2 + KZ[sb]**2 + 1)
            ph = pk[:, sb]
            def sdotk(s):
                return np.array([KZ[sb] * s[0] + (KX[sb] - 1j * KY[sb]) * s[1], (KX[sb] + 1j * KY[sb]) * s[0] - KZ[sb] * s[1]])
            Hph = np.concatenate([-sdotk(ph[0:2]) + ph[2:4], sdotk(ph[2:4]) + ph[0:2]])
            out[:, sb] = 0.5 * (ph + Hph / E)
    psi = np.fft.ifftn(out.reshape(4, n, n, n), axes=(1, 2, 3))
    psi /= np.sqrt(np.sum(np.abs(psi)**2))
    v, rho = dirac_velocity(psi)
    x = W.x1
    X, Y, _ = np.meshgrid(x, x, x, indexing="ij")
    mu = 0.5 * np.sum(X * v[1] * rho - Y * v[0] * rho)
    S = 0.5 * np.real(np.sum(np.conj(psi) * np.tensordot(SIGZ4, psi, axes=(1, 0))))
    return 2 * mu / S


res = {}
gc = g_factor(3.0, 48, 0.5, False)
print(f"continuum Dirac, sigma = 3: g = {gc:.5f}", flush=True)
res["continuum"] = gc
for th, n in ((0.5, 48), (0.25, 96), (0.125, 192)):
    g = g_factor(3.0, n, th, True)
    res[str(th)] = g
    print(f"walk theta = {th}: g = {g:.5f}; g/g_cont - 1 = {g / gc - 1:+.4f}; -(theta^2/3)/(1+theta^2/3) = {-(th**2 / 3) / (1 + th**2 / 3):+.4f}  ({time.time() - t0:.0f} s)", flush=True)
json.dump(res, open(os.path.join(OUT, "zz_gfactor.json"), "w"), indent=1)
