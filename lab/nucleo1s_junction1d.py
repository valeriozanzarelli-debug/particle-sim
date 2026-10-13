"""nucleo1s_junction1d.py - DEVELOPMENT prototype (phase 4, point 2; no exam): a live, exactly unitary junction between a
coarse zig-zag lattice (cell a, tick a) and a fine one (cell a/2, tick a/2), in 1-D, as the building block of the 3-D bcc
onion of docs/NUCLEO_1S.md section 13.3.

Channel counting. In one coarse tick the coarse side sends 1 right-mover into the junction and wants 1 left-mover back;
the fine side, in the same time (2 fine ticks), sends 2 left-movers and wants 2 right-movers. 1 + 2 in, 1 + 2 out: a 3x3
unitary J closes it exactly. Choice ("transmit the smooth part, reflect the rest"):
    coarse right-mover -> the 2 fine slots with the expected relative phase  (1, p)/sqrt2
    the fine combination (f1 + p* f2)/sqrt2 -> coarse left-mover
    the orthogonal fine combination (rough, short-wavelength part) -> reflected back into the fine side
with p = exp(-i w_ref a/2), w_ref = m (rest energy, a constant of nature; w_ref = 0 tried, slightly worse), and a
junction coin of angle thj on the pair (coarse right-mover, smooth fine mode): thj = 0 (plain pass) reflects a lot
(R = 0.34 at k = 0.4, a = 0.5: the buffered amplitude misses one fine tick of mass rotation); thj = m a/2 (exactly that
missing rotation, no knob) gives R = 0.008 and R ~ a^4. A one-tick buffer
keeps causality (both fine slots are collected before J is applied). In 3-D the same counting gives 1 coarse channel vs
4 faces x 2 ticks = 8 fine slots per coarse face (J is 9x9 per face and colour).
Units hbar = c = m = 1."""
import numpy as np
import json, os, sys

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def coin(th):
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -1j * s], [-1j * s, c]])


def run(a=0.5, k0=0.3, sig=40.0, wref=1.0, Nc=None, Nf=None, steps=None, m=1.0, x0=None, thj=0.0):
    Nc = Nc or int(40 * sig / a)
    Nf = Nf or int(40 * sig / a)
    Cc, Cf = coin(m * a), coin(m * a / 2)
    # positive-energy right-moving Bloch vector of the coarse walk at k0: W(k) = S(k) C, S = diag(e^{-ik a}, e^{ik a})
    W = np.diag([np.exp(-1j * k0 * a), np.exp(1j * k0 * a)]) @ Cc
    ev, V = np.linalg.eig(W)
    om = -np.angle(ev) / a
    v = V[:, np.argmin(np.abs(om - np.sqrt(m**2 + k0**2)))]
    x = (np.arange(Nc) - Nc + 0.5) * a                       # coarse cells, x < 0
    x0 = -10 * sig if x0 is None else x0
    env = np.exp(-(x - x0)**2 / (4 * sig**2) + 1j * k0 * x)
    Rc, Lc = env * v[0], env * v[1]
    nrm = np.sqrt(np.sum(abs(Rc)**2 + abs(Lc)**2)); Rc /= nrm; Lc /= nrm
    Rf = np.zeros(Nf, complex); Lf = np.zeros(Nf, complex)
    p = np.exp(-1j * wref * a / 2)
    # J: smooth fine mode s = (f1 + p* f2)/sqrt2 and rough mode d = (f1 - p* f2)/sqrt2;
    # [g_s, cL] = C(thj) [cR, s] (a junction coin of angle thj; thj = 0: plain pass), rough mode reflected
    cj, sj = np.cos(thj), np.sin(thj)
    es = np.array([1, p]) / np.sqrt(2); ed = np.array([1, -p]) / np.sqrt(2)
    J = np.zeros((3, 3), complex)
    # input basis vectors: cR, f1, f2 ; s = es^H f, d = ed^H f
    for col, vin in enumerate(np.eye(3)):
        cR_, f = vin[0], vin[1:]
        s_, d_ = np.vdot(es, f), np.vdot(ed, f)
        gs, cL = cj * cR_ - 1j * sj * s_, -1j * sj * cR_ + cj * s_
        J[0, col] = cL; J[1:, col] = gs * es + d_ * ed
    assert np.allclose(J.conj().T @ J, np.eye(3), atol=1e-14)
    coll = []                     # fine left-movers that reached the junction (this coarse tick)
    inj = [0j, 0j]                # fine right-movers to inject in the next two fine ticks
    buf_c = 0j                    # coarse left-mover waiting to enter
    vg = k0 / np.sqrt(m**2 + k0**2)
    steps = steps or int((-x0 + 8 * sig) / vg / (a / 2))

    norms = []
    for t in range(steps):
        if t % 2 == 0:            # coarse tick
            Rc, Lc = Cc[0, 0] * Rc + Cc[0, 1] * Lc, Cc[1, 0] * Rc + Cc[1, 1] * Lc
            cR = Rc[-1]
            Rc = np.concatenate([[0j], Rc[:-1]])
            Lc = np.concatenate([Lc[1:], [0j]])
            f1, f2 = (coll + [0j, 0j])[:2]
            out = J @ np.array([cR, f1, f2])
            Lc[-1] = out[0]; inj = [out[1], out[2]]; coll = []
        Rf, Lf = Cf[0, 0] * Rf + Cf[0, 1] * Lf, Cf[1, 0] * Rf + Cf[1, 1] * Lf
        coll.append(Lf[0])
        Rf = np.concatenate([[inj[t % 2]], Rf[:-1]])
        Lf = np.concatenate([Lf[1:], [0j]])
        if t % 50 == 0 or t == steps - 1:
            buffered = sum(abs(z)**2 for z in coll) + (abs(inj[1])**2 if t % 2 == 0 else 0)
            norms.append(np.sum(abs(Rc)**2 + abs(Lc)**2) + np.sum(abs(Rf)**2 + abs(Lf)**2) + buffered)
    Pref = float(np.sum(abs(Rc)**2 + abs(Lc)**2))
    return dict(a=a, k0=k0, wref=wref, thj=thj, omega=float(np.sqrt(m**2 + k0**2)), R=Pref,
                norm_drift=float(max(abs(np.array(norms) - 1))), norms=norms)


if __name__ == "__main__":
    res = []
    for a in (0.5, 0.25):
        for thj in (0.0, a / 2):           # thj = m a/2: the mass angle of the fine tick spent in the buffer
            for k0 in (0.1, 0.2, 0.4, 0.8):
                r = run(a=a, k0=k0, wref=1.0, sig=25.0, thj=thj)
                r.pop("norms"); res.append(r); print(r, flush=True)
    json.dump(res, open(os.path.join(OUT, "nucleo1s_junction1d.json"), "w"), indent=1)
