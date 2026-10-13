"""nucleo1s_junction3d.py - phase 5, point 2: the live junction between a coarse 3-D zig-zag lattice (cell a, tick a)
and a fine one (cell a/2, tick a/2) across a PLANE x = 0 (transverse directions periodic), far from any nucleus.
The 3-D walk is the split-step Dirac walk of the same family as zz_core.FastDirac3D: per tick coin exp(-i m dt beta),
then the shifts along x, y, z (P_i+ = (1 + alpha_i)/2 moves +1 cell, P_i- moves -1 cell).
Junction (per coarse transverse cell and per coarse tick): 1 coarse right-mover (2-dim P_x+ space) + 8 fine left-mover
slots (2 x 2 transverse fine cells x 2 fine ticks, each 2-dim P_x-) in; 1 coarse left-mover + 8 fine right-mover slots
out: a 9x9 (x 2 spin) unitary, built as in the 1-D prototype:
    smooth fine mode  s = sum_slots conj(e_slot) f_slot,  e = (1/sqrt8) * (1 for the first fine tick, p for the second),
    p = exp(-i m a/2);   [g_s, cL] = cos(th) [cR, s] - i sin(th) beta [s, cR],  th = m a/2 (the mass rotation missed
    while the fine slots are collected: no knob);   rough part f - e s  -> reflected into the fine side through beta.
Only the faces are tested (no edges/corners of a closed shell). Units hbar = c = m = 1."""
import os, sys, json, time
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex)
Z2 = np.zeros((2, 2)); I2 = np.eye(2)
BETA = np.block([[I2, Z2], [Z2, -I2]]).astype(complex)
ALP = [np.block([[Z2, s], [s, Z2]]) for s in (s1, s2, s3)]
PP = [(np.eye(4) + A) / 2 for A in ALP]
PM = [(np.eye(4) - A) / 2 for A in ALP]


def coin(th):
    return np.cos(th) * np.eye(4) - 1j * np.sin(th) * BETA


def apply(M, psi):
    return np.einsum("ij,...j->...i", M, psi)


def shift_T(psi, ax):
    """transverse shift (periodic) along axis ax = 1 (y) or 2 (z)"""
    i = ax
    return np.roll(apply(PP[i], psi), 1, axis=ax) + np.roll(apply(PM[i], psi), -1, axis=ax)


def run(a=0.25, k=(0.4, 0.0, 0.0), sig=12.0, Ny=1, Nz=1, m=1.0, thj=None, steps=None, nside=30):
    thj = m * a / 2 if thj is None else thj
    kx, ky, kz = k
    Cc, Cf = coin(m * a), coin(m * a / 2)
    # positive-energy Bloch vector of the coarse walk W(k) = S_z S_y S_x C
    S = [PP[i] * np.exp(-1j * k[i] * a) + PM[i] * np.exp(1j * k[i] * a) for i in range(3)]
    W = S[2] @ S[1] @ S[0] @ Cc
    ev, V = np.linalg.eig(W)
    om = -np.angle(ev) / a
    E = np.sqrt(m**2 + kx**2 + ky**2 + kz**2)
    v = V[:, np.argmin(np.abs(om - E))]
    Nxc = int(nside * sig / a); Nxf = int(nside * sig / (a / 2))
    x = (np.arange(Nxc) - Nxc + 0.5) * a
    y = np.arange(Ny) * a; z = np.arange(Nz) * a
    x0 = -nside * sig / 3
    env = np.exp(-(x - x0)**2 / (4 * sig**2) + 1j * kx * x)
    psi_c = (env[:, None, None, None] * np.exp(1j * ky * y)[None, :, None, None] * np.exp(1j * kz * z)[None, None, :, None]
             * v[None, None, None, :])
    psi_c /= np.sqrt(np.sum(np.abs(psi_c)**2))
    psi_f = np.zeros((Nxf, 2 * Ny, 2 * Nz, 4), complex)
    p = np.exp(-1j * m * a / 2)
    cj, sj = np.cos(thj), np.sin(thj)
    coll = []
    inj = [np.zeros((2 * Ny, 2 * Nz, 4), complex)] * 2
    vg = kx / E
    steps = steps or int((-x0 + 8 * sig) / vg / (a / 2))
    norms = []
    lost = 0.0

    def block(f):          # (2Ny, 2Nz, 4) -> (Ny, Nz, 2, 2, 4): the 4 fine cells under each coarse cell
        return f.reshape(Ny, 2, Nz, 2, 4).transpose(0, 2, 1, 3, 4)

    def unblock(g):
        return g.transpose(0, 2, 1, 3, 4).reshape(2 * Ny, 2 * Nz, 4)
    for t in range(steps):
        if t % 2 == 0:     # coarse tick
            psi_c = apply(Cc, psi_c)
            Pp, Pm = apply(PP[0], psi_c), apply(PM[0], psi_c)
            cR = Pp[-1].copy()
            lost += np.sum(np.abs(Pm[0])**2)            # leaves through the far left edge of the box (counted)
            Pp = np.concatenate([np.zeros_like(Pp[:1]), Pp[:-1]])
            Pm = np.concatenate([Pm[1:], np.zeros_like(Pm[:1])])
            f1, f2 = (coll + [np.zeros((2 * Ny, 2 * Nz, 4), complex)] * 2)[:2]
            F1, F2 = block(f1), block(f2)
            s_ = (F1.sum((2, 3)) + np.conj(p) * F2.sum((2, 3))) / np.sqrt(8)
            gs = cj * cR - 1j * sj * apply(BETA, s_)
            cL = -1j * sj * apply(BETA, cR) + cj * s_
            e1, e2 = 1 / np.sqrt(8), p / np.sqrt(8)
            G1 = e1 * gs[:, :, None, None, :] + apply(BETA, F1 - e1 * s_[:, :, None, None, :])
            G2 = e2 * gs[:, :, None, None, :] + apply(BETA, F2 - e2 * s_[:, :, None, None, :])
            inj = [unblock(G1), unblock(G2)]; coll = []
            Pm[-1] += cL
            psi_c = Pp + Pm
            psi_c = shift_T(shift_T(psi_c, 1), 2)
        psi_f = apply(Cf, psi_f)
        Pp, Pm = apply(PP[0], psi_f), apply(PM[0], psi_f)
        coll.append(Pm[0].copy())
        lost += np.sum(np.abs(Pp[-1])**2)                # leaves through the far right edge of the box (counted)
        Pp = np.concatenate([inj[t % 2][None], Pp[:-1]])
        Pm = np.concatenate([Pm[1:], np.zeros_like(Pm[:1])])
        psi_f = shift_T(shift_T(Pp + Pm, 1), 2)
        if t % 100 == 0 or t == steps - 1:
            buf = sum(np.sum(np.abs(c)**2) for c in coll) + (np.sum(np.abs(inj[1])**2) if t % 2 == 0 else 0)
            norms.append(float(np.sum(np.abs(psi_c)**2) + np.sum(np.abs(psi_f)**2) + buf + lost))
    R = float(np.sum(np.abs(psi_c)**2))
    return dict(a=a, k=list(k), thj=thj, R=R, lost_edges=float(lost), norm_drift=float(max(abs(np.array(norms) - 1))), steps=steps,
                Ny=Ny, Nz=Nz, sig=sig)


if __name__ == "__main__":
    cases = json.load(open(os.path.join(OUT, "nucleo1s_f5_prereg.json")))["config"]["junction"]["cases"]
    outp = os.path.join(OUT, "nucleo1s_junction3d.json")
    res = json.load(open(outp)) if os.path.exists(outp) else []
    done = {json.dumps(r["case"]) for r in res}
    T0 = time.time()
    for c in cases:
        if json.dumps(c) in done:
            continue
        if time.time() - T0 > 60:
            print("time budget: run again"); break
        r = run(a=c["a"], k=tuple(c["k"]), Ny=c["Ny"], Nz=c["Nz"], thj=c.get("thj"), sig=c["sig"])
        r["R_conservative"] = r["R"] + r["lost_edges"]
        r["case"] = c; res.append(r); print(r, flush=True)
        json.dump(res, open(outp, "w"), indent=1)
