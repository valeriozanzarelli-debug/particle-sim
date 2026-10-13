"""emiss_ritz.py - exact levels of the 3-D rule (48^3, learned nucleus) by filter diagonalization, with the BLIND
criterion of the v4 pre-registration (lab/out/emiss_prereg2.json).
Run: python lab/emiss_ritz.py <state> [learned|bare]  -> lab/out/emiss_ck/level4_48_<state>.npy, lab/out/emiss_levels4_48.json"""
import sys, os, json, time
import numpy as np
from emiss_core import CK, OUT, norm2, ovl
from emiss_v4_dev import make_atom, start_state, S_LEARNED

NT = 1200
TOL = 0.0015


def ritz(D, psi0, lo, hi, step=0.005):
    a = D.a
    oms = np.arange(lo, hi + 1e-9, step)
    acc = np.zeros((len(oms),) + psi0.shape, complex)
    p = psi0.copy()
    for t in range(NT + 1):
        if t % 2 == 0:
            h = np.sin(np.pi * t / NT)**2
            ph = h * np.exp(1j * oms * t * a)
            for k in range(len(oms)):
                acc[k] += ph[k] * p
        if t < NT:
            D.tick(p, t)
    # orthonormal basis of the filtered waves
    A = acc.reshape(len(oms), -1)
    G = A.conj() @ A.T
    ev, V = np.linalg.eigh(G)
    keep = ev > 1e-10 * ev.max()
    Q = (V[:, keep] / np.sqrt(ev[keep])).T @ A          # rows orthonormal
    # two-tick rule on the basis
    UQ = np.empty_like(Q)
    for i in range(Q.shape[0]):
        v = Q[i].reshape(psi0.shape).copy()
        D.tick(v, 0)
        D.tick(v, 1)
        UQ[i] = v.ravel()
    H = Q.conj() @ UQ.T                                 # H_ij = <q_i|U2|q_j>
    lam, C = np.linalg.eig(H)
    vecs = C.T @ Q                                      # Ritz vectors (rows)
    out = []
    for l, v in zip(lam, vecs):
        v = v / np.linalg.norm(v)
        w = -np.angle(l) / (2 * a)
        Uv = v.reshape(psi0.shape).copy()
        D.tick(Uv, 0); D.tick(Uv, 1)
        res = np.linalg.norm(Uv.ravel() - l * v)
        wt = abs(np.vdot(v, psi0.ravel()))**2
        out.append((w, res, wt, v))
    return out


if __name__ == "__main__":
    st = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else "learned"
    t0 = time.time()
    D = make_atom(s=S_LEARNED if mode == "learned" else (1.0, 1.0))
    psi0 = start_state(D, st)
    psi0 /= np.sqrt(norm2(psi0))
    lo, hi = (0.82, 0.94) if st.startswith("1s") else (0.94, 1.01)
    R = ritz(D, psi0, lo, hi)
    R = [r for r in R if lo - 0.01 < r[0] < hi + 0.01]
    best = max(R, key=lambda r: r[2])
    clus = [r for r in R if abs(r[0] - best[0]) < TOL]
    phi = sum(r[3] * np.vdot(r[3], psi0.ravel()) for r in clus).reshape(psi0.shape)
    phi /= np.sqrt(norm2(phi))
    # Rayleigh frequency and residual of the chosen state
    Up = phi.copy(); D.tick(Up, 0); D.tick(Up, 1)
    l = np.vdot(phi.ravel(), Up.ravel())
    w_state = -np.angle(l) / (2 * D.a)
    res_state = np.linalg.norm(Up.ravel() - l * phi.ravel())
    tag = f"level4_48_{st.replace('/', '_')}" + ("" if mode == "learned" else "_bare")
    np.save(os.path.join(CK, tag + ".npy"), phi)
    jf = os.path.join(OUT, "emiss_levels4_48.json")
    js = json.load(open(jf)) if os.path.exists(jf) else {}
    js[st + ("" if mode == "learned" else "_bare")] = dict(
        omega=float(w_state), residual=float(res_state), weight_in_start=float(abs(ovl(phi, psi0))**2),
        cluster=[dict(omega=float(r[0]), residual=float(r[1]), weight=float(r[2])) for r in clus],
        ritz=sorted([dict(omega=float(r[0]), residual=float(r[1]), weight=float(r[2])) for r in R if r[2] > 1e-3], key=lambda d: -d["weight"])[:10],
        seconds=time.time() - t0)
    json.dump(js, open(jf, "w"), indent=1)
    print(f"{st} ({mode}): omega {w_state:.6f}, residual {res_state:.1e}, weight in start {js[st + ('' if mode == 'learned' else '_bare')]['weight_in_start']:.4f}, "
          f"cluster {[(round(r[0], 5), round(r[2], 3)) for r in clus]}  ({time.time() - t0:.0f} s)", flush=True)
