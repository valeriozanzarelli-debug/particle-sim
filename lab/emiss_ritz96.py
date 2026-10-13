"""emiss_ritz96.py - exact levels by filter diagonalization for large boxes (memory-lean, checkpointed), v5 settings
(lab/out/emiss_prereg3.json).  Same blind criterion as emiss_ritz.py.
Run (repeat until exit 0): python lab/emiss_ritz96.py <state> [n] [NT] [step] [lo hi]
 -> lab/out/emiss_ck/level5_<n>_<state>.npy, lab/out/emiss_levels5.json"""
import sys, os, json, time
import numpy as np
from emiss_core import CK, OUT, norm2, ovl
from emiss_v4_dev import make_atom, start_state, S_LEARNED

TOL = 0.0015
t0 = time.time()
st = sys.argv[1]
n = int(sys.argv[2]) if len(sys.argv) > 2 else 96
NT = int(sys.argv[3]) if len(sys.argv) > 3 else 1200
step = float(sys.argv[4]) if len(sys.argv) > 4 else 0.005
if len(sys.argv) > 6:
    lo, hi = float(sys.argv[5]), float(sys.argv[6])
else:
    lo, hi = (0.88, 0.92) if st.startswith("1s") else (0.955, 0.995)
tag = f"{n}_{st.replace('/', '_')}_{NT}_{step}"
ckf = os.path.join(CK, f"ritz5_{tag}.npz")
D = make_atom(n=n, s=S_LEARNED)
psi0 = start_state(D, st)
psi0 /= np.sqrt(norm2(psi0))
a = D.a
oms = np.arange(lo, hi + 1e-9, step)
if os.path.exists(ckf):
    z = np.load(ckf)
    acc, p, tst = z["acc"], z["p"], int(z["t"])
    print(f"resumed at tick {tst}")
else:
    acc = np.zeros((len(oms),) + psi0.shape, complex)
    p = psi0.copy()
    tst = 0
for t in range(tst, NT + 1):
    if t % 2 == 0:
        h = np.sin(np.pi * t / NT)**2
        ph = h * np.exp(1j * oms * t * a)
        for k in range(len(oms)):
            acc[k] += ph[k] * p
    if t < NT:
        D.tick(p, t)
        if time.time() - t0 > 500 and t + 1 < NT:
            np.savez(ckf, acc=acc, p=p, t=t + 1)
            print(f"checkpoint at tick {t + 1}/{NT}; run again", flush=True)
            sys.exit(3)
del p
A = acc.reshape(len(oms), -1)
G = A.conj() @ A.T
ev, V = np.linalg.eigh(G)
keep = ev > 1e-10 * ev.max()
Q = (V[:, keep] / np.sqrt(ev[keep])).T @ A
del acc, A
m = Q.shape[0]
H = np.empty((m, m), complex)
for j in range(m):
    v = Q[j].reshape(psi0.shape).copy()
    D.tick(v, 0)
    D.tick(v, 1)
    H[:, j] = Q.conj() @ v.ravel()
lam, C = np.linalg.eig(H)
C = C / np.linalg.norm(C, axis=0)
c0 = Q.conj() @ psi0.ravel()                    # <q_i|psi0>
wts = np.abs(C.conj().T @ c0)**2
oml = -np.angle(lam) / (2 * a)
order = np.argsort(-wts)
best = order[0]
clus = [i for i in range(m) if abs(oml[i] - oml[best]) < TOL]
vecs = {i: (C[:, i] @ Q).reshape(psi0.shape) for i in set(clus) | set(order[:4])}
phi = sum(vecs[i] * np.vdot(vecs[i].ravel(), psi0.ravel()) for i in clus)
phi /= np.sqrt(norm2(phi))
Up = phi.copy(); D.tick(Up, 0); D.tick(Up, 1)
l = np.vdot(phi.ravel(), Up.ravel())
w_state, res_state = -np.angle(l) / (2 * a), np.linalg.norm(Up.ravel() - l * phi.ravel())
for i in order[:4]:
    np.save(os.path.join(CK, f"ritzvec5_{tag}_{i}.npy"), vecs[i])
np.save(os.path.join(CK, f"level5_{n}_{st.replace('/', '_')}" + ("" if (NT, step) == (1200, 0.005) else f"_{NT}_{step}") + ".npy"), phi)
jf = os.path.join(OUT, "emiss_levels5.json")
js = json.load(open(jf)) if os.path.exists(jf) else {}
js[tag] = dict(omega=float(w_state), residual=float(res_state), weight_in_start=float(abs(ovl(phi, psi0))**2),
               cluster=[dict(omega=float(oml[i]), weight=float(wts[i])) for i in clus],
               top=[dict(index=int(i), omega=float(oml[i]), weight=float(wts[i])) for i in order[:6]], seconds=time.time() - t0)
json.dump(js, open(jf, "w"), indent=1)
if os.path.exists(ckf):
    os.remove(ckf)
print(f"{tag}: omega {w_state:.6f}, residual {res_state:.1e}, cluster {[(round(oml[i], 5), round(wts[i], 3)) for i in clus]}, "
      f"top {[(round(oml[i], 5), round(wts[i], 3)) for i in order[:4]]}  ({time.time() - t0:.0f} s)", flush=True)
