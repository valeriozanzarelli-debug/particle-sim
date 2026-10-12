"""colla_modes.py JOB K [--wall S] - C2b: lowest normal modes of a relaxed AN Skyrmion (colla_ckpt_JOB.pt),
generalised eigenproblem  K xi = omega^2 (2M) xi  on the tangent bundle (K = Hessian of E incl. the Lagrange term,
M = kinetic metric). LOBPCG (scipy) with an FFT preconditioner P (-Delta + m^2 + 0.2)^-1 / 2, resumable in chunks
(block X saved in lab/out/colla_modes_ckpt_JOB.npz). Writes lab/out/colla_modes_JOB.json when converged
(residual |K x - w2 2M x| / (|2M x| max(w2, 0.1 m^2)) < 1e-2 for all requested modes, or iteration cap 150)."""
import json, os, sys, time, math, argparse, numpy as np, torch
from scipy.sparse.linalg import LinearOperator, lobpcg
import colla_core as cc, knot_core as kc, knot2_core as k2
torch.set_num_threads(int(os.environ.get('COLLA_THREADS', '2')))
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
H = 0.2
ap = argparse.ArgumentParser(); ap.add_argument("job"); ap.add_argument("k", type=int); ap.add_argument("--wall", type=float, default=500)
a = ap.parse_args()
t0 = time.time()
AN = json.load(open(os.path.join(OUT, "colla_units.json")))["AN"]; m = AN["m"]
phi = torch.load(os.path.join(OUT, f"colla_ckpt_{a.job}.pt"), weights_only=False)["phi"]
shape = phi.shape; n = phi.numel()
K = cc.hvp_factory(phi, H, m=m)
g = k2.derivs_fb(phi, H)
CBIG = 50.0


def T(x):
    return torch.from_numpy(np.ascontiguousarray(x)).reshape(shape)


def Aop(X):
    X = np.atleast_2d(X.T).T if X.ndim == 1 else X
    out = np.empty_like(X)
    for j in range(X.shape[1]):
        x = T(X[:, j]); xt = kc.tangent(phi, x)
        out[:, j] = (K(xt) + CBIG * (x - xt)).detach().numpy().ravel()
    return out


def Bop(X):
    X = np.atleast_2d(X.T).T if X.ndim == 1 else X
    out = np.empty_like(X)
    for j in range(X.shape[1]):
        x = T(X[:, j]); xt = kc.tangent(phi, x)
        out[:, j] = (2 * cc.mass_apply(phi, xt, H, 1.0, g) + (x - xt)).numpy().ravel()
    return out


L = shape[1]
q = 2 * np.pi * np.fft.fftfreq(L, d=H)
lap = (2 - 2 * np.cos(q * H)) / H**2
D = lap[:, None, None] + lap[None, :, None] + lap[None, None, :] + m * m + 0.2


def Pop(X):
    X = np.atleast_2d(X.T).T if X.ndim == 1 else X
    out = np.empty_like(X)
    for j in range(X.shape[1]):
        x = X[:, j].reshape(shape)
        y = np.real(np.fft.ifftn(np.fft.fftn(x, axes=(1, 2, 3)) / D[None], axes=(1, 2, 3))) / 2
        out[:, j] = kc.tangent(phi, torch.from_numpy(y)).numpy().ravel()
    return out


A = LinearOperator((n, n), matvec=Aop, matmat=Aop, dtype=np.float64)
Bm = LinearOperator((n, n), matvec=Bop, matmat=Bop, dtype=np.float64)
Pm = LinearOperator((n, n), matvec=Pop, matmat=Pop, dtype=np.float64)
ck = os.path.join(OUT, f"colla_modes_ckpt_{a.job}.npz")
if os.path.exists(ck):
    d = np.load(ck); X = d["X"]; its = int(d["its"])
else:
    rng = np.random.default_rng(3)
    X = Pop(rng.standard_normal((n, a.k))); its = 0
log = open(os.path.join(OUT, "colla_logs", f"modes_{a.job}.log"), "a")
done = False
while time.time() - t0 < a.wall - 120 and not done:
    lam, X = lobpcg(A, X, B=Bm, M=Pm, largest=False, maxiter=15, tol=1e-30, verbosityLevel=0)
    its += 15
    R = Aop(X) - Bop(X) * lam[None, :]
    rel = np.linalg.norm(R, axis=0) / (np.linalg.norm(Bop(X), axis=0) * np.maximum(lam, 0.1 * m * m))
    s = f"its {its}: omega^2/m^2 = {np.round(lam / m**2, 4).tolist()}  max rel res {rel.max():.2e}  t {time.time()-t0:.0f}s"
    print(s, flush=True); log.write(s + "\n"); log.flush()
    np.savez(ck, X=X, its=its)
    done = rel.max() < 1e-2 or its >= 150
if done:
    # mode character: overlap of each eigenvector with translation/rotation generators
    A_, S_, c = cc.generators(phi, H)
    tr = [kc.tangent(phi, d_) for d_ in kc.derivs(phi, H)]
    gens = [x / math.sqrt(float((x * x).sum())) for x in tr + A_ + S_]
    ov = []
    for j in range(X.shape[1]):
        v = kc.tangent(phi, T(X[:, j])); v = v / v.norm()
        ov.append([float(abs((v * gg).sum())) for gg in gens])
    res = dict(job=a.job, m=m, iterations=its, omega2=lam.tolist(), omega2_over_m2=(lam / m**2).tolist(),
               rel_res=rel.tolist(), overlap_trans_iso_spin=ov)
    json.dump(res, open(os.path.join(OUT, f"colla_modes_{a.job}.json"), "w"), indent=1)
    print("DONE", flush=True)
