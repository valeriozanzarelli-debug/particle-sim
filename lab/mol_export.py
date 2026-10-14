"""mol_export.py - study mol_: export of a molecule that passed its exams, for the renderer (same spirit as the helium
export /mnt/project-files/ripresa/dati/elio/). Tier "esatto" (ECG): electron positions sampled from |Psi|^2 with
Metropolis (all electrons together, one fixed spin assignment), density in cylindrical coordinates, and the vibrational
levels with infrared transition dipoles (what a photon can do to the molecule) from the curve and dipole function.
Usage: python lab/mol_export.py lih"""
import os, sys, json, itertools, time
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
import numba as nb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mol_ecg as M

DEST = os.environ.get("MOL_DEST", "/mnt/project-files/debiti/dati/molecole")
UP, DN = 0, 1


def perm_amplitudes(chi, N, sigma):
    """A[phi chi](r, sigma) = sum_p sign(p) chi(sigma_p) phi(r_p): list of (p, weight) with weight != 0"""
    chi = chi.reshape((2,) * N)
    out = []
    for p in itertools.permutations(range(N)):
        sgn = np.linalg.det(np.eye(N)[list(p)])
        w = sgn * chi[tuple(sigma[i] for i in p)]
        if abs(w) > 1e-12:
            out.append((p, w))
    return np.array([o[0] for o in out], np.int64), np.array([o[1] for o in out])


@nb.njit(cache=True)
def psi_batch(r, A, sz, c, perms, wts):
    """r (nw, N, 3); A (K, N, N); sz (K, N) z centres; returns Psi (nw,)"""
    nw, N = r.shape[0], r.shape[1]
    K = A.shape[0]
    out = np.zeros(nw)
    d = np.zeros((N, 3))
    for a in range(nw):
        tot = 0.0
        for q in range(perms.shape[0]):
            for k in range(K):
                for i in range(N):
                    pi = perms[q, i]
                    d[i, 0] = r[a, pi, 0]; d[i, 1] = r[a, pi, 1]; d[i, 2] = r[a, pi, 2] - sz[k, i]
                ex = 0.0
                for i in range(N):
                    for j in range(N):
                        ex += A[k, i, j] * (d[i, 0] * d[j, 0] + d[i, 1] * d[j, 1] + d[i, 2] * d[j, 2])
                tot += wts[q] * c[k] * np.exp(-ex)
        out[a] = tot
    return out


def metropolis(A, sz, c, perms, wts, N, centre_z, rng, nw=800, nstep=1500, burn=500, every=20, step=0.3, log=print,
               init=None):
    """single-electron moves; 10 % of the moves are long jumps (step 1.5) so electrons can change nucleus"""
    if init is None:
        r = rng.normal(0, 1.0, (nw, N, 3)); r[..., 2] += rng.choice(centre_z, (nw, N))
    else:
        r = init(rng, nw)
    p = psi_batch(r, A, sz, c, perms, wts) ** 2
    out = []; acc = 0
    t0 = time.time()
    for it in range(nstep):
        e = rng.integers(N)                         # move one electron at a time
        rn = r.copy(); rn[:, e] += rng.normal(0, step if rng.random() > 0.1 else 1.5, (nw, 3))
        pn = psi_batch(rn, A, sz, c, perms, wts) ** 2
        ok = rng.random(nw) < pn / np.maximum(p, 1e-300)
        r[ok] = rn[ok]; p[ok] = pn[ok]; acc += ok.mean()
        if it >= burn and it % every == 0:
            out.append(r.copy())
        if it % 50 == 0:
            log(f"  step {it} acc {acc / (it + 1):.3f} {time.time() - t0:.0f}s")
    return np.concatenate(out).astype(np.float32), acc / nstep


def export_lih():
    import mol_lih as L
    res = L.load()
    an = res["analysis"]; ver = res["verdict"]
    C = res["curve"]
    Rs = np.array(sorted(v["R"] for v in C.values()))
    Rb = float(Rs[np.argmin(np.abs(Rs - an["R_e"]))])
    z = np.load(L.ck_path(f"R{Rb:.3f}"))
    A, s = z["A"], z["s"]
    sy = L.sys_lih(Rb)
    E, zsum, c = L.dipole_and_energy(sy, A, s)
    # coefficients for normalised basis vs raw: dipole_and_energy works on raw (unnormalised) functions -> c is raw
    chi = M.spin_singlet_pairs(2)
    sigma = (UP, DN, UP, DN)
    perms, wts = perm_amplitudes(chi, 4, sigma)
    rng = np.random.default_rng(2026)
    def init(rng, nw):                              # per spin: one electron on Li (1s), the other on Li or H
        r = rng.normal(0, 0.4, (nw, 4, 3))
        for a in range(nw):
            for pair in ((0, 2), (1, 3)):
                i, j = pair if rng.random() < 0.5 else pair[::-1]
                r[a, j] *= 3.0; r[a, j, 2] += Rb if rng.random() < 0.7 else 0.0
        return r
    S, acc = metropolis(A, s[..., 2].copy(), c, perms, wts, 4, [0.0, Rb], rng, init=init)
    d = os.path.join(DEST, "LiH"); os.makedirs(d, exist_ok=True)
    np.save(os.path.join(d, "lih_elettroni.npy"), S)
    # checks from the samples: <sum z_i> and the dipole
    zs = float(S[..., 2].sum(1).mean())
    mu_s = Rb - zs
    # density in cylindrical coordinates (rho, z), electrons per bohr^3, all 4 electrons
    rr = np.hypot(S[..., 0], S[..., 1]).ravel(); zz = S[..., 2].ravel()
    re_ = np.linspace(0, 4, 81); ze = np.linspace(-3, Rb + 3, 2 * int((Rb + 6) / 0.05) // 2 + 1)
    H, _, _ = np.histogram2d(rr, zz, [re_, ze])
    vol = np.pi * (re_[1:] ** 2 - re_[:-1] ** 2)[:, None] * np.diff(ze)[None]
    dens = H / len(S) / vol
    np.savez_compressed(os.path.join(d, "lih_densita.npz"), rho_edges=re_, z_edges=ze, densita=dens.astype(np.float32))
    # vibrational levels and infrared transitions (photon absorption by the nuclei's motion)
    from scipy.interpolate import CubicSpline
    Es = np.array([C[f"{R:.3f}"]["E"] for R in Rs]); mus = np.array([C[f"{R:.3f}"]["mu_au"] for R in Rs])
    cs = CubicSpline(Rs, Es); cm = CubicSpline(Rs, mus)
    n = 1200; g = np.linspace(Rs[0], Rs[-1], n); h = g[1] - g[0]; ii = np.arange(n)
    Tm = np.where(ii[:, None] == ii[None], np.pi ** 2 / 3, 2.0 / np.where(ii[:, None] == ii[None], 1, (ii[:, None] - ii[None]) ** 2))
    Tm = Tm * (-1.0) ** (ii[:, None] - ii[None]) / (2 * an["mu_red_me"] * h * h)
    w, U = np.linalg.eigh(Tm + np.diag(cs(g)))
    HART_CM = 219474.6313632; AU_D = 2.541746473
    lev = []
    nv = 4
    for v in range(nv):
        lev.append(dict(v=v, E_cm=float((w[v] - w[0]) * HART_CM), mu_D=float(np.sum(U[:, v] ** 2 * cm(g)) * AU_D)))
    trans = []
    for v in range(nv):
        for vp in range(v + 1, nv):
            md = float(np.sum(U[:, v] * U[:, vp] * cm(g)))
            trans.append(dict(da=v, a=vp, nu_cm=float((w[vp] - w[v]) * HART_CM), dipolo_transizione_D=abs(md) * AU_D))
    B_e = 1 / (2 * an["mu_red_me"] * an["R_e"] ** 2) * HART_CM
    man = dict(
        cosa="LiH (idruro di litio), stato fondamentale X1Sigma+, 4 elettroni tutti correlati, per il renderer a bagliori",
        tier="esatto",
        livello="ECG a 4 elettroni con la massa della regola (studio mol_, docs/MOLECOLE.md); nuclei fermi (Born-Oppenheimer)",
        geometria=dict(unita="bohr (a0)", nuclei=[dict(Z=3, nome="Li", pos=[0, 0, 0]), dict(Z=1, nome="H", pos=[0, 0, Rb])],
                       R_e_trovata=an["R_e"], nota="la distanza e' quella del minimo della curva calcolata (punto di griglia piu' vicino)"),
        energia=dict(E_Eh=E, D_e_eV=an["D_e_eV"], D0_eV=an["D0_eV"]),
        file={
            "lih_elettroni.npy": "float32 (n,4,3): n 'fotografie' delle posizioni dei 4 elettroni insieme, da |Psi|^2 con "
                                 "Metropolis. Elettroni 0 e 2 spin su, 1 e 3 spin giu'. Due stanno stretti sul litio (guscio 1s), "
                                 "due formano il legame e stanno soprattutto sull'idrogeno (Li+ H-: il dipolo).",
            "lih_densita.npz": "rho_edges, z_edges, densita (elettroni/bohr^3, tutti e 4) in coordinate cilindriche (rho, z): "
                               "la molecola e' simmetrica attorno all'asse z (stato Sigma)."},
        fotoni=dict(
            nota="Livelli di vibrazione dei nuclei sulla curva calcolata e righe infrarosse (assorbimento di un fotone che fa "
                 "vibrare la molecola). Gli stati elettronici eccitati (A1Sigma+, ~3,3 eV) NON sono ancora calcolati qui.",
            livelli_vibrazione=lev, transizioni_IR=trans, B_e_cm=float(B_e),
            esempio="un fotone di circa {:.0f} cm-1 (lunghezza d'onda {:.1f} micron) porta la molecola da v=0 a v=1".format(
                trans[0]["nu_cm"], 1e4 / trans[0]["nu_cm"])),
        controlli=dict(dipolo_da_campioni_au=mu_s, dipolo_esatto_au=float(Rb - zsum), accettazione=float(acc), n_campioni=int(len(S)),
                       esami={k: dict(PASS=v["PASS"]) for k, v in ver.items()}),
        note=["Stato stazionario: per animare si pescano fotografie diverse a ogni fotogramma.",
              "Debiti: nuclei fermi (Born-Oppenheimer), base ECG scelta a caso e raffinata (un risolutore, non i camminatori della regola), "
              "Coulomb istantaneo e antisimmetria imposti (debiti del thread principale)."],
        script="lab/mol_export.py (repo particle-sim, studio mol_)")
    json.dump(man, open(os.path.join(d, "manifest.json"), "w"), indent=1, ensure_ascii=False)
    print(json.dumps(man["controlli"], indent=1))


if __name__ == "__main__":
    {"lih": export_lih}[sys.argv[1]]()
