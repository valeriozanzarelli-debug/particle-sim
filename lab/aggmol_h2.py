"""aggmol_h2.py - phase 2 of study aggmol_: H2 with two electrons (pre-registration lab/out/aggmol_h2_prereg.json).
  gate    G1: H atom and H2+ with the same ECG machine
  curve   Born-Oppenheimer curve on the pre-registered R grid (checkpoints: lab/out/aggmol_h2.json, bases .npz not committed)
  hf      H5: restricted Hartree-Fock of H2 on the learned a0/4 grid machine (64^3, box +-8 a0)
  vib     R_e, omega_e, vibrational levels, D_e, D_0
Usage: timeout 590 python lab/aggmol_h2.py JOB"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggmol_ecg as E

OUT = E.OUT
JS = os.path.join(OUT, "aggmol_h2.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
RGRID = [0.8, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 2.0, 2.4, 3.0, 4.0, 5.0, 6.0]
EH_EV, EH_CM = 27.211386245988, 219474.6313632
MP = 1836.15267343
M = E.M_KIN


def save():
    json.dump(res, open(JS, "w"), indent=1, default=float)


def bfile(R):
    return os.path.join(OUT, f"aggmol_h2_basis_R{R:.2f}.npz")


def job_gate():
    import aggmol_exact as ex
    out = {}
    e, B, t = E.build(1, 0.0, 14, ntrial=40, nref=6, budget=200)
    out["H_atom"] = dict(E=e, exact=-M / 2, err=e + M / 2)
    e, B, t = E.build(1, 2.0, 40, ntrial=40, nref=6, budget=300)
    key = M * ex.energy(M * 2.0)[0]
    out["H2plus_R2"] = dict(E_el=e, exact=key, err=e - key)
    out["PASS"] = bool(abs(out["H_atom"]["err"]) <= 1e-6 and abs(out["H2plus_R2"]["err"]) <= 1e-6)
    res["G1"] = out
    save()
    print(out, flush=True)


def job_curve():
    t0 = time.time()
    cur = res.setdefault("curve", {})
    order = [1.4] + [r for r in RGRID if r < 1.4][::-1] + [r for r in RGRID if r > 1.4]
    for R in order:
        key = f"{R:.2f}"
        if key in cur:
            continue
        if time.time() - t0 > 330:
            print("checkpoint: run again", flush=True); return
        if R == 1.4:
            e, B, t = E.build(2, R, 180, ntrial=40, nref=8, budget=480)
        else:
            Rn = 1.4 if R in (1.35, 1.45) else [r for r in order[:order.index(R)] if abs(r - R) == min(abs(q - R) for q in order[:order.index(R)])][0]
            st = np.load(bfile(Rn))
            s = st["s"] * (R / Rn)
            e, B, t = E.build(2, R, 180, ntrial=40, nref=4, budget=220, state=dict(A=st["A"].copy(), s=s))
        np.savez(bfile(R), A=B.A, s=B.s)
        cur[key] = dict(R=R, E_el=e, E_tot=e + 1 / R, K=len(B.A), sec=t)
        save()
        print(f"R={R}: E_tot {e + 1 / R:.9f} Eh ({t:.0f} s)", flush=True)
    print("CURVE DONE", flush=True)


def job_hf():
    import aggmol_core as mc
    import scipy.fft as sfft
    B = 1024; h = B * mc.TH; N = 64; R = 1.4
    rec = mc.level_rec(B); C, offs, g0 = mc.coupling_rec(B)
    T = mc.kinetic_grid(rec, N, h)
    cs = [(0, 0, R * mc.A0 / 2 / h), (0, 0, -R * mc.A0 / 2 / h)]
    fn = mc.MultiField(N, h, cs)
    Vn = {s: fn.V[s].copy() for s in fn.subs}
    M2 = 2 * N
    k = np.rint(np.fft.fftfreq(N) * N).astype(int) % M2
    # Coulomb kernel on a 4N grid (zero padding of the 2N grid, cutoff R_c = box length): isolated Hartree potential
    M4 = 2 * M2; hh = h / 2
    q1 = 2 * np.pi * np.fft.fftfreq(M4, d=hh)
    QX, QY, QZ = np.meshgrid(q1, q1, q1, indexing="ij", sparse=True)
    q2 = QX ** 2 + QY ** 2 + QZ ** 2
    Rc = N * h
    with np.errstate(divide="ignore", invalid="ignore"):
        Kq = np.where(q2 > 0, 4 * np.pi * (1 - np.cos(np.sqrt(q2) * Rc)) / q2, 2 * np.pi * Rc ** 2)
    del q2, QX, QY, QZ

    def hartree(phi):
        F = sfft.fftn(phi)
        G = np.zeros((M2, M2, M2), complex); G[np.ix_(k, k, k)] = F
        up = np.real(sfft.ifftn(G)) * 8                      # phi on the 2N grid (exact band-limited interpolation)
        n = up ** 2 / h ** 3                                 # density of one electron per volume (lambda_C^-3)
        big = np.zeros((M4,) * 3); big[:M2, :M2, :M2] = n
        # Kq is the continuous Fourier transform of the kernel, so the volume element is already inside it
        # (bug found after the first production run: an extra hh^3 made the Hartree potential 5026x too big; declared)
        V = np.real(sfft.ifftn(sfft.fftn(big) * Kq))[:M2, :M2, :M2] * mc.ALPHA
        # the 2N grid starts at coarse index 0: sub-lattice s = V[s0::2, s1::2, s2::2]; the density sits in the corner
        # [0, M2) of the padded box, i.e. the potential is evaluated on the same points
        return {s: V[s[0]::2, s[1]::2, s[2]::2].copy() for s in fn.subs}

    t0 = time.time()
    ck = os.path.join(OUT, "aggmol_h2_hf_phi.npy")           # checkpoint (not committed)
    st = res.get("H5_run", {})
    if os.path.exists(ck) and st:
        phi = np.load(ck); Eold, hist = st["Eold"], st["hist"]
    else:
        phi = mc.guess_mol(N, h, cs, "1sg").reshape((N,) * 3); phi /= np.linalg.norm(phi)
        Eold, hist = 0.0, []
    for it in range(len(hist), 40):
        if time.time() - t0 > 380:
            np.save(ck, phi); res["H5_run"] = dict(Eold=Eold, hist=hist); save()
            print("checkpoint: run again", flush=True); return
        VH = hartree(phi)
        fn.V = {s: Vn[s] + VH[s] for s in fn.subs}
        mol = mc.Mol(T, fn, gamma0=g0.real, sym=(1, 1, 1))
        X, Es, r, conv, its = mc.solve_until(mol, phi.ravel()[:, None], 120, rtol=1e-9, log=lambda *a: None)
        phi_new = X[:, 0].reshape((N,) * 3); phi_new /= np.linalg.norm(phi_new)
        if np.sum(phi_new * phi) < 0:
            phi_new = -phi_new
        phi = phi_new
        # energy: 2 <h> + J,   <h> = <f> - J  ->  E = 2 eps - J
        fn.V = VH
        J = float(np.sum(phi * fn.apply(phi).real))
        eps = Es[0]
        Etot = (2 * eps - J) / mc.ALPHA ** 2 + 1 / R
        hist.append(Etot)
        print(f"  SCF {it}: E_HF = {Etot:.9f} Eh (eps {eps / mc.ALPHA ** 2:.6f}, J {J / mc.ALPHA ** 2:.6f})", flush=True)
        if abs(Etot - Eold) < 1e-8:
            break
        Eold = Etot
    res.pop("H5_run", None)
    res["H5"] = dict(E_HF_grid=Etot, literature_HF_limit=-1.133629, err_mEh=1e3 * (Etot + 1.133629), history=hist,
                     eps_Eh=eps / mc.ALPHA ** 2, J_Eh=J / mc.ALPHA ** 2, PASS=bool(abs(Etot + 1.133629) <= 1.6e-3))
    save()
    print(res["H5"], flush=True)


def job_vib():
    from scipy.interpolate import CubicSpline
    from scipy.linalg import eigh_tridiagonal
    cur = res["curve"]
    Rs = np.array(sorted(v["R"] for v in cur.values()))
    Et = np.array([cur[f"{R:.2f}"]["E_tot"] for R in Rs])
    EH = -M / 2
    sel = (Rs >= 1.2) & (Rs <= 1.7)
    c = np.polyfit(Rs[sel] - 1.4, Et[sel], 4)
    rr = np.roots(np.polyder(c)); rr = rr[np.abs(rr.imag) < 1e-9].real + 1.4; Re = float(rr[(rr > 1.2) & (rr < 1.7)][0])
    Emin = float(np.polyval(c, Re - 1.4)); k = float(np.polyval(np.polyder(c, 2), Re - 1.4))
    mu = MP / 2
    we = np.sqrt(k / mu)
    De = 2 * EH - Emin
    # 1-D nuclear problem (J = 0) on the spline; beyond 6 a0 the curve is flat to within 1e-4 Eh (asymptote 2 E_H)
    cs = CubicSpline(Rs, Et)
    x = np.linspace(0.8, 6.0, 6000); dx = x[1] - x[0]
    V = cs(x)
    d = 1 / (mu * dx ** 2) + V; o = -0.5 / (mu * dx ** 2) * np.ones(len(x) - 1)
    lev = eigh_tridiagonal(d, o, select="i", select_range=(0, 3))[0]
    nu01 = lev[1] - lev[0]; zpe = lev[0] - Emin
    out = dict(R_e=Re, E_min=Emin, k_e=k, omega_e_cm=we * EH_CM, omega_e_eV=we * EH_EV, D_e_eV=De * EH_EV,
               nu01_cm=nu01 * EH_CM, ZPE_eV=zpe * EH_EV, D_0_eV=(2 * EH - lev[0]) * EH_EV, levels_Eh=lev.tolist())
    ex = dict(R_e=1.401, D_e_eV=4.747, omega_e_cm=4401.0, nu01_cm=4161.0, D_0_eV=4.478)
    E14 = cur["1.40"]["E_tot"]
    V = {
        "H0": dict(E=E14, err_mEh=1e3 * (E14 + 1.174475714), PASS=bool(abs(E14 + 1.174475714) <= 5e-5)),
        "H1": dict(R_e=Re, rel=Re / 1.401 - 1, PASS=bool(abs(Re / 1.401 - 1) <= 3e-3)),
        "H2": dict(D_e_eV=De * EH_EV, diff_meV=1e3 * (De * EH_EV - 4.747), PASS=bool(abs(De * EH_EV - 4.747) <= 0.0434)),
        "H3": dict(omega_e_cm=we * EH_CM, rel=we * EH_CM / 4401 - 1, PASS=bool(abs(we * EH_CM / 4401 - 1) <= 0.01)),
        "H4": dict(nu01_cm=nu01 * EH_CM, rel=nu01 * EH_CM / 4161 - 1, D_0_eV=out["D_0_eV"], D0_diff_meV=1e3 * (out["D_0_eV"] - 4.478),
                   PASS=bool(abs(nu01 * EH_CM / 4161 - 1) <= 0.01 and abs(out["D_0_eV"] - 4.478) <= 0.0434)),
    }
    res["vib"] = out; res["verdict"] = V
    save()
    print(json.dumps(out, indent=1)); print(json.dumps(V, indent=1))


if __name__ == "__main__":
    dict(gate=job_gate, curve=job_curve, hf=job_hf, vib=job_vib)[sys.argv[1]]()
