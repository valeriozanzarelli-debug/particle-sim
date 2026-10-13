"""emiss_twa.py - the OPTIONS for spontaneous emission, compared on a toy that can be run thousands of times:
a two-level atom (frequency w0, dipole d along z) at the centre of a 1-D zig-zag light line (the light rule along
one axis: right movers u = (E_z - B_y)/sqrt2 move +1 cell per tick, left movers v = (E_z + B_y)/sqrt2 move -1;
exact c along an axis).  Periodic line long enough that nothing comes back during the run.

  QE   one-excitation quantum (Wigner-Weisskopf on the lattice): amplitude b of |e,0>, complex photon field
       sourced by the transition current.  Exact for this model in the rotating-wave approximation.  -> rate A.
  SC   semiclassical (Maxwell-Bloch): pure |e> never decays; a superposition decays as dp/dt = -A p (1-p).
  HTWA hybrid truncated Wigner: SC plus the light's VACUUM fluctuations (each mode of the line sampled with
       energy hbar w/2), the atom kept as a pure two-level wave.  Prediction (rate equations, this file's doc):
       dp/dt = (A/2)(1 - 2p) - A p (1 - p): starts at A/2 (half the true rate), and the GROUND state is pumped
       to p* = 1 - 1/sqrt2 = 0.293.
  STWA spin truncated Wigner (dTWA, Schachenmayer et al. PRX 5, 011022 (2015)): the atom too is sampled (classical
       Bloch vector s with s_z = +-1, s_x, s_y = +-1 random), plus the vacuum noise.  The atom is then no longer a
       wave of the rule but an ensemble of classical spins.

Usage: python lab/emiss_twa.py [n_traj]  ->  lab/out/emiss_twa.json
"""
import sys, json, os, time
import numpy as np
import numba as nb

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

# toy parameters (light cell b = 1, atom sub-steps per light tick)
N = 8192
W0 = 2 * np.pi / 20          # 20 cells per wavelength
SUB = 8
T_TICKS = 3000
ED = 0.06                    # e*d (coupling); sets A


@nb.njit(cache=True)
def run_traj(mode, u0, v0, ce0, cg0, s0, nt, w0, ed, sub, rec_every):
    """mode: 0 QE, 1 SC/HTWA (pure state), 2 STWA (classical spin). Light: real u, v (or complex for QE).
    Atom at site c = N/2.  Field at the atom E = (u + v)/sqrt2 (E_z).  Source: dE_z = -J (J = current/cell)."""
    n = u0.shape[0]
    c = n // 2
    u = u0.copy()
    v = v0.copy()
    un = np.empty_like(u)
    vn = np.empty_like(v)
    ce = ce0
    cg = cg0
    sx, sy, sz = s0[0], s0[1], s0[2]
    s2 = 1.0 / np.sqrt(2.0)
    a = 1.0 / sub
    nrec = nt // rec_every + 1
    rec = np.empty(nrec)
    k = 0
    jprev = 0.0 + 0.0j
    jprev2 = 0.0 + 0.0j
    for t in range(nt):
        if t % rec_every == 0:
            if mode == 0:
                rec[k] = abs(ce)**2
            elif mode == 1:
                rec[k] = abs(ce)**2
            else:
                rec[k] = 0.5 * (1 + sz)
            k += 1
        # free light tick: shift u right, v left
        for i in range(n):
            un[i] = u[(i - 1) % n]
            vn[i] = v[(i + 1) % n]
        u, un = un, u
        v, vn = vn, v
        # predictor of this tick's current
        if abs(jprev2) > 0:
            r = jprev / jprev2
            if abs(r) > 1:
                r = r / abs(r)
            jp = jprev * r
        else:
            jp = jprev
        # self field of a deposit J at site c: dE_z = -J -> E changes by -J (u, v each by -J/sqrt2 -> E by -J)
        E = (u[c] + v[c]) * s2 - 0.5 * jp
        if mode == 0:
            d0 = np.conj(cg) * ce          # transition dipole amplitude (cg = e^{-i0 t} fixed ground, sector)
        elif mode == 1:
            d0 = 2 * (np.conj(cg) * ce).real
        else:
            d0 = sx
        for s in range(sub):
            if mode == 0:
                # vacuum sector |e,0>: free phase, + coupling from the photon sector (ground amplitude 1)
                ce = ce * np.exp(-1j * w0 * a) - 1j * a * ed * E
            elif mode == 1:
                # exact 2x2 step with H = w0 |e><e| + ed E sigma_x
                Er = E.real
                om = np.sqrt((w0 / 2)**2 + (ed * Er)**2)
                ph = np.exp(-1j * w0 / 2 * a)
                co = np.cos(om * a)
                si = np.sin(om * a) / om if om > 0 else a
                # H' = H - w0/2 = (w0/2) sz + ed E sx  (sz = +1 on e)
                ne = (co - 1j * (w0 / 2) * si) * ce - 1j * ed * Er * si * cg
                ng = (co + 1j * (w0 / 2) * si) * cg - 1j * ed * Er * si * ce
                ce, cg = ne * ph, ng * ph
            else:
                Er = E.real
                # ds/dt = Omega x s, Omega = (2 ed E, 0, w0)
                ox, oz = 2 * ed * Er, w0
                om = np.sqrt(ox * ox + oz * oz)
                nx, nz = ox / om, oz / om
                ang = om * a
                co, si = np.cos(ang), np.sin(ang)
                dot = nx * sx + nz * sz
                cx = 0 * sx - nz * sy           # n x s with n = (nx, 0, nz)
                cy = nz * sx - nx * sz
                cz = nx * sy
                sx, sy, sz = (sx * co + cx * si + nx * dot * (1 - co),
                              sy * co + cy * si,
                              sz * co + cz * si + nz * dot * (1 - co))
        if mode == 0:
            d1 = np.conj(cg) * ce
        elif mode == 1:
            d1 = 2 * (np.conj(cg) * ce).real
        else:
            d1 = sx
        J = -(d1 - d0) * ed        # e * d(dipole)/dt over one light tick (b = 1); STWA: <r> = d s_x
        if mode == 2:
            J = J * 1.0
        dE = -0.5 * (J - jp)          # corrector: the self field of this tick's actual deposit
        if mode == 0:
            ce = ce - 1j * ed * dE
        elif mode == 1:
            dr = dE.real
            ne = ce - 1j * ed * dr * cg
            ng = cg - 1j * ed * dr * ce
            nn = np.sqrt(abs(ne)**2 + abs(ng)**2)
            ce, cg = ne / nn, ng / nn
        else:
            ang = 2 * ed * dE.real
            co, si = np.cos(ang), np.sin(ang)
            sy, sz = sy * co - sz * si, sz * co + sy * si
        jprev2 = jprev
        jprev = J
        u[c] -= J * s2
        v[c] -= J * s2
    if mode == 2:
        rec[k] = 0.5 * (1 + sz)
    else:
        rec[k] = abs(ce)**2
    return rec


def vacuum(rng, complex_field=False):
    """real right/left movers with energy w_k/2 per mode (light cell b = 1, energy 1/2 sum (u^2 + v^2))"""
    k = 2 * np.pi * np.fft.rfftfreq(N)
    out = []
    for _ in range(2):
        amp = np.sqrt(np.abs(k) / 2 / N)       # <N|alpha_k|^2> * ... -> mode energy w/2
        al = amp * (rng.normal(size=k.size) + 1j * rng.normal(size=k.size)) / np.sqrt(2)
        al[0] = 0
        f = np.fft.irfft(al * N, n=N)
        out.append(f.astype(np.complex128))
    return out[0], out[1]


def traj_seed(mode_tag, i):
    return np.random.default_rng([20261013, mode_tag, i])


if __name__ == "__main__":
    ntraj = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    rec_every = 10
    ckf = os.path.join(OUT, "emiss_ck", f"twa_{ntraj}.npz")
    nrec = T_TICKS // rec_every + 1
    t0 = time.time()
    z = np.zeros(N, np.complex128)
    if os.path.exists(ckf):
        ck = dict(np.load(ckf))
    else:
        ck = {"HTWA_e": np.zeros(nrec), "HTWA_g": np.zeros(nrec), "STWA_e": np.zeros(nrec), "STWA_g": np.zeros(nrec),
              "done": np.zeros(4, int)}
    for mi, (key, tag) in enumerate((("HTWA_e", 1), ("HTWA_g", 2), ("STWA_e", 3), ("STWA_g", 4))):
        while ck["done"][mi] < ntraj:
            i = int(ck["done"][mi])
            rng = traj_seed(tag, i)
            u, v = vacuum(rng)
            if key.startswith("HTWA"):
                ce0, cg0 = (1.0, 0.0) if key.endswith("e") else (0.0, 1.0)
                ck[key] += run_traj(1, u, v, ce0 + 0j, cg0 + 0j, np.zeros(3), T_TICKS, W0, ED, SUB, rec_every)
            else:
                s0 = np.array([rng.choice([-1, 1]), rng.choice([-1, 1]), 1.0 if key.endswith("e") else -1.0])
                ck[key] += run_traj(2, u, v, 0j, 0j, s0, T_TICKS, W0, ED, SUB, rec_every)
            ck["done"][mi] += 1
            if time.time() - t0 > 540:
                np.savez(ckf, **ck)
                print(f"checkpoint {ck['done']}; run again", flush=True)
                sys.exit(3)
    np.savez(ckf, **ck)
    res = {}
    tt = np.arange(nrec) * rec_every
    rq = run_traj(0, z, z, 1.0 + 0j, 1.0 + 0j, np.zeros(3), T_TICKS, W0, ED, SUB, rec_every)
    m = (tt > 200) & (rq > 0.05)
    A = -np.polyfit(tt[m], np.log(rq[m]), 1)[0]
    res["A_QE"] = A
    res["QE"] = rq.tolist()
    for p0 in (1.0, 0.5):
        res[f"SC_p{p0}"] = run_traj(1, z, z, np.sqrt(p0) + 0j, np.sqrt(1 - p0) + 0j, np.zeros(3), T_TICKS, W0, ED, SUB, rec_every).tolist()
    for key in ("HTWA_e", "HTWA_g", "STWA_e", "STWA_g"):
        res[key] = (ck[key] / ntraj).tolist()
    pe = ck["HTWA_e"] / ntraj
    m1 = (tt >= 50) & (tt <= 300)
    res["T1_htwa_initial_rate_over_A"] = float(-np.polyfit(tt[m1], pe[m1], 1)[0] / A)
    res["T2_htwa_ground_p_end"] = float(ck["HTWA_g"][-1] / ntraj)
    ps = ck["STWA_e"] / ntraj
    m3 = (tt >= 50) & (tt <= 1500)
    res["T3_stwa_rate_over_A"] = float(-np.polyfit(tt[m3], np.log(ps[m3]), 1)[0] / A)
    res["T4_sc_pure_e_dev"] = float(np.max(np.abs(np.array(res["SC_p1.0"]) - 1)))
    res["t"] = tt.tolist()
    res["params"] = dict(N=N, W0=W0, SUB=SUB, T=T_TICKS, ED=ED, ntraj=ntraj)
    # rate-equation prediction for HTWA (no fit): dp/dt = (A/2)(1-2p) - A p (1-p)
    for lab, p0 in (("e", 1.0), ("g", 0.0)):
        p, out = p0, []
        for k in range(nrec):
            out.append(p)
            for _ in range(rec_every):
                p += (A / 2 * (1 - 2 * p) - A * p * (1 - p))
        res[f"HTWA_{lab}_rate_eq"] = out
    json.dump(res, open(os.path.join(OUT, "emiss_twa.json" if ntraj >= 1000 else "emiss_twa_dev.json"), "w"))
    print({k: v for k, v in res.items() if k.startswith("T") or k == "A_QE"})
