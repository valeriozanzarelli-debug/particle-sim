"""supercrit_run_g2.py - exams G2 and D1: the dived 1s let go in time.
Initial state: the quasi-bound eigenstate of the rule in a box closed at the middle of the Coulomb barrier
(r_w), embedded in a large box (r < 70) with an absorbing sponge for r > 40. Measured: P_in(t) = weight inside
r_w (decays like exp(-Gamma t)); the wave at r = 30 (outgoing negative-energy electron = the emitted positron,
charge +e); its spectrum (peak = resonance energy, width = Gamma). Lattice a = 0.002 (G2 config).
Writes lab/out/supercrit_g2.json and lab/out/supercrit_g2_snaps.npz (not in git). Usage: [Z ...]"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supercrit_core import SWalk, V_shell, ALPHA, quasi_bound, evolve_probe, sponge, fit_breit_wigner
from supercrit_prereg import P

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "supercrit_g2.json")
res = json.load(open(OUT)) if os.path.exists(OUT) else {}
g1 = json.load(open(os.path.join(HERE, "out", "supercrit_g1.json")))
a = P["configs"]["G2"]["a"]
Zlist = [int(x) for x in sys.argv[1:]] or P["configs"]["G2"]["Z"]
RMAX, RS0, RS1, RPROBE = 70.0, 40.0, 70.0, P["configs"]["D1"]["detector_r"]
t0 = time.time()
for Z in Zlist:
    if str(Z) in res:
        continue
    Za = Z * ALPHA
    Vf = lambda r: V_shell(r, Za, 0.031)
    ref = g1[str(Z)][f"{a:.6g}"]
    Gph, Eph = ref["Gamma"], ref["E_res"]
    om, chi_box, rw = quasi_bound(lambda rmax: SWalk(-1, a, rmax, Vf), ref["E_box"])
    w = SWalk(-1, a, RMAX, Vf)
    chi0 = np.zeros((2, w.N), complex)
    chi0[:, :chi_box.shape[1]] = chi_box
    chi0 /= np.sqrt(np.sum(np.abs(chi0)**2))
    mask = sponge(w.r, RS0, RS1, strength=0.02)
    T = min(25.0 / Gph, 2500.0)
    nt = int(T / a)
    nsnap = 12
    chi, rec, probe, snaps = evolve_probe(w, chi0, nt, mask, rw, RPROBE, rec_every=50,
                                          snap_every=max(nt // (4 * nsnap), 1), nsnap=nsnap)
    t = rec[:, 0] * a
    Pin = rec[:, 1]
    sel = (Pin < 0.7) & (Pin > 0.05)
    slope = np.polyfit(t[sel], np.log(Pin[sel]), 1)[0]
    Gt = -slope
    # spectrum of the outgoing wave at r = 30 (g component); quasi-energy omega: chi ~ exp(-i omega t)
    g = (probe[:, 0] + probe[:, 1]) / np.sqrt(2)
    tt = np.arange(nt) * a
    sel_t = tt > (RPROBE / 1.0) * 0.5
    om_grid = np.linspace(Eph - 8 * Gph, Eph + 8 * Gph, 801)
    # direct Fourier sum on the grid (fine resolution), A(omega) = sum g(t) exp(+i omega t)
    gs, ts = g[sel_t][::5], tt[sel_t][::5]
    A = np.array([np.sum(gs * np.exp(1j * o * ts)) for o in om_grid])
    S = np.abs(A)**2
    k = int(np.argmax(S))
    half = S[k] / 2
    above = np.where(S > half)[0]
    fwhm = om_grid[above[-1]] - om_grid[above[0]]
    # Lorentzian fit
    from scipy.optimize import curve_fit
    L = lambda x, A0, x0, G, c: A0 * (G / 2)**2 / ((x - x0)**2 + (G / 2)**2) + c
    pL, _ = curve_fit(L, om_grid, S, p0=[S[k], om_grid[k], max(fwhm, 1e-9), 0.0])
    # charge bookkeeping before absorption: P_in + outside (r < 40) + absorbed = 1
    res[str(Z)] = dict(a=a, r_w=rw, E_box=om, Gamma_phase=Gph, E_phase=Eph, Gamma_time=Gt,
                       ratio_time_phase=Gt / Gph, spec_peak=float(pL[1]), spec_fwhm=float(abs(pL[2])),
                       spec_fwhm_raw=float(fwhm), T=T, nticks=nt, P_in_end=float(Pin[-1]),
                       norm_end=float(rec[-1, 2]))
    print(Z, res[str(Z)], f"{time.time() - t0:.0f}s", flush=True)
    np.savez_compressed(os.path.join(HERE, "out", f"supercrit_g2_{Z}.npz"), t=t, Pin=Pin, tot=rec[:, 2],
                        snaps=snaps, r=w.r, snap_dt=max(nt // (4 * nsnap), 1) * a, om=om_grid, S=S)
    json.dump(res, open(OUT, "w"), indent=1)
