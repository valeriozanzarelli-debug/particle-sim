"""zstar_film.py - side task of ZSTAR: export the giro-1 film (SUPERCRITICO.md, G2/D1, Z = 200, shell nucleus
R = 0.031 lambda_C, a = 0.002) frame by frame for a slow-motion video.

Same run as lab/supercrit_run_g2.py (same initial state: the rule's quasi-bound state in a box closed at the middle
of the Coulomb barrier; big box r < 70, sponge beyond r = 40), recorded on a radial grid of bins:
  rho[k, b]  radial probability density of the electron (per lambda_C, already x 4 pi r^2: integral over r = weight)
  j[k, e]    radial probability current through the sphere at edge e (units c), AVERAGED over the frame interval
             (exact lattice flux: per tick |chi+_i|^2 - |chi-_{i+1}|^2 through the link i -> i+1)
  Q[k, e]    cumulative weight inside the edge radius e: Q(r, t) = int_0^r rho  (continuity: dQ/dt = -j exactly)
Charge reading (hole theory, SUPERCRITICO.md 5.5): the state is the dived 1s; what flows out is the emitted
positron (+e), what stays is the charged vacuum (-e). Exported numbers are electron PROBABILITIES (multiply by -e).

Writes lab/out/zstar_film.json (compact: values with 4 significant digits) and PNG previews in
/mnt/project-files/problemi-aperti/zstar_film_*.png.  Usage: python zstar_film.py   (~2-4 min, 1 core)
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supercrit_core import SWalk, V_shell, ALPHA, quasi_bound, sponge, njit
from zz_core import MC2_EV

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "zstar_film.json")
PNGDIR = "/mnt/project-files/problemi-aperti"
Z, R_NUC, A = 200, 0.031, 0.002
RMAX, RS0, RS1 = 70.0, 40.0, 70.0
T_FILM, NFRAMES = 240.0, 300


@njit(cache=True)
def _film(chip, chim, M, refl0, refl1, mask, nticks, every, bin_of, nb, edge_link, ne, rho_out, j_out):
    N = chip.shape[0]
    p = np.empty(N, np.complex128)
    q = np.empty(N, np.complex128)
    acc = np.zeros(ne)
    k = 0
    for t in range(nticks):
        # flux through the chosen links during this tick's shift
        for e in range(ne):
            i = edge_link[e]
            acc[e] += (chip[i].real**2 + chip[i].imag**2) - (chim[i + 1].real**2 + chim[i + 1].imag**2)
        p[1:] = chip[:-1]
        p[0] = refl0 * chim[0]
        q[:-1] = chim[1:]
        q[N - 1] = refl1 * chip[N - 1]
        for i in range(N):
            a_ = M[i, 0, 0] * p[i] + M[i, 0, 1] * q[i]
            b_ = M[i, 1, 0] * p[i] + M[i, 1, 1] * q[i]
            chip[i] = a_ * mask[i]
            chim[i] = b_ * mask[i]
        if (t + 1) % every == 0:
            for b in range(nb):
                rho_out[k, b] = 0.0
            for i in range(N):
                b = bin_of[i]
                if b >= 0:
                    rho_out[k, b] += chip[i].real**2 + chip[i].imag**2 + chim[i].real**2 + chim[i].imag**2
            for e in range(ne):
                j_out[k, e] = acc[e] / every
                acc[e] = 0.0
            k += 1
    return chip, chim


def r4(x):
    """round to 4 significant digits (float32 precision is ~7: this is the 'compact' part)"""
    x = np.asarray(x, dtype=np.float64)
    out = np.zeros_like(x)
    nz = x != 0
    mag = np.floor(np.log10(np.abs(x[nz])))
    out[nz] = np.round(x[nz] / 10**(mag - 3)) * 10**(mag - 3)
    return [float(f"{v:.4g}") for v in out.ravel()]


if __name__ == "__main__":
    t0 = time.time()
    g1 = json.load(open(os.path.join(HERE, "out", "supercrit_g1.json")))
    ref = g1[str(Z)][f"{A:.6g}"]
    Za = Z * ALPHA
    Vf = lambda r: V_shell(r, Za, R_NUC)
    om, chi_box, rw = quasi_bound(lambda rmax: SWalk(-1, A, rmax, Vf), ref["E_box"])
    w = SWalk(-1, A, RMAX, Vf)
    chi0 = np.zeros((2, w.N), complex)
    chi0[:, :chi_box.shape[1]] = chi_box
    chi0 /= np.sqrt(np.sum(np.abs(chi0)**2))
    mask = sponge(w.r, RS0, RS1, strength=0.02)
    # radial bins: fine near the nucleus (the dived 1s has size ~ 0.1-0.3 lambda_C), coarser out to the sponge
    edges = np.concatenate([np.arange(0.0, 2.0, 0.025), np.arange(2.0, 40.0 + 1e-9, 0.125)])
    nb = len(edges) - 1
    bin_of = np.searchsorted(edges, w.r, side="right") - 1
    bin_of[(w.r >= edges[-1])] = -1
    # link i -> i+1 crossing each inner edge (edge between r_i and r_{i+1})
    edge_link = np.clip(np.searchsorted(w.r, edges[1:]) - 1, 0, w.N - 2)
    nticks = int(round(T_FILM / A))
    every = nticks // NFRAMES
    rho = np.zeros((NFRAMES, nb))
    jj = np.zeros((NFRAMES, nb))
    chip, chim = chi0[0].copy(), chi0[1].copy()
    # frame 0 = initial state
    rho0 = np.zeros(nb)
    np.add.at(rho0, bin_of[bin_of >= 0], (np.abs(chip)**2 + np.abs(chim)**2)[bin_of >= 0])
    _film(chip, chim, w.M, w.refl0, w.refl1, mask, nticks, every, bin_of, nb, edge_link.astype(np.int64), nb, rho, jj)
    jj /= A                                          # per tick -> per unit time (tick = a / c)
    rho = np.vstack([rho0[None], rho])
    jj = np.vstack([np.zeros((1, nb)), jj])
    dt_frame = every * A
    t = np.arange(NFRAMES + 1) * dt_frame
    dr = np.diff(edges)
    Q = np.cumsum(rho, axis=1)                       # weight inside each outer edge
    dens = rho / dr[None, :]                         # per lambda_C
    # continuity check: Q(edge, t_k) - Q(edge, t_{k-1}) = - j * dt_frame  (j averaged over the interval)
    cont = np.max(np.abs((Q[1:] - Q[:-1]) + jj[1:] * dt_frame))
    i_w = int(np.searchsorted(edges, rw))
    i30 = int(np.searchsorted(edges, 30.0))
    P_in = Q[:, i_w - 1]
    passed30 = np.concatenate([[0.0], np.cumsum(jj[1:, i30 - 1] * dt_frame)])
    sel = (P_in < 0.7) & (P_in > 0.05)
    Gamma_t = -np.polyfit(t[sel], np.log(P_in[sel]), 1)[0]
    print(f"Z={Z} a={A}: frames {NFRAMES + 1} every {dt_frame:.3f}, r_w = {rw:.3f}, continuity err {cont:.1e}, "
          f"Gamma(time, pre-reg window) = {Gamma_t:.5f}  (phase {ref['Gamma']:.5f})  {time.time() - t0:.0f}s")
    man = {
        "formato": "zstar_film v1 (JSON): stesse parole dell'export dell'atomo (lab/atomo_export_leggi.py su "
                   "origin/master: manifest + campi per fotogramma + scalari per fotogramma), ma 1-D radiale e in "
                   "un solo file JSON; valori arrotondati a 4 cifre significative",
        "studio": "SUPERCRITICO.md 5.5 (G2/D1), esportato da ZSTAR.md", "Z": Z, "nucleo": "guscio R = 0.031 lambda_C (12 fm)",
        "regola": "zz_core RadialWalk (kappa = -1, s1/2) con V a guscio; passo a = 0.002 lambda_C",
        "unita": {"r": "lambda_C = 386.16 fm", "t": "lambda_C / c = 1.288e-21 s", "rho": "probabilita' per lambda_C "
                  "(gia' integrata sull'angolo solido)", "j": "probabilita' per (lambda_C/c) attraverso la sfera di "
                  "raggio r (= c x densita' di flusso radiale)", "Q": "probabilita' entro r", "E": "m c^2 = 511 keV"},
        "carica": "numeri dell'elettrone (probabilita'): carica = -e x numero. Il peso che esce e' il positrone "
                  "emesso (+e), quello che resta il vuoto carico (-e)",
        "fotogrammi": NFRAMES + 1, "dt_fotogramma": dt_frame, "tick_per_fotogramma": every,
        "griglia": {"r_bordi": "len nb+1", "r_centri": "len nb", "nb": nb},
        "campi": {"rho": {"forma": [NFRAMES + 1, nb], "dove": "centri"},
                  "j": {"forma": [NFRAMES + 1, nb], "dove": "bordi esterni r_bordi[1:]", "nota": "media sull'intervallo (t_{k-1}, t_k]; j[0] = 0"},
                  "Q": {"forma": [NFRAMES + 1, nb], "dove": "bordi esterni r_bordi[1:]"}},
        "risonanza": {"E_r": ref["E_res"], "Gamma": ref["Gamma"], "E_box": float(om), "r_w": float(rw),
                      "Gamma_dal_tempo": float(Gamma_t), "Gamma_eV": ref["Gamma"] * MC2_EV},
        "spugna": "assorbe per r > 40 (fuori dalla griglia esportata)",
        "controllo_continuita_max": float(cont),
    }
    data = {"manifest": man, "r_bordi": r4(edges), "r_centri": r4(0.5 * (edges[1:] + edges[:-1])), "t": r4(t),
            "rho": [r4(x) for x in dens], "j": [r4(x) for x in jj], "Q": [r4(x) for x in Q],
            "fotogrammi_scalari": {"t": r4(t), "P_in_rw": r4(P_in), "Q_entro_40": r4(Q[:, -1]),
                                   "passato_r30": r4(passed30)}}
    s = json.dumps(data, separators=(",", ":"))
    open(OUT, "w").write(s)
    print("wrote", OUT, f"{len(s) / 1e6:.2f} MB")
    # previews
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(PNGDIR, exist_ok=True)
    rc = 0.5 * (edges[1:] + edges[:-1])
    picks = [0, int(round(10 / dt_frame)), int(round(30 / dt_frame)), int(round(60 / dt_frame)),
             int(round(120 / dt_frame)), NFRAMES]
    for n, k in enumerate(picks):
        fig, ax = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
        ax[0].semilogy(rc, np.maximum(dens[k], 1e-9), color="#2e86c1")
        ax[0].set_ylim(1e-7, 20)
        ax[0].set_ylabel("ρ(r) [1/λ_C]")
        ax[1].plot(rc, jj[k], color="#c0392b")
        ax[1].set_ylabel("j(r) [c]  (flusso radiale)")
        ax[1].axhline(0, color="k", lw=0.5)
        ax[2].plot(rc, Q[k], color="#27ae60")
        ax[2].set_ylim(0, 1.02)
        ax[2].set_ylabel("Q(<r)")
        ax[2].set_xlabel("r [λ_C]")
        ax[2].set_xscale("log")
        ax[2].set_xlim(0.02, 40)
        ax[0].set_title(f"Z = {Z}: tuffo dell'1s, t = {t[k]:.1f} λ_C/c  (P entro r_w = {P_in[k]:.3f}; "
                        f"1/Γ = {1 / ref['Gamma']:.1f})")
        fig.tight_layout()
        fig.savefig(os.path.join(PNGDIR, f"zstar_film_{n + 1}.png"), dpi=90)
        plt.close(fig)
    print("previews done", f"{time.time() - t0:.0f}s")
