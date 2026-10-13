"""dilat_ione.py - numeri per un esperimento con uno ione che simula Dirac (Gerritsma et al., Nature 463, 68 (2010)).
Nessun esame: previsione. Parametri dello ione dal lavoro di Gerritsma: eta = 0.06, Omega~ = 2 pi 68 kHz,
c_sim = 2 eta Omega~ x0 (x0 = dimensione dello stato fondamentale), mc^2 = hbar Omega con Omega <= 2 pi 13 kHz.
Calcolo: continuo di Dirac (stesso motore della regola) + la regola a a = 1/16 come controllo.
Scrive lab/out/dilat_ione.json."""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dilat_core import (two_packets, u_plus_dirac, u_plus_lattice, spinor_state, overlap_spec)

HERE = os.path.dirname(os.path.abspath(__file__))
eta, Omt, Om = 0.06, 2 * np.pi * 68e3, 2 * np.pi * 13e3
c_sim = 2 * eta * Omt                      # in x0 / s
mc = Om / c_sim                            # in hbar / x0
T_unit = 1 / Om                            # s
sig_p = 0.5 / mc                           # stato fondamentale: sigma_p = hbar / (2 x0), in unita' mc
Dl = np.sqrt(2) * sig_p                    # convenzione Smith-Ahmadi
print("mc = %.3f hbar/x0, unita' di tempo %.1f us, sigma_p = %.3f mc, Dl = %.3f mc" % (mc, T_unit * 1e6, sig_p, Dl))

p = 2 * np.pi * np.fft.fftfreq(2 ** 15, d=2048 / 2 ** 15)
a = 1 / 16
k = 2 * np.pi * np.fft.fftfreq(int(1024 / a), d=a)
rows = []
for p1, p2, al, ph in [(0.5, 0.5 + Dl, np.pi / 8, np.pi), (0.5, 0.5 + Dl, np.pi / 8, 0.0), (0.0, Dl, np.pi / 4, np.pi),
                       (0.5, 0.5 + 2 * Dl, np.pi / 8, np.pi)]:
    def beta(g, kk, u_fun, kind):
        # <beta> per uno stato di energia positiva = ritmo dell'orologio (Hellmann-Feynman), da delta -> 0
        u = u_fun
        O = overlap_spec(spinor_state(g, u), kk, [1.0 if kind == "dirac" else 16], 1.0, 1e-6, kind, a if kind != "dirac" else None)[0]
        return -np.angle(O) / 1e-6
    g, g1, g2 = two_packets(p, np.cos(al), np.sin(al), ph, p1, p2, Dl)
    u = u_plus_dirac(p, 1.0)
    bs, b1, b2 = beta(g, p, u, "dirac"), beta(g1, p, u, "dirac"), beta(g2, p, u, "dirac")
    gl, gl1, gl2 = two_packets(k, np.cos(al), np.sin(al), ph, p1, p2, Dl)
    ul = u_plus_lattice(k, a, 1.0)
    bl = beta(gl, k, ul, "lattice") - (np.cos(al) ** 2 * beta(gl1, k, ul, "lattice") + np.sin(al) ** 2 * beta(gl2, k, ul, "lattice"))
    dq = bs - (np.cos(al) ** 2 * b1 + np.sin(al) ** 2 * b2)
    # visibilita' con echo di Loschmidt (avanti con Omega + d/2, indietro con Omega - d/2), d = 0.5 mc^2, t = 5..20 unita'
    ts = [5.0, 10.0, 20.0]
    dl = 0.5
    Vs = np.abs(overlap_spec(spinor_state(g, u), p, ts, 1.0, dl, "dirac"))
    O1 = overlap_spec(spinor_state(g1, u), p, ts, 1.0, dl, "dirac")
    O2 = overlap_spec(spinor_state(g2, u), p, ts, 1.0, dl, "dirac")
    Vm = np.abs(np.cos(al) ** 2 * O1 + np.sin(al) ** 2 * O2)
    N1 = (1 / dq) ** 2 * 3 if dq != 0 else np.inf     # colpi per 1 sigma ~ (rumore di proiezione <= 1)/sqrt(N), 3 misure
    rows.append(dict(p1=p1, p2=p2, alpha=al, phi=ph, beta_sup=bs, beta_1=b1, beta_2=b2, dq=dq, dq_regola_a16=bl,
                     colpi_per_5sigma=25 * N1, t_us=[t * T_unit * 1e6 for t in ts], V_sup_echo=Vs.tolist(),
                     V_mix_echo=Vm.tolist(), V2_sup=(Vs ** 2).tolist(), V2_mix=(Vm ** 2).tolist()))
    print("p=(%.2f,%.2f) al=%.3f ph=%.2f  <beta> sup %.4f mix %.4f  dQ %.4f (regola %.4f)  colpi 5sig %.0f  V2 sup %s mix %s" % (
        p1, p2, al, ph, bs, np.cos(al) ** 2 * b1 + np.sin(al) ** 2 * b2, dq, bl, 25 * N1,
        np.round(Vs ** 2, 4), np.round(Vm ** 2, 4)))
json.dump(dict(parametri=dict(eta=eta, Omega_tilde_Hz=Omt / 2 / np.pi, Omega_Hz=Om / 2 / np.pi, mc_hbar_su_x0=mc,
                              unita_tempo_us=T_unit * 1e6, sigma_p_mc=sig_p, Dl_mc=Dl), righe=rows),
          open(os.path.join(HERE, "out", "dilat_ione.json"), "w"), indent=1, default=float)
