"""ciclo_p1_ledger.py - Phase 1 of CICLO_LIVELLI: the INFLUENCE LEDGER (registro delle influenze).
For each term: weight relative to the dominant one at the scale where the engine has it, scaling exponent y with the
level size (measured = from ciclo_p1_map.json / committed results of other studies, analytic = from dimensions; marked),
and gravity as a tracked term from the start, with its crossover radius against cohesion DERIVED from constants and from
the engine's own solid-H2 numbers (aggmol phase 3 pair curve, read-only). Writes lab/out/ciclo_p1_ledger.json.
Run: python lab/ciclo_p1_ledger.py"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
# constants (SI, CODATA 2018) -- inputs allowed by the rules (physical constants, nuclear masses)
G = 6.67430e-11
HBAR = 1.054571817e-34
C = 299792458.0
ME = 9.1093837015e-31
MP = 1.67262192369e-27
E2 = 2.307077552e-28            # e^2 / (4 pi eps0)  [J m]
ALPHA = 7.2973525693e-3
A0 = 5.29177210903e-11
EH = 4.3597447222071e-18         # hartree [J]
EV = 1.602176634e-19
LC = HBAR / (ME * C)             # reduced Compton length


def fcc_shells(nmax=12):
    """distances (in units of the nearest-neighbour distance) and multiplicities of the fcc lattice"""
    pts = []
    r = range(-nmax, nmax + 1)
    for i in r:
        for j in r:
            for k in r:
                if (i + j + k) % 2 == 0 and (i, j, k) != (0, 0, 0):
                    pts.append(np.sqrt(i * i + j * j + k * k))
    d = np.sqrt(2.0)
    rr = np.round(np.array(pts) / d, 9)
    u, cnt = np.unique(rr, return_counts=True)
    keep = u <= nmax / d * 0.999
    return u[keep], cnt[keep]


def h2_solid():
    """static fcc solid of the ENGINE's isotropic H2-H2 curve (aggmol phase 3, 10 points 5-12 a0, C6/C8 tail fitted
    there): cohesive energy per molecule and volume per molecule at the minimum. No zero-point motion (declared)."""
    d = json.load(open(os.path.join(OUT, "aggmol_h2h2.json")))
    R = np.array([5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0])
    E = np.array(d["iso"]["E_int"])
    C6, C8 = d["tail"]["C6_fit"], d["tail"]["C8_fit"]
    sp = CubicSpline(R, E)

    def V(r):
        r = np.asarray(r, float)
        out = np.where(r <= 12.0, sp(np.clip(r, 5.0, 12.0)), -C6 / r ** 6 - C8 / r ** 8)
        return np.where(r < 5.0, np.nan, out)
    u, cnt = fcc_shells(14)

    def coh(dnn):
        r = u * dnn
        e = 0.5 * np.sum(cnt * V(r))
        rc = u[-1] * dnn                                     # tail beyond the shells (continuum, density 4/a^3)
        rho = np.sqrt(2) / dnn ** 3
        e += 0.5 * rho * 4 * np.pi * (-C6 / (3 * rc ** 3) - C8 / (5 * rc ** 5))
        return float(e)
    m = minimize_scalar(coh, bounds=(5.6, 8.0), method="bounded")
    dnn = m.x
    v = dnn ** 3 / np.sqrt(2)                               # volume per molecule (a0^3)
    return dict(d_nn_a0=dnn, eps_Eh=-m.fun, eps_eV=-m.fun * EH / EV, v_a0_3=v, v_m3=v * A0 ** 3,
                rho_kg_m3=2 * MP / (v * A0 ** 3), pair_well_Eh=float(-E.min()), C6=C6)


def R_cross(eps_J, v_m3, m_kg):
    """mountain argument: weight of a column of height R on the surface of a body of radius R, per molecule,
    m g R = eps with g = (4 pi / 3) G rho R  ->  R = sqrt(3 eps v / (4 pi G m^2))"""
    return float(np.sqrt(3 * eps_J * v_m3 / (4 * np.pi * G * m_kg ** 2)))


def main():
    mp = json.load(open(os.path.join(OUT, "ciclo_p1_map.json")))
    sl = mp["slopes_c4"]
    cp = mp["couplings"]
    sol = h2_solid()
    # crossover radii
    Rx_h2 = R_cross(sol["eps_eV"] * EV, sol["v_m3"], 2 * MP)
    # atoms crushed by gravity (energy per proton = 1 Rydberg, volume per atom (4 pi/3) a0^3): engine's own 1s
    Rx_atom = R_cross(0.5 * EH, 4 * np.pi / 3 * A0 ** 3, MP)
    # exam-only (NOT engine inputs): water ice, sublimation 0.49 eV / molecule, rho = 917 kg/m^3
    v_ice = 18.015 * 1.66053907e-27 / 917.0
    Rx_ice = R_cross(0.49 * EV, v_ice, 18.015 * 1.66053907e-27)
    R_round_obs = 2.0e5          # exam: Mimas (198 km) is round, Proteus (~210 km) is not: ~200 km for icy bodies
    fy_implied = (R_round_obs / Rx_ice) ** 2
    sigma_implied = fy_implied * 0.49 * EV / v_ice
    grav = dict(
        pair_pp_over_coulomb=G * MP ** 2 / E2, pair_ep_over_coulomb=G * MP * ME / E2,
        exponent_body="+2 (energy per particle G M m / R ~ R^2 against a constant cohesive energy)",
        h2_solid_engine=sol, R_cross_h2_engine_m=Rx_h2, R_cross_atoms_crushed_m=Rx_atom,
        R_cross_ice_exam_m=Rx_ice, R_round_ice_observed_m=R_round_obs, f_y_implied_ice=fy_implied,
        yield_stress_implied_Pa=sigma_implied,
        note="R_round = sqrt(f_y) R_x; f_y (yield energy / cohesive energy) needs plasticity (dislocations, grain "
             "boundaries): NOT in the engine (debt). The implied f_y and yield stress are exam numbers, not inputs.")
    # ledger rows: term, regime, dominant term there, weight at the regime's natural scale, exponent y, source
    a0 = A0
    rows = [
        dict(term="massa (riposo) mc^2", regime="L < lambda_C (Dirac, z=1)", dominant="moto a c: hbar c / L",
             weight_at="L/lambda_C (1 a L = lambda_C)", y=+1, kind="rilevante",
             source="misurato: fase a riposo per passo = theta*B esatta a 1e-12 a ogni livello (agg_zz_3d_A)"),
        dict(term="cinetico hbar^2 k^2/2m", regime="lambda_C < L (Schroedinger, z=2)", dominant="se stesso",
             weight_at="1", y=0, kind="marginale (punto fisso)",
             source=f"misurato: m_kin/m al livello = 1 entro {mp['exams']['P1a']['value']:.1e} su tutti i livelli liberi"),
        dict(term="relativistico -k^4/8m^3", regime="Schroedinger", dominant="cinetico",
             weight_at="(lambda_C/L)^2/4", y=-2, kind="irrilevante",
             source=f"misurato 1-D: pendenza {sl['1D_ladder']:.3f} (scala b=2), {sl['1D_auto']:.3f} (auto); 3-D libero: "
                    f"NON determinato dai dati ({sl['3D_A']:.2f}); catena T: portato esatto (1,0000)"),
        dict(term="Coulomb -Z alpha hbar c / r", regime="Schroedinger", dominant="cinetico",
             weight_at="2 L / a0 (=1 a L = a0/2)", y=+1, kind="rilevante",
             source="analitico (il campo di Galerkin e' -Z alpha/r ristretto, nessun parametro); l'atomo nasce dove si "
                    "incrociano: 1s a a0/16 a 0,0018 % dalla regola fine"),
        dict(term="Darwin (1/8m^2) lap V", regime="atomo", dominant="Coulomb",
             weight_at="(lambda_C/L)^2/8 = alpha^2/8 a L = a0", y=-2, kind="irrilevante",
             source=f"misurato: g_k {cp['gk']['values']} a B=256,512,1024, y = {cp['gk']['y']:.3f}"),
        dict(term="campo-moto (1/4) grad V . grad", regime="atomo", dominant="Coulomb",
             weight_at="(lambda_C/L)^2/4", y=-2, kind="irrilevante",
             source=f"misurato: g_x {cp['gx']['values']}, y = {cp['gx']['y']:.3f}"),
        dict(term="spin-orbita (1/4m^2) sigma.(grad V x p)", regime="atomo", dominant="Coulomb",
             weight_at="(lambda_C/L)^2/4", y=-2, kind="irrilevante",
             source=f"misurato: g_so {cp['gso']['values']}, y = {cp['gso']['y']:.3f}"),
        dict(term="contatto (cuspide degli s, cella h)", regime="atomo", dominant="Coulomb",
             weight_at="~ (Z alpha)^2 theta^2 ln theta (fine), ~ h/a0 sugli s (grosso)", y="~ -2 (log)",
             kind="irrilevante ma il piu' lento", source="NUCLEO_1S (legge d0 + c a^(2 gamma)); AGGLOMERATO_COMPLETO s.4.3"),
        dict(term="QED (Lamb, g-2)", regime="atomo", dominant="Coulomb", weight_at="alpha^3 ln(1/alpha)", y="0 (relativo)",
             kind="assente nel motore", source="debito #2"),
        dict(term="legame chimico (H2: D_e 4,75 eV)", regime="L ~ a0 .. 3 a0 (molecola)", dominant="se stesso",
             weight_at="1", y="0 per molecola", kind="nuovo grado di liberta' (gli atomi)",
             source="AGGLOMERATO_DUE_ATOMI s.12 (H2: R_e -0,007 %, D_e +0,7 meV)"),
        dict(term="dispersione C6/R^6 (H2-H2)", regime="fra molecole", dominant="repulsione di scambio a corto raggio",
             weight_at=f"buca {sol['pair_well_Eh'] * 1e3:.3f} mEh a 6,53 a0 (pair), coesione statica fcc "
                       f"{sol['eps_eV'] * 1e3:.1f} meV/molecola", y="-6 (coppia vs distanza); 0 per molecola in un corpo",
             kind="marginale (per molecola)", source="aggmol fase 3 (C6 10,86, buca 0,106 mEh) + somma fcc (questo studio)"),
        dict(term="legame a idrogeno", regime="acqua", dominant="-", weight_at="-", y="0 per molecola", kind="non ancora nel motore",
             source="serve l'acqua (MOLECOLE.md)"),
        dict(term="gravita' G m^2 / r", regime="tutti", dominant="Coulomb (atomo) / coesione (corpo)",
             weight_at=f"{G * MP ** 2 / E2:.2e} (p-p), {G * MP * ME / E2:.2e} (e-p) nell'atomo", y="0 fra due particelle; +2 "
             "per particella in un corpo neutro", kind="rilevante per corpi grandi",
             source=f"analitico; incrocio derivato: H2 solido del motore R_x = {Rx_h2 / 1e3:.0f} km, atomi schiacciati "
                    f"{Rx_atom / 1e3:.0f} km"),
    ]
    res = dict(ledger=rows, gravity=grav)
    json.dump(res, open(os.path.join(OUT, "ciclo_p1_ledger.json"), "w"), indent=1, default=float)
    print(json.dumps(sol, indent=1))
    print(f"R_x H2 (engine) = {Rx_h2 / 1e3:.0f} km; atoms crushed = {Rx_atom / 1e3:.0f} km; ice (exam) = {Rx_ice / 1e3:.0f} km; "
          f"f_y implied (200 km) = {fy_implied:.2e}, yield stress {sigma_implied / 1e6:.1f} MPa")


if __name__ == "__main__":
    main()
