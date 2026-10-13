"""multi_e_summary.py - giudizio contro la pre-registrazione (lab/out/multi_e_prereg.json): raccoglie i risultati delle corse
e scrive lab/out/multi_e_results.json (PASS / FAIL / NON FATTO / FERMO per gate ed esami, con i numeri)."""
import numpy as np
from multi_e_core import load_json, save_json, EH_EV, NUC_MASS

K = load_json("multi_e_kernel.json"); D = load_json("multi_e_disp.json"); E = load_json("multi_e_ecg.json")
M = load_json("multi_e_mf.json"); L = load_json("multi_e_lattice.json"); Q = load_json("multi_e_dmc.json")
R = {}


def put(name, verdict, **kw):
    R[name] = dict(verdict=verdict, **kw)


put("G0", "PASS" if K["G0"]["PASS"] else "FAIL", max_eigphase_diff=K["G0"]["max_eigphase_diff"],
    band_diff=K["G0"]["max_pauli_band_phase_diff"], learned_odd=K["G0"]["learned_stencil_odd_rel"])
g1 = E["G1"]
put("G1", "FAIL", parts={k: dict(PASS=v["PASS"], **{q: v[q] for q in ("err", "diff", "worst_rel") if q in v}) for k, v in g1.items() if isinstance(v, dict)},
    diagnosis_posthoc=dict(Hooke=E.get("G1_diag_posthoc", {}).get("Hooke"),
                           MC={k: v for k, v in E.get("G1_diag_posthoc", {}).get("MC", {}).items() if k != "recs"}))
put("G2", "RINVIATO (DV6)")
g3 = M.get("G3", {})
put("G3", "PASS" if all(g3.get(a, {}).get("PASS") for a in ("He", "H-", "Li", "Be", "Ne")) else
    ("PARZIALE: Ne non finito" if all(g3.get(a, {}).get("PASS") for a in ("He", "H-", "Li", "Be")) else "FAIL"),
    rel_err={a: g3.get(a, {}).get("rel_err") for a in ("He", "H-", "Li", "Be", "Ne")},
    random_starts={a: v["max_diff"] for a, v in M.get("starts", {}).items()})
c3 = M.get("carrier", {}).get("He3S", {})
put("G5", "FAIL (al livello dell'arrotondamento, vedi diagnosi)" if c3 and not c3.get("G5_monotone") else ("PASS" if c3 else "NON FATTO"),
    worst_step=c3.get("worst_step"), diagnosis="passi non monotoni solo per r < 4.2e-7 a0, ampiezza <= 1.4e-15: il rapporto e' "
    "c0 + c2 r^2 vicino al nucleo (le due cuspidi uguali cancellano il termine lineare), piatto sotto l'arrotondamento")
put("G4", "NON VALUTABILE (serve la base ECG, ferma per G1)")
put("R1", "PASS" if K["R1"]["PASS"] else "FAIL", fine_E50=K["R1"]["fine_energy_vs_beta_mc2"]["50"], pauli_rel=K["R1"]["rel_err"])
put("R2", "PASS" if K["R2"]["PASS"] else "FAIL", beta_min=K["R2"]["summary"]["B16"]["beta_min"], D=K["R2"]["D"],
    D_real=K["R2"]["D_real"], kurtosis=K["R2"]["kurtosis"], chosen=K["R2"]["chosen_walker_level"])
put("R3", "PASS" if K["R3"]["PASS"] else "FAIL", m0=K["R3"]["m0"], sigma=K["R3"]["sigma_seeds"])
det = L.get("det", {})
put("R4a", "PASS" if det.get("PASS_a") else "FAIL", imag_rel=det.get("rel_imag_vs_lobpcg"), real_rel=det.get("rel_realtime_vs_lobpcg"),
    diagnosis_posthoc={k: L.get("diag_R4a", {}).get(k) for k in ("imag_rel", "real_rel", "c2_imag", "c2_real")})
wk = {b: L.get(f"walk_beta{b}") for b in ("0.0002", "0.0001")}
put("R4b", "PASS" if all(w and w.get("PASS_b") for w in wk.values()) else "FAIL",
    **{f"beta{b}": dict(E=w["E_mixed_walkers"], sigma=w["sigma"], det=w["E_mixed_det"], z=w["z"]) for b, w in wk.items() if w})
put("relearned_level", "IL RISCHIO 4 SI E' AVVERATO: nessun livello Pauli passa il controllo a 3e-5; la macchina sceglie dirac4 B=2",
    chosen=D["relearn"]["A"]["chosen"], best_pauli_check=D["relearn"]["A"]["best_pauli_rejected_check"])
for ex in ("R5", "R6", "M1 He 1^1S", "M2 He 2^3S", "M3 H-", "E1", "E2", "D1", "D2", "D6"):
    put(ex, "FERMO (G1 FAIL: produzione ECG sospesa)" if ex[0] in "RME" or ex == "D1" or ex == "D2" else "NON FATTO")

# DMC dell'elio (riportato, non giudicato: D1 confronta con l'ECG, ferma)
he = {k: v for k, v in Q.items() if k.startswith("He_dmc")}
if he:
    taus = sorted(v["tau"] for v in he.values())
    pts = [(he[f"He_dmc_tau{t:g}"]["tau"], he[f"He_dmc_tau{t:g}"]["E"], he[f"He_dmc_tau{t:g}"]["sigma"]) for t in taus]
    rep = dict(points=pts)
    if len(pts) >= 3:
        t, e, s = map(np.array, zip(*pts))
        Wm = np.diag(1 / s ** 2)
        A = np.stack([np.ones_like(t), t], 1)
        cov = np.linalg.inv(A.T @ Wm @ A); c = cov @ A.T @ Wm @ e
        chi2 = float(np.sum((e - A @ c) ** 2 / s ** 2))
        m0, mA = K["R3"]["m0"], K["R3"]["m_kin_theta_A_used"]
        E0 = c[0] * m0 / mA
        sE0 = np.sqrt(cov[0, 0]) * m0 / mA
        sm = abs(c[0]) * K["R3"]["sigma_seeds"]
        mu = NUC_MASS[2] / (NUC_MASS[2] + 1)
        IE = ((-2.0 * m0) - E0) * mu * EH_EV
        sIE = np.hypot(sE0, sm) * mu * EH_EV
        rep.update(E_tau0_thetaA=float(c[0]), sigma_tau0=float(np.sqrt(cov[0, 0])), chi2_linear=chi2,
                   E_theta0=float(E0), sigma_theta0=float(np.hypot(sE0, sm)), target=-2.903724377,
                   diff_mEh=float((E0 + 2.903724377) * 1e3), IE_eV=float(IE), sigma_IE_eV=float(sIE), IE_target=24.587389)
    put("He_DMC_report_only", "RIPORTATO (non pre-registrato come esame: D1 e E1 vanno con l'ECG)", **rep)
save_json("multi_e_results.json", R)
for k, v in R.items():
    print(f"{k:22s} {v['verdict']}")
if "He_DMC_report_only" in R:
    print(R["He_DMC_report_only"])
