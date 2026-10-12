"""agg_zz_export.py - exports the learned hierarchy (levels, b, w, descriptor, check errors, stencils) to
assets/atomo/agglomerato_zz.json for the Rust lab (run the atom at the Bohr scale with the learned coarse rule).
Run after agg_zz_1d.py, agg_zz_3d.py A/B, agg_zz_atom.py, agg_zz_zoom.py."""
import json, os, pickle, datetime
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
ASSET = os.path.join(HERE, "..", "assets", "atomo", "agglomerato_zz.json")
A0 = 137.035999084


def cplx(a):
    a = np.asarray(a)
    return {"re": np.round(a.real, 15).tolist(), "im": np.round(a.imag, 15).tolist()}


def chain(name, sub):
    st = pickle.load(open(os.path.join(OUT, f"agg_zz_3d_{name}.pkl"), "rb"))
    th = A0 / 16 / sub
    levels = []
    for l in st["levels"]:
        levels.append(dict(
            B_fine_cells=l["B"], b=l["b"], w=l["w"], descriptor=l["desc"], components=l["C"],
            cell_lambdaC=l["B"] * th, cell_over_a0=l["B"] * th / A0, tick_hbar_over_mc2=l["B"] * th,
            observable_state_size_a0=1 / (l["k_obs_over_m"] * A0),
            check=dict(packet=l["check"]["packet"], energy_rel=l["check"]["energy_rel"], rest_rel=l["check"]["rest_rel"],
                       bz_rel90=l["check"].get("bz_rel90"), neg_gen=l["check"].get("neg_gen")),
            rest_phase_phi0=l["phi0"],
            offsets=l["offs"].tolist(),
            generator_G=cplx(l["K"]),
            spin_part_rel=l["disp"].get("spin_part_rel"),
            n_candidates_tried=len(l["candidates"])))
    return dict(theta_fine=th, fine_cell_lambdaC=th, levels=levels, stopped_no_candidate=st.get("stopped_no_candidate"))


d1 = json.load(open(os.path.join(OUT, "agg_zz_1d.json")))
atom = json.load(open(os.path.join(OUT, "agg_zz_atom.json"))) if os.path.exists(os.path.join(OUT, "agg_zz_atom.json")) else {}
zoom = json.load(open(os.path.join(OUT, "agg_zz_zoom.json"))) if os.path.exists(os.path.join(OUT, "agg_zz_zoom.json")) else {}
E = dict(
    what="Agglomerato della regola zig-zag: gerarchia di livelli appresi automaticamente (docs/AGGLOMERATO_ZIGZAG.md)",
    written=datetime.date.today().isoformat(),
    units="hbar = c = m_e = 1: lunghezze in lambda_C = 3.8616e-13 m, tempi in hbar/(m c^2) = 1.2881e-21 s, energie in m c^2",
    how_to_use=("Un livello con descrittore 'pauli2' agisce su phi (2 componenti = spinore di Pauli, la parte grande dello "
                "stato a energia positiva) su un reticolo cubico di passo h = cell_lambdaC. Un passo dura tau = h/c "
                "(= tick_hbar_over_mc2). Regola: in Fourier phi(k) -> e^{i phi0} exp(-i G(k)) phi(k), con "
                "G(k) = sum_j G_j e^{i k.j}, j = offsets (in celle), G_j = generator_G[a][b][j] (matrice 2x2, Hermitiana in "
                "totale). Con il potenziale: e^{-i V tau/2} [regola] e^{-i V tau/2} (V = media di cella di -alpha/r). "
                "Per gli stati legati: H = G/tau + V (Hermitiano). La parte di spin di G e' ~1e-10: si puo' usare "
                "(G_00 + G_11)/2 come stencil scalare (Schroedinger con correzione relativistica appresa)."),
    zoom_levels=("livello <-> zoom: l'osservabile fissa la cella. Stato 1s (taglia a0): cella <= a0/16 -> livello B=256 "
                 "della catena A; stati n=2 (taglia 2 a0): cella <= a0/8 -> B=512; dentro il nucleo (lambda_C) serve la "
                 "regola fine (zoom ibrido, catena B: fine theta=0.268 nella cella del protone, B=32 fuori)."),
    chain_A_3d=chain("A", 256),
    chain_B_3d=chain("B", 32),
    chain_1d=dict(theta_fine=d1["theta"], ladder=[{k: r[k] for k in ("B", "b", "w", "taps", "exam_err", "m_kin_over_m", "c4_over_1_8m3")} for r in d1["ladder"]],
                  auto=[{k: r[k] for k in ("B", "b", "w", "taps", "exam_err", "m_kin_over_m", "c4_over_1_8m3")} for r in d1["auto"]]),
    hydrogen=atom,
    zoom=zoom,
)
json.dump(E, open(ASSET, "w"), indent=1, default=float)
print("written", os.path.abspath(ASSET), os.path.getsize(ASSET) // 1024, "kB")
