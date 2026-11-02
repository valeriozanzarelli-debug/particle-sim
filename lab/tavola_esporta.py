"""tavola_esporta.py - writes the periodic-table exports of the Laboratorio Atomo (docs/ATOMO.md §10.2) from the
studies' own wavefunctions. Format: assets/esportazioni/<nome>/manifest.json + campioni.f32 (n x n_el x 3 float32 LE, bohr).

    python3 lab/tavola_esporta.py elio            # ~2 min: He 1^1S ECG (lab/out/multi_e_ecg_he1S_Z2.npz), 16000 pairs
    python3 lab/tavola_esporta.py h2 <cartella>   # H2 at R = 1.4: rebuilds Debiti's ECG with lab/h2v_wf.py into <cartella>
                                                  # (~12 min, resumable: run again), then 16000 pairs
Nothing physical is chosen here: the samples are Metropolis draws from |Psi|^2 of the study's state.
"""
import os, sys, json, time
import numpy as np
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "lab"))
OUT = os.path.join(REPO, "lab", "out")
ESP = os.path.join(REPO, "assets", "esportazioni")
N = 16000


def scrivi(nome, campioni, man):
    d = os.path.join(ESP, nome)
    os.makedirs(d, exist_ok=True)
    c = np.ascontiguousarray(campioni.astype("<f4"))
    c.tofile(os.path.join(d, "campioni.f32"))
    man["campioni"] = {"file": "campioni.f32", "n": int(len(c)), "unita": "bohr", "ordine": "n x n_elettroni x 3, float32 little endian"}
    json.dump(man, open(os.path.join(d, "manifest.json"), "w"), indent=1, ensure_ascii=False)
    print(nome, c.shape, man.get("controlli_campioni"))


def elio():
    # same recipe as the helium glow export (MOLTI_ELETTRONI, phase-5 mass m0): ECG coefficients rebuilt from the basis
    cwd = os.getcwd(); os.chdir(os.path.join(REPO, "lab"))
    from multi_e_ecg import System, SVM
    os.chdir(cwd)
    m_A = json.load(open(os.path.join(OUT, "multi_e_kernel.json")))["R3"]["m_kin_theta_A_used"]
    m0 = json.load(open(os.path.join(OUT, "multi_e_v5.json")))["checks"]["R3v5"]["m0"]
    A = np.load(os.path.join(OUT, "multi_e_ecg_he1S_Z2.npz"))["A"]
    s = System(2, 1, 2.0, m_A); sv = SVM(s, np.random.default_rng(1), 2.0, 20.0, 30); sv.A = A
    E_A = sv.rebuild(); E, c = sv.state(); c = c * sv.nset
    E_m0 = E_A * m0 / m_A
    A = A * (m0 / m_A) ** 2          # lengths scale by m_A/m0 (exact)

    def psi(r1, r2):
        d11 = np.sum(r1 * r1, 1); d22 = np.sum(r2 * r2, 1); d12 = np.sum(r1 * r2, 1)
        ex = -(A[:, 0, 0][None] * d11[:, None] + A[:, 1, 1][None] * d22[:, None] + (A[:, 0, 1] + A[:, 1, 0])[None] * d12[:, None])
        return np.exp(ex) @ c
    rng = np.random.default_rng(2026)
    nw = N; r1 = rng.normal(0, 0.6, (nw, 3)); r2 = rng.normal(0, 0.6, (nw, 3)); p = psi(r1, r2) ** 2
    for it in range(400):
        n1 = r1 + rng.normal(0, 0.35, r1.shape); n2 = r2 + rng.normal(0, 0.35, r2.shape); pn = psi(n1, n2) ** 2
        a = rng.random(nw) < pn / np.maximum(p, 1e-300)
        r1[a] = n1[a]; r2[a] = n2[a]; p[a] = pn[a]
    S = np.stack([r1, r2], 1)          # 16000 independent walkers after 400 moves
    r = np.linalg.norm(S, axis=2); r12 = np.linalg.norm(S[:, 0] - S[:, 1], axis=1)
    scrivi("He", S, {
        "formato": "particle-sim/esportazione-v1", "nome": "elio", "simbolo": "He", "tipo": "atomo", "Z": 2, "carica": 0,
        "livello": "esatto",
        "metodo": "ECG (200 gaussiane correlate) sull'operatore della regola agglomerata (massa di fase 5), stato 1¹S",
        "fonte": "lab/out/multi_e_ecg_he1S_Z2.npz (commit f7d085f), campioni Metropolis di lab/tavola_esporta.py elio (seme 2026)",
        "nuclei": [{"Z": 2, "pos": [0.0, 0.0, 0.0]}], "n_elettroni": 2, "simmetria": "sferica",
        "energia_Eh": float(E_m0), "riferimento_Eh": -2.903724377,
        "esami": [
            {"id": "M1v5", "cosa": "energia dell'elio -2,9037213 contro -2,9037244 Eh", "PASS": True},
            {"id": "E1v5", "cosa": "prima ionizzazione 24,58814 eV contro 24,58739", "PASS": True},
            {"id": "Q1v5", "cosa": "banda della regola contro Dirac cinetico sull'elio (0,9997)", "PASS": True},
            {"id": "he_v E1-E5", "cosa": "video a bagliori: densita', <r>, <r12>, buca di Coulomb, mappa g (VIDEO_ELIO)", "PASS": True},
            {"id": "D6'", "cosa": "cuspide del misto DMC 0,447 contro 0,50 (non tocca l'ECG)", "PASS": False}],
        "controlli_campioni": {"r_medio": float(r.mean()), "r_medio_rif": 0.929472, "r12_medio": float(r12.mean()), "r12_medio_rif": 1.42207},
        "note": ["stato stazionario: nessun moto nel tempo; i bagliori sono coppie pescate a ogni fotogramma (moto = campionatore, non traiettoria)",
                 "singoletto: parte spaziale simmetrica", "correzione relativistica solo nell'energia, non nella funzione d'onda"]})


def h2(cartella):
    import h2v_wf as W
    os.makedirs(cartella, exist_ok=True)
    W.OUT = cartella; W.JS = os.path.join(cartella, "h2v_wf.json")
    W.res = json.load(open(W.JS)) if os.path.exists(W.JS) else {}
    base = os.path.join(cartella, "h2v_basis_R1.40.npz")
    if not os.path.exists(base):
        e = W.run(1.4, None, time.time(), 500)
        if e is None:
            print("checkpoint: run again"); return
        ref = json.load(open(os.path.join(OUT, "aggmol_h2.json")))["curve"]["1.40"]["E_tot"]
        print("E_tot", e + 1 / 1.4, "Debiti", ref, "diff", e + 1 / 1.4 - ref)
    os.environ["H2V_SCR"] = cartella
    import h2v_core as C, h2v_sim as Sm
    wf = C.WF(1.4, path=base)
    rng = np.random.default_rng(2026)
    tgt = lambda X: wf.psi(X) ** 2
    X = Sm.init_walkers(4000, 1.4, rng)
    X, p, acc = Sm.metropolis(X, tgt, 0.35, rng, 200)
    out = []
    for k in range(4):
        X, p, acc = Sm.metropolis(X, tgt, 0.35, rng, 25, p)
        out.append(X.copy())
    c = np.concatenate(out)
    scrivi("H2", c, {
        "formato": "particle-sim/esportazione-v1", "nome": "H2", "simbolo": "H2", "tipo": "molecola", "Z": 0, "carica": 0,
        "livello": "esatto",
        "metodo": "ECG a due centri (180 gaussiane correlate) sulla regola a molti elettroni del thread principale, stato fondamentale a R = 1,4 a0 (nuclei fermi)",
        "fonte": "lab/aggmol_h2.py (Debiti, aggmol fase 2) ricostruita con lab/h2v_wf.py (E identica, diff 2e-13 Eh); 16000 coppie Metropolis di lab/tavola_esporta.py h2 (seme 2026)",
        "nuclei": [{"Z": 1, "pos": [0.0, 0.0, 0.7]}, {"Z": 1, "pos": [0.0, 0.0, -0.7]}],
        "n_elettroni": 2, "simmetria": "assiale_z", "energia_Eh": -1.1744677, "riferimento_Eh": -1.174475714,
        "esami": [
            {"id": "H0", "cosa": "E(R = 1,4) -1,1744677 contro -1,1744757 Eh", "PASS": True},
            {"id": "H1", "cosa": "lunghezza di legame 1,40090 contro 1,401 a0", "PASS": True},
            {"id": "H2", "cosa": "energia di legame 4,7477 contro 4,747 eV", "PASS": True},
            {"id": "H3", "cosa": "quanto vibrazionale 4419 contro 4401 cm-1", "PASS": True},
            {"id": "H4", "cosa": "salto 0-1 4167 contro 4161 cm-1, D0 4,4772 eV", "PASS": True},
            {"id": "H5", "cosa": "Hartree-Fock sulla griglia a0/4 +0,94 mEh", "PASS": True},
            {"id": "G1", "cosa": "cancello della macchina ECG (soglia 1e-6 troppo stretta)", "PASS": False},
            {"id": "V3", "cosa": "palline contro densita' (Pearson 0,99990)", "PASS": True},
            {"id": "V4", "cosa": "elettroni da parti opposte 0,5749 (HF 0,5): la correlazione", "PASS": True}],
        "controlli_campioni": {"P_opposti": float(np.mean(c[:, 0, 2] * c[:, 1, 2] < 0)), "P_opposti_rif": 0.574871665125156,
                               "abs_z": float(np.mean(np.abs(c[:, :, 2]))), "abs_z_rif": 0.7977708301775607, "accettazione": float(acc)},
        "note": ["nuclei fermi a R = 1,4 a0 (Born-Oppenheimer; la vibrazione non e' in questa esportazione)",
                 "bagliori = coppie pescate da |Psi|^2 (il moto e' il campionatore, non una traiettoria)",
                 "singoletto: parte spaziale simmetrica"]})


if __name__ == "__main__":
    if sys.argv[1] == "elio":
        elio()
    elif sys.argv[1] == "h2":
        h2(sys.argv[2])
