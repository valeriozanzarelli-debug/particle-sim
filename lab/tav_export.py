"""tav_export.py - esportazione per la modalita' 'tavola periodica' del programma (stesso formato della cartella
dell'elio, /mnt/project-files/ripresa/dati/elio/): manifest.json, posizioni degli elettroni campionate da |psi|^2,
densita' e dati dei livelli per l'assorbimento/riemissione di un fotone.
Uso: python lab/tav_export.py li | be   (solo per elementi che passano: vedi lab/out/tav_results.json)"""
import os, sys, json, time, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import tav_ecg as T
from tav_analyze import load_svm, PRE, REF, NUC, ALPHA, AU_S, HC_EVNM

DEST = {"li": "/mnt/project-files/debiti/dati/litio", "be": "/mnt/project-files/debiti/dati/berillio"}
NOME = {"li": "litio", "be": "berillio"}
SIGMA = {"p3": (0, 0, 1), "s4": (0, 1, 0, 1), "p4": (0, 0, 1, 1), "s2": (0, 1)}


def radial_densities(sv, r, chunk=4000):
    """densita' di un elettrone (somma su tutti) esatta dalla funzione d'onda: rho(r, theta) = rho0(r) + rho2(r) P2(cos theta)
    (rho2 = 0 per gli stati S). Elettroni / bohr^3."""
    E, c = sv.state()
    K = len(c); N = sv.s.N
    I, J = np.triu_indices(K)
    w = c[I] * c[J] * np.where(I == J, 1.0, 2.0)
    rho0 = np.zeros_like(r); rho2 = np.zeros_like(r)
    for i0 in range(0, len(I), chunk):
        sl = slice(i0, i0 + chunk)
        A = sv.A[I[sl]]; B = sv.A[J[sl]]; u = sv.u[I[sl]]; v = sv.u[J[sl]]
        for cP, Pm in sv.s.perms:
            Bp = np.einsum("ji,njk,kl->nil", Pm, B, Pm)
            C = A + Bp; Ci = np.linalg.inv(C)
            S0 = (np.pi ** N / np.linalg.det(C)) ** 1.5
            Sig = 0.5 * Ci
            for i in range(N):
                s2 = Sig[:, i, i]
                g = (2 * np.pi * s2[:, None]) ** -1.5 * np.exp(-r[None] ** 2 / (2 * s2[:, None]))     # (n, nr)
                if sv.s.L == 0:
                    f0 = S0[:, None] * g; f2 = 0
                else:
                    vp = np.einsum("ji,nj->ni", Pm, v)
                    Su = np.einsum("nij,nj->ni", Sig, u); Sv = np.einsum("nij,nj->ni", Sig, vp)
                    uSv = np.einsum("ni,ni->n", u, Sv)
                    gu = Su[:, i]; gv = Sv[:, i]
                    a = (uSv - gu * gv / s2)[:, None]; b = (gu * gv / s2 ** 2)[:, None] * r[None] ** 2
                    f0 = S0[:, None] * g * (a + b / 3); f2 = S0[:, None] * g * (2 * b / 3)
                rho0 += cP * (w[sl] @ f0)
                if sv.s.L == 1:
                    rho2 += cP * (w[sl] @ f2)
    return rho0, rho2


def transition_density(svS, svP, r, chunk=4000):
    """<S| sum_i delta(r_i - r) |P0> = t(r) cos(theta) (radiale t, elettroni/bohr^3), pesi di spin incrociati"""
    E1, c1 = svS.state(); E2, c2 = svP.state()
    KS, KP = len(c1), len(c2); N = svS.s.N
    I, J = np.meshgrid(np.arange(KS), np.arange(KP), indexing="ij"); I = I.ravel(); J = J.ravel()
    w = c1[I] * c2[J]
    t = np.zeros_like(r)
    cw = T.cross_weights(svS.s.spin, svP.s.spin)
    for i0 in range(0, len(I), chunk):
        sl = slice(i0, i0 + chunk)
        A = svS.A[I[sl]]; B = svP.A[J[sl]]; v = svP.u[J[sl]]
        for cP, Pm in cw:
            Bp = np.einsum("ji,njk,kl->nil", Pm, B, Pm); vp = np.einsum("ji,nj->ni", Pm, v)
            C = A + Bp; Ci = np.linalg.inv(C)
            S0 = (np.pi ** N / np.linalg.det(C)) ** 1.5
            Sig = 0.5 * Ci
            Sv = np.einsum("nij,nj->ni", Sig, vp)
            for i in range(N):
                s2 = Sig[:, i, i]
                g = (2 * np.pi * s2[:, None]) ** -1.5 * np.exp(-r[None] ** 2 / (2 * s2[:, None]))
                t += cP * (w[sl] @ (S0[:, None] * g * (Sv[:, i] / s2)[:, None] * r[None]))
    return t


def export(el, n_samp=100000):
    t0 = time.time()
    res = json.load(open(os.path.join(T.OUT, "tav_results.json")))
    key = "Li" if el == "li" else "Be"
    R = res[key]
    if not R["passa_tavola"]:
        raise SystemExit(f"{key} non passa gli esami della tavola: niente esportazione (vedi documento)")
    g, i, p = ("li2S", "liP1S", "li2P") if el == "li" else ("be1S", "beP2S", "be1P")
    out = DEST[el]; os.makedirs(out, exist_ok=True)
    nm = NOME[el]
    svS = load_svm(g); svP = load_svm(p)
    sc = T.M0 / T.M_T_THETA_A                      # theta -> 0: lunghezze x m_T/m0, A x (m0/m_T)^2
    rng = np.random.default_rng(2026)
    files = {}
    samp_info = {}
    for tag, sv in (("fondamentale", svS), ("eccitato", svP)):
        E, c = sv.state()
        wv = T.Wave(sv.s, sv.A * sc ** 2, sv.u, c)
        sig = SIGMA[sv.s.spin]
        X, info = T.metropolis(wv, sig, n_samp, rng, nw=2000, burn=300, every=10, scale=1.5)
        fn = f"{nm}_{tag}.npy"
        np.save(os.path.join(out, fn), X)
        rr = np.linalg.norm(X, axis=2)
        buca = {}
        for a in range(X.shape[1]):
            for b in range(a + 1, X.shape[1]):
                same = sig[a] == sig[b]
                d = np.linalg.norm(X[:, a] - X[:, b], axis=1)
                dsh = np.linalg.norm(X[:, a] - X[rng.permutation(len(X)), b], axis=1)   # stesse distribuzioni, scorrelate
                buca.setdefault("stesso_spin" if same else "spin_opposti", []).append(
                    dict(coppia=[a, b], frazione_sotto_0p3=float(np.mean(d < 0.3)), scorrelata=float(np.mean(dsh < 0.3))))
        samp_info[tag] = dict(spin_etichette=["su" if s_ == 0 else "giu" for s_ in sig], **info,
                              r2_somma_campioni=float(np.mean(np.sum(rr ** 2, 1))),
                              r_medio_per_elettrone=np.mean(rr, 0).tolist(), n=int(len(X)), buca=buca)
        files[fn] = (f"float32 (n,{sv.s.N},3): n 'fotografie' delle {sv.s.N} posizioni degli elettroni insieme, estratte da "
                     f"|psi|^2 con Metropolis (stato {tag}). Etichette di spin per colonna: {samp_info[tag]['spin_etichette']}.")
        print(tag, samp_info[tag], f"{time.time() - t0:.0f}s", flush=True)
    r = np.linspace(0, 20, 401)
    rS0, _ = radial_densities(svS, r * sc)          # densita' alla massa m_T sui raggi che a m0 corrispondono a r
    rP0, rP2 = radial_densities(svP, r * sc)
    tr = transition_density(svS, svP, r * sc)
    svI = load_svm(i)
    rI0, _ = radial_densities(svI, r * sc)
    j3 = sc ** 3                                    # rho(m0)(r) = sc^3 rho(mT)(r sc)
    dens = dict(r=r, rho_fondamentale=rS0 * j3, rho_eccitato_0=rP0 * j3, rho_eccitato_2=rP2 * j3, transizione_t=tr * j3,
                rho_ione=rI0 * j3)
    # penetrazione (idea geometrica, sez. 'Idee'): quanta carica dell'elettrone esterno (rho_atomo - rho_ione) sta dentro
    # il raggio del nocciolo (dove 4 pi r^2 rho_ione ha il minimo dopo il picco interno o, per Li+, dove contiene 99 %)
    cumI = np.cumsum(4 * np.pi * r ** 2 * dens["rho_ione"]) * (r[1] - r[0])
    Nion = cumI[-1]
    rc = float(r[np.searchsorted(cumI, 0.99 * Nion)]) if el == "li" else float(r[np.searchsorted(cumI, (Nion - 1) * 0.99)])
    m = r <= rc
    pen = {k: float(np.trapezoid((4 * np.pi * r ** 2 * (dens[k] - dens["rho_ione"]))[m], r[m])) for k in ("rho_fondamentale", "rho_eccitato_0")}
    np.savez_compressed(os.path.join(out, f"{nm}_densita.npz"), **dens)
    integ = {k: float(np.trapezoid(4 * np.pi * r ** 2 * v, r)) for k, v in dens.items() if k in ("rho_fondamentale", "rho_eccitato_0")}
    dz = float(np.trapezoid(4 * np.pi / 3 * r ** 3 * dens["transizione_t"], r))      # <S|sum z|P0> dalla densita' di transizione
    files[f"{nm}_densita.npz"] = ("r (bohr), rho_fondamentale(r) (elettroni/bohr^3, somma su tutti, sferica), "
                                  "rho_eccitato_0/_2: stato P (M=0) rho(r,theta) = rho0 + rho2 P2(cos theta) (theta dall'asse z), "
                                  "rho_ione: densita' dello ione (Li+ / Be+); transizione_t: densita' di transizione <S|sum delta(r_i-r)|P0> = t(r) cos(theta). Durante l'assorbimento "
                                  "la nuvola e' rho_S + 2 a(t) t(r) cos(theta) cos(omega t): un'oscillazione di dipolo lungo z (esatta dalla funzione d'onda).")
    mu = NUC[3 if el == "li" else 4] / (NUC[3 if el == "li" else 4] + 1)
    dE_eV = R["line_eV"]; lam = HC_EVNM / dE_eV
    levels = dict(
        unita="energie in Eh (theta->0, massa infinita) e eV (massa ridotta); lunghezze in bohr",
        fondamentale=dict(nome="1s2 2s 2S" if el == "li" else "1s2 2s2 1S", E_Eh=R["E_theta0"]["ground"], g=2 if el == "li" else 1),
        eccitato=dict(nome="1s2 2p 2P" if el == "li" else "1s2 2s2p 1P", E_Eh=R["E_theta0"]["P"], g=6 if el == "li" else 3),
        ione=dict(nome="Li+ 1s2 1S" if el == "li" else "Be+ 1s2 2s 2S", E_Eh=R["E_theta0"]["ion"]),
        riga=dict(dE_eV=dE_eV, lambda_nm_vuoto=lam, f_lunghezza=R["f_L"], f_velocita=R["f_V"], dipolo_z_au=abs(R["z_len"]),
                  A_s=1 / (R["tau_P_ns"] * 1e-9), vita_ns=R["tau_P_ns"],
                  dipolo_da_densita_transizione=abs(dz),
                  regola="assorbe solo fotoni con energia dE (larghezza naturale hbar/tau); dopo ~vita_ns riemette lo stesso fotone in una direzione a caso (dipolo: intensita' ~ sin^2 rispetto a z)"),
        ionizzazione_eV=R["IP_eV"],
        esami=R["esami"])
    json.dump(levels, open(os.path.join(out, f"{nm}_livelli.json"), "w"), indent=1, ensure_ascii=False, default=float)
    files[f"{nm}_livelli.json"] = "livelli (fondamentale, eccitato P, ione), riga ottica (dE, lambda, forza d'oscillatore, vita), ionizzazione, esiti degli esami"
    man = dict(
        cosa=f"{nm.capitalize()} dalla regola (ECG con la massa imparata della regola, theta->0): stato fondamentale e primo stato P, per la tavola periodica",
        sorgente=dict(basi=f"lab/out/tav_{g}.npz, tav_{p}.npz, tav_{i}.npz (repo particle-sim, non tracciati)",
                      massa_m0=T.M0, massa_mT=T.M_T_THETA_A, scala_lunghezze=1 / sc,
                      risultati="lab/out/tav_results.json", documento="docs/TAVOLA_LITIO_BERILLIO.md"),
        unita="bohr (a0); nucleo nell'origine, fermo",
        file=files, campioni=samp_info,
        controlli=dict(integrale_rho_fondamentale=integ["rho_fondamentale"], integrale_rho_eccitato=integ["rho_eccitato_0"],
                       r2_esatto_fondamentale=R["r2"], densita_al_nucleo_esatta=R["rho0"],
                       dipolo_da_funzioni=abs(R["z_len"]) * (T.M_T_THETA_A / T.M0), dipolo_da_densita=abs(dz),
                       raggio_nocciolo_bohr=rc, penetrazione_dentro_nocciolo=dict(S=pen["rho_fondamentale"], P=pen["rho_eccitato_0"]),
                       buca_di_fermi=samp_info["fondamentale"].get("buca")),
        note=["Stati stazionari: per animarli si pescano fotografie diverse a ogni fotogramma (come per l'elio).",
              "Assorbimento: mescolare per un tempo breve le fotografie del fondamentale e dell'eccitato con il peso |a(t)|^2, "
              "o usare la densita' di transizione per l'oscillazione a frequenza omega = dE/hbar.",
              "Senza struttura fine (spin-orbita assente dalla regola, DEBT-4): una sola riga al centroide."],
        script="lab/tav_export.py (repo)")
    json.dump(man, open(os.path.join(out, "manifest.json"), "w"), indent=1, ensure_ascii=False, default=float)
    print(json.dumps(man["controlli"], indent=1), f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    export(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 100000)
