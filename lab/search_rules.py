"""search_rules.py - ricerca automatica di regole locali luce + "struttura" che formano oggetti
stabili (docs/RICERCA_REGOLE.md, parte B).

Spazio delle regole (2+1 D, c = 1), un candidato = vettore theta di 8 numeri:
  campi: A_mu (luce, U(1)), n = (n1, n2, n3), |n| = 1, campo di struttura (vettore unitario, modello sigma O(3));
  la luce ruota (n1, n2) = psi (carica g), n3 e' neutro.
  energia statica (A_0 = 0):
    E = sum h^2 [ 1/2 |D n|^2                            (gradienti covarianti, D = d - i g A su psi)
                 + V(n)                                   (potenziale, sotto)
                 + kappa^2/2 |D_x n x D_y n|^2            (termine di Skyrme, "baby-Skyrme")
                 + 1/2 (1 + zeta (1-n3)) B^2              (Maxwell con costante dielettrica che dipende
                                                           dalla struttura: idea "dielettrico di colore")
                 + a4 B^4 + a6 B^6 ]                      (elettrodinamica non lineare
                                                           L_EM = -F - 4 a4 F^2 - 8 a6 F^3, F=(B^2-E^2)/2)
    V = mu1^2 (1-n3) + mu2^2 (1-n3)^2 + nu n3^2  - V_vuoto
  theta = (g, zeta, mu1, mu2, nu, kappa, a4, a6) nella scatola BOX.
Contiene come casi particolari: Maxwell (g=0, phi banale), Maxwell non lineare (Born-Infeld/
Euler-Heisenberg troncati), baby-Skyrme (g=0), modello sigma O(3) con gauge (Schroers 1995),
baby-Skyrme con gauge (Gladikowski-Piette-Schroers 1996), Higgs abeliano (nu>0: vuoto all'equatore,
vortici di Nielsen-Olesen nel limite di London).
Nota: nel codice n = phi/|phi| (phi libero, solo la direzione conta): il modulo non e' un grado di liberta'.

Comandi:
  python lab/search_rules.py prereg           scrive il pre-registro (fitness e vincoli FISSATI) + hash
  timeout 590 python lab/search_rules.py cma  (riprende dal checkpoint, si ferma da solo dopo ~400 s)
  timeout 590 python lab/search_rules.py rand
  python lab/search_rules.py known            valuta le teorie note dentro lo spazio
  timeout 590 python lab/search_rules.py dyn FILE.json ETICHETTA [pert]  dinamica vera, 100 tempi di attraversamento
"""
import hashlib
import json
import os
import sys
import time

import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)

NAMES = ["g", "zeta", "mu1", "mu2", "nu", "kappa", "a4", "a6"]
BOX = np.array([[0.0, 3.0], [-0.45, 3.0], [0.0, 1.5], [0.0, 1.5], [-1.0, 1.0], [0.0, 1.5], [-0.3, 0.3], [0.0, 0.3]])
N_GRID, H = 48, 0.3
L_BOX = N_GRID * H

PREREG = {
    "data": "2026-10-09, scritto PRIMA di qualsiasi valutazione della ricerca",
    "spazio": "theta = " + ", ".join(NAMES) + "; scatola " + json.dumps(BOX.tolist()),
    "griglia": f"{N_GRID}x{N_GRID}, passo h={H}, bordo: phi fissato al vuoto, A libero",
    "vincoli_duri": {
        "C0": "buona posizione: energia EM convessa in B (se a4<0: 144 a4^2 < 120 a6)",
        "C1": "luce libera invariata: nel vuoto la parte quadratica in A e' esattamente Maxwell: "
              "massa del fotone g*rho_vac < 1e-6 (rho_vac = |psi| nel vuoto), vuoto non degenere con rho>0 se g>0",
        "C2a": "esiste una configurazione ferma localizzata: rilassamento (flusso di Newton arrestato, 1500 passi) stazionario "
               "(variazione relativa di E negli ultimi 300 passi < 5e-3 e |forza|max < 0.1) e "
               ">= 90% dell'energia entro il disco di raggio L/4 = 3.6 dal baricentro (Q=2,3: 3.6 sqrt(Q), max L/2-1)",
        "C2b": "robusta: perturbazione casuale 5% + nuovo rilassamento => stessa energia entro 1% e stessa carica "
               "(provato solo se C2a e C3 passano; altrimenti conta come non superato)",
        "C3": "carica topologica intera conservata: |Q| finale = |Q| iniziale (=1) e oggetto risolto (raggio rms "
              "della densita' topologica > 2h)  [vuoto al polo: grado di n = phi/|phi| (Berg-Luscher); "
              "vuoto all'equatore: avvolgimento al bordo]",
    },
    "C2_finalisti": "dinamica hamiltoniana vera (RK4), oggetto fermo, 100 tempi di attraversamento (2 r_E/c), "
                    ">= 90% dell'energia iniziale nel disco; ripetuta con perturbazione 5%",
    "punteggi_morbidi": {
        "a_Derrick": "S_a = tanh(E''(1)/E) * exp(-|E'(1)|/(0.1 E))  (scala x->x/s, termini separati)",
        "b_spin": "S_b = 1 se la carica e' in pi2(S^2) (termine di Hopf theta ammesso => spin qualsiasi, anche 1/2; "
                  "Wilczek-Zee 1983), 0 altrimenti",
        "c_spettro": "S_c = (1/2) #{Q in 2,3 : stazionario, localizzato, |Q| conservato e risolto, E_Q < Q E_1 (1+0.005)}",
        "d_luce": "S_d = min(1, 4 f_EM), f_EM = energia elettromagnetica / energia totale (oggetto Q=1)",
    },
    "pesi": {"a": 0.2, "b": 0.2, "c": 0.3, "d": 0.3},
    "fitness": "F = [tutti i vincoli duri] * (0.2 S_a + 0.2 S_b + 0.3 S_c + 0.3 S_d)",
    "rango_per_la_selezione": "n_vincoli_superati (0..5) + F  (solo per ordinare; il risultato riportato e' F)",
    "budget": "CMA-ES: 10 generazioni x 8 = 80 valutazioni, seme 2026; ricerca casuale uniforme: 80, seme 2027",
    "test_di_predizione": "solo DOPO la ricerca, sui migliori, senza ritoccare: una massa fissa la scala; "
                          "rapporti E_2/E_1, E_3/E_1, r*m_q confrontati con mu/e, tau/e, p/e, Delta/N",
}
WEIGHTS = PREREG["pesi"]


def decode(u):
    """u in [0,1]^8 -> dict di parametri."""
    t = BOX[:, 0] + np.clip(u, 0, 1) * (BOX[:, 1] - BOX[:, 0])
    return dict(g=t[0], zeta=t[1], mu1=t[2], mu2=t[3], nu=t[4], kappa=t[5], a4=t[6], a6=t[7])


def Vn3(p3, P):
    return P["mu1"] ** 2 * (1 - p3) + P["mu2"] ** 2 * (1 - p3) ** 2 + P["nu"] * p3 ** 2


def vacuum(P):
    """vuoto sulla sfera: minimo di V(n3), n3 in [-1, 1]; rho = |psi| = sqrt(1 - n3^2)."""
    z = np.linspace(-1, 1, 200001)
    v = Vn3(z, P)
    k = int(np.argmin(v))
    Vmin, p3 = float(v[k]), float(z[k])
    others = z[np.abs(v - Vmin) < 1e-10]
    degenerate_rho = bool(np.any(np.abs(others) < 1 - 1e-4))
    rho = float(np.sqrt(max(0.0, 1 - p3 ** 2)))
    if abs(abs(p3) - 1) < 1e-9:
        rho, p3 = 0.0, float(np.sign(p3))
    e = 1e-4
    if rho == 0.0:
        def Vr(r):
            return Vn3(p3 * np.sqrt(1 - r ** 2), P)
        m2q = (Vr(e) - 2 * Vr(0.0) + Vr(-e)) / e ** 2
        m23 = None
    else:
        th0 = np.arccos(p3)

        def Vt(t):
            return Vn3(np.cos(t), P)
        m23 = (Vt(th0 + e) - 2 * Vt(th0) + Vt(th0 - e)) / e ** 2
        m2q = 0.0  # modo di Goldstone (mangiato dal fotone se g > 0)
    return dict(Vmin=Vmin, rho=rho, p3=p3, degenerate_rho=degenerate_rho, m2_charged=float(m2q),
                m2_neutral=None if m23 is None else float(m23))


def _torch():
    import torch
    torch.set_num_threads(1)
    torch.set_default_dtype(torch.float64)
    return torch


def params_tensor(Plist, torch):
    return {k: torch.tensor([P[k] for P in Plist]).reshape(-1, 1, 1) for k in Plist[0]}


def energy_terms(phi, A, Pt, vmin, torch, h=H):
    """phi (b,3,N,N), A (b,2,N,N), Pt: parametri (b,1,1), vmin (b,1,1). Ritorna somme per b e densita'."""
    g = Pt["g"]
    phi = phi / torch.sqrt((phi ** 2).sum(1, keepdim=True) + 1e-30)  # solo la direzione n conta
    th = g * h * A[:, 0, :-1, :]
    c, s = torch.cos(th), torch.sin(th)
    n1, n2 = phi[:, 0, 1:, :], phi[:, 1, 1:, :]
    Dx = torch.stack([n1 * c + n2 * s - phi[:, 0, :-1, :], -n1 * s + n2 * c - phi[:, 1, :-1, :],
                      phi[:, 2, 1:, :] - phi[:, 2, :-1, :]], 1) / h
    th = g * h * A[:, 1, :, :-1]
    c, s = torch.cos(th), torch.sin(th)
    n1, n2 = phi[:, 0, :, 1:], phi[:, 1, :, 1:]
    Dy = torch.stack([n1 * c + n2 * s - phi[:, 0, :, :-1], -n1 * s + n2 * c - phi[:, 1, :, :-1],
                      phi[:, 2, :, 1:] - phi[:, 2, :, :-1]], 1) / h
    eS2x = 0.5 * (Dx ** 2).sum(1)
    eS2y = 0.5 * (Dy ** 2).sum(1)
    # termine di Skyrme sulla plaquette: media dei 4 prodotti (lato basso/alto in x) x (lato sinistro/destro in y)
    # (la sola combinazione "avanti-avanti" ha un modo spurio a scacchiera che srotola il solitone)
    eSk = 0.0
    for xs in (slice(0, -1), slice(1, None)):
        for ys in (slice(0, -1), slice(1, None)):
            cr = torch.cross(Dx[:, :, :, ys], Dy[:, :, xs, :], dim=1)
            eSk = eSk + 0.125 * Pt["kappa"] ** 2 * (cr ** 2).sum(1)
    p3 = phi[:, 2]
    eV = (Pt["mu1"] ** 2 * (1 - p3) + Pt["mu2"] ** 2 * (1 - p3) ** 2
          + Pt["nu"] * p3 ** 2 - vmin)
    B = (A[:, 0, :-1, :-1] + A[:, 1, 1:, :-1] - A[:, 0, :-1, 1:] - A[:, 1, :-1, :-1]) / h
    p3p = 0.25 * (p3[:, :-1, :-1] + p3[:, 1:, :-1] + p3[:, :-1, 1:] + p3[:, 1:, 1:])
    eM = 0.5 * (1 + Pt["zeta"] * (1 - p3p)) * B ** 2
    e4 = Pt["a4"] * B ** 4
    e6 = Pt["a6"] * B ** 6
    h2 = h * h

    def sm(x):
        return h2 * x.sum((1, 2))
    T = dict(sigma=sm(eS2x) + sm(eS2y), V=sm(eV), Sk=sm(eSk), M=sm(eM), B4=sm(e4), B6=sm(e6))
    return T, dict(eS2x=eS2x, eS2y=eS2y, eSk=eSk, eV=eV, eM=eM + e4 + e6, Dx=Dx, Dy=Dy, B=B)


def density_map(d, i):
    e = d["eV"][i].detach().numpy().copy()
    e[:-1, :] += d["eS2x"][i].detach().numpy()
    e[:, :-1] += d["eS2y"][i].detach().numpy()
    e[:-1, :-1] += d["eSk"][i].detach().numpy() + d["eM"][i].detach().numpy()
    return e


def topo_charge(phi):
    """grado di n = phi/|phi| (Berg-Luscher, esatto sul reticolo) e densita'."""
    n = phi / np.maximum(np.linalg.norm(phi, axis=0), 1e-12)
    a, b, c, d = n[:, :-1, :-1], n[:, 1:, :-1], n[:, 1:, 1:], n[:, :-1, 1:]

    def omega(x, y, z):
        num = np.einsum("i...,i...->...", x, np.cross(y, z, axis=0))
        den = 1 + (x * y).sum(0) + (y * z).sum(0) + (z * x).sum(0)
        return 2 * np.arctan2(num, den)
    q = (omega(a, b, c) + omega(a, c, d)) / (4 * np.pi)
    return q.sum(), q


def winding_boundary(phi):
    ps = phi[0] + 1j * phi[1]
    ring = np.concatenate([ps[0, :], ps[1:, -1], ps[-1, -2::-1], ps[-2:0:-1, 0]])
    ring = np.where(np.abs(ring) < 1e-12, 1e-12, ring)
    dth = np.angle(ring[np.r_[1:len(ring), 0]] / ring)
    return dth.sum() / (2 * np.pi)


def grid():
    x = (np.arange(N_GRID) - (N_GRID - 1) / 2) * H
    return np.meshgrid(x, x, indexing="ij")


def init_config(Q, vac, r0=1.5):
    X, Y = grid()
    R = np.hypot(X, Y)
    TH = np.arctan2(Y, X)
    if vac["rho"] == 0.0:
        mod = 1.0
        sgn = 1.0 if vac["p3"] >= 0 else -1.0
        f = np.pi * np.exp(-R / (r0 * np.sqrt(Q)))
        edge = np.minimum(L_BOX / 2 - np.abs(X), L_BOX / 2 - np.abs(Y))
        f = f * np.clip(edge / 1.0, 0, 1)
        phi = mod * np.stack([np.sin(f) * np.cos(Q * TH), np.sin(f) * np.sin(Q * TH), sgn * np.cos(f)])
    else:
        prof = np.tanh(R / (r0 * np.sqrt(Q)))
        phi = np.stack([vac["rho"] * prof * np.cos(Q * TH), vac["rho"] * prof * np.sin(Q * TH),
                        vac["p3"] * np.ones_like(R)])
    return phi, np.zeros((2, N_GRID, N_GRID))


def flow(phi0, A0, Plist, vmins, torch, steps=1500, lbfgs_iter=0):
    """flusso di Newton arrestato, in parallelo su un lotto di configurazioni (ognuna con i suoi parametri).
    x'' = -grad E / h^2 ; se E sale la velocita' si azzera; se esplode si torna indietro e si dimezza dt."""
    b, N = len(Plist), N_GRID
    Pt = params_tensor(Plist, torch)
    vmin = torch.tensor(np.array(vmins, float)).reshape(-1, 1, 1)
    bmask = np.zeros((N, N), bool)
    bmask[0, :] = bmask[-1, :] = bmask[:, 0] = bmask[:, -1] = True
    fmask = torch.tensor(~bmask, dtype=torch.float64)[None, None]
    phi = torch.tensor(np.array(phi0))
    A = torch.tensor(np.array(A0))
    vphi, vA = torch.zeros_like(phi), torch.zeros_like(A)
    w2 = np.array([(8 / H ** 2) * (2 + 4 * P["kappa"] ** 2 + 6 * abs(P["a4"]) + 15 * P["a6"] + 2 * abs(P["zeta"]))
                   + 4 * (P["mu1"] ** 2 + P["mu2"] ** 2 + abs(P["nu"])) for P in Plist])
    dt = torch.tensor(0.4 / np.sqrt(w2)).reshape(-1, 1, 1, 1)
    nblow = np.zeros(b, int)

    def EG(phi, A):
        phi = phi.detach().requires_grad_(True)
        A = A.detach().requires_grad_(True)
        T, _ = energy_terms(phi, A, Pt, vmin, torch)
        Eb = sum(T.values())
        gphi, gA = torch.autograd.grad(Eb.sum(), [phi, A])
        return Eb.detach(), gphi * fmask / H ** 2, gA / H ** 2
    phi = phi / torch.sqrt((phi ** 2).sum(1, keepdim=True))
    E, gphi, gA = EG(phi, A)
    E80 = E.clone()
    for it in range(steps):
        if it == int(0.8 * steps):
            E80 = E.clone()
        vphi = vphi - dt * gphi
        vA = vA - dt * gA
        phi_new = phi + dt * vphi
        phi_new = phi_new / torch.sqrt((phi_new ** 2).sum(1, keepdim=True))
        vphi = vphi - (vphi * phi_new).sum(1, keepdim=True) * phi_new  # velocita' tangente alla sfera
        A_new = A + dt * vA
        E_new, g1, g2 = EG(phi_new, A_new)
        bad = ~torch.isfinite(E_new) | (E_new > E + 0.01 * E.abs() + 1e-9)
        up = E_new > E
        if bad.any():
            bb = bad.reshape(-1, 1, 1, 1)
            phi_new = torch.where(bb, phi, phi_new)
            A_new = torch.where(bb, A, A_new)
            dt = torch.where(bb, dt * 0.5, dt)
            nblow += bad.numpy()
            E_new, g1, g2 = EG(phi_new, A_new)
        stop = (up | bad).reshape(-1, 1, 1, 1)
        vphi = torch.where(stop, torch.zeros_like(vphi), vphi)
        vA = torch.where(stop, torch.zeros_like(vA), vA)
        phi, A, E, gphi, gA = phi_new, A_new, E_new, g1, g2
    if lbfgs_iter > 0:  # rifinitura: L-BFGS sulla somma delle energie (problemi indipendenti => Hessiano a blocchi)
        xphi = phi.clone().requires_grad_(True)
        xA = A.clone().requires_grad_(True)
        phib = phi.clone()
        opt = torch.optim.LBFGS([xphi, xA], lr=1, max_iter=lbfgs_iter, tolerance_grad=1e-12,
                                tolerance_change=1e-15, history_size=50, line_search_fn="strong_wolfe")
        fm = fmask.bool()

        def closure():
            opt.zero_grad()
            T, _ = energy_terms(torch.where(fm, xphi, phib), xA, Pt, vmin, torch)
            Et = sum(T.values()).sum()
            Et.backward()
            return Et
        try:
            opt.step(closure)
            phi_l = torch.where(fm, xphi, phib).detach()
            phi_l = phi_l / torch.sqrt((phi_l ** 2).sum(1, keepdim=True))
            E_l, g1, g2 = EG(phi_l, xA.detach())
            okl = torch.isfinite(E_l) & (E_l <= E + 1e-9)
            kk = okl.reshape(-1, 1, 1, 1)
            phi = torch.where(kk, phi_l, phi)
            A = torch.where(kk, xA.detach(), A)
            gphi = torch.where(kk, g1, gphi)
            gA = torch.where(kk, g2, gA)
            E = torch.where(okl, E_l, E)
        except Exception:
            pass
    # forza tangente (la componente radiale di grad E e' nulla per costruzione)
    fmax = torch.maximum(gphi.abs().amax((1, 2, 3)), gA.abs().amax((1, 2, 3)))
    with torch.no_grad():
        T, d = energy_terms(phi, A, Pt, vmin, torch)
    drift = ((E80 - E).abs() / E.abs().clamp_min(1e-12)).numpy()
    out = []
    for i in range(b):
        out.append(dict(drift=float(drift[i]), phi=phi[i].numpy(), A=A[i].numpy(), T={k: float(v[i]) for k, v in T.items()},
                        E=float(sum(v[i] for v in T.values())), fmax=float(fmax[i]), dens=density_map(d, i),
                        B=d["B"][i].numpy(), nblow=int(nblow[i])))
    return out


def analyse(res, Q0, vac):
    X, Y = grid()
    E = res["E"]
    dens = res["dens"]
    if not np.isfinite(E) or E <= 0:
        return dict(ok=False, E=float(E))
    w = np.maximum(dens, 0)
    xc, yc = (w * X).sum() / w.sum(), (w * Y).sum() / w.sum()
    Rr = np.hypot(X - xc, Y - yc)
    Rloc = min(L_BOX / 4 * np.sqrt(Q0), L_BOX / 2 - 1.0)
    frac = float(dens[Rr < Rloc].sum() / dens.sum())
    rE = float(np.sqrt((w * Rr ** 2).sum() / w.sum()))
    if vac["rho"] == 0.0:
        Q, q = topo_charge(res["phi"])
        Xq, Yq = X[:-1, :-1] + H / 2, Y[:-1, :-1] + H / 2
        rq = float(np.sqrt(abs((q * ((Xq - xc) ** 2 + (Yq - yc) ** 2)).sum() / q.sum()))) if abs(q.sum()) > 0.5 else 0.0
    else:
        Q = winding_boundary(res["phi"])
        rq = rE
    T = res["T"]
    Ep = 2 * T["V"] - 2 * (T["Sk"] + T["M"]) - 6 * T["B4"] - 10 * T["B6"]
    Epp = 2 * T["V"] + 6 * (T["Sk"] + T["M"]) + 42 * T["B4"] + 110 * T["B6"]
    fEM = (T["M"] + T["B4"] + T["B6"]) / E
    return dict(ok=True, E=E, Q=float(Q), Qint=int(round(abs(Q))), rq=rq, rE=rE, frac_in=frac, fmax=res["fmax"],
                Ep=Ep, Epp=Epp, fEM=fEM, T=T, flux=float(res["B"].sum() * H * H), centroid=[float(xc), float(yc)],
                nblow=res["nblow"], Bmax=float(np.abs(res["B"]).max()))


def prescreen(u):
    P = decode(np.asarray(u, float))
    out = dict(u=[float(v) for v in u], P={k: float(v) for k, v in P.items()})
    passed = {"C0": bool(not (P["a4"] < 0 and 144 * P["a4"] ** 2 >= 120 * P["a6"]))}
    vac = vacuum(P)
    out["vac"] = vac
    out["m_photon"] = float(P["g"] * vac["rho"])
    passed["C1"] = bool(out["m_photon"] < 1e-6 and not (P["g"] > 1e-6 and vac["degenerate_rho"]))
    out["topology"] = "pi2(S2)" if vac["rho"] == 0.0 else "pi1(S1)"
    Bs = []
    if abs(P["a4"]) > 0 or P["a6"] > 0:
        bb = np.logspace(-4, 3, 4000)
        rel = (np.abs(P["a4"]) * bb ** 2 + P["a6"] * bb ** 4) / 0.5
        Bs.append(float(bb[np.argmax(rel > 0.01)]))
    if P["g"] > 0 and vac["m2_charged"] > 0:
        Bs.append(float(vac["m2_charged"] / P["g"]))
    out["E_star"] = min(Bs) if Bs else None
    out["passed"] = passed
    return out, P, vac


def converged(a):
    return bool(a.get("ok") and a["drift"] < 5e-3 and a["fmax"] < 0.1)


def evaluate_population(Ulist, seed=0, keep_fields=False, steps=1500):
    """fase 1: Q=1 per tutti i candidati che passano C0, C1.  fase 2 (solo per chi passa C2a e C3):
    Q=1 perturbato del 5% (C2b) e Q=2, Q=3 (spettro)."""
    torch = _torch()
    t0 = time.time()
    outs, jobs = [], []
    for i, u in enumerate(Ulist):
        o, P, vac = prescreen(u)
        outs.append(o)
        if o["passed"]["C0"] and o["passed"]["C1"]:
            phi, A = init_config(1, vac)
            jobs.append((i, 1, P, vac, phi, A))
    res = flow([j[4] for j in jobs], [j[5] for j in jobs], [j[2] for j in jobs], [j[3]["Vmin"] for j in jobs],
               torch, steps) if jobs else []
    n1 = len(jobs)
    for j, r in zip(jobs, res):
        i, Q, P, vac = j[:4]
        a = analyse(r, Q, vac)
        a["drift"] = r["drift"]
        outs[i]["raw1"] = (a, r)
        o = outs[i]
        o["passed"]["C2a"] = bool(converged(a) and a["frac_in"] >= 0.9)
        o["passed"]["C3"] = bool(a.get("ok") and a["Qint"] == 1 and abs(abs(a["Q"]) - 1) < 0.05 and a["rq"] > 2 * H)
    # fase 2
    rng = np.random.default_rng(1000 + seed)
    jobs2 = []
    for i, o in enumerate(outs):
        if "raw1" in o and o["passed"]["C2a"] and o["passed"]["C3"]:
            a1, r1 = o["raw1"]
            ph = r1["phi"] * (1 + 0.05 * rng.standard_normal(r1["phi"].shape))
            ph[:, 0, :], ph[:, -1, :], ph[:, :, 0], ph[:, :, -1] = (r1["phi"][:, 0, :], r1["phi"][:, -1, :],
                                                                    r1["phi"][:, :, 0], r1["phi"][:, :, -1])
            Ap = r1["A"] + 0.05 * (np.abs(r1["A"]).max() + 1e-3) * rng.standard_normal(r1["A"].shape)
            P = decode(np.asarray(o["u"]))
            jobs2.append((i, "pert", P, o["vac"], ph, Ap))
            for Q in (2, 3):
                phi, A = init_config(Q, o["vac"])
                jobs2.append((i, Q, P, o["vac"], phi, A))
    res2 = flow([j[4] for j in jobs2], [j[5] for j in jobs2], [j[2] for j in jobs2],
                [j[3]["Vmin"] for j in jobs2], torch, steps) if jobs2 else []
    for j, r in zip(jobs2, res2):
        i, Q, P, vac = j[:4]
        a = analyse(r, 1 if Q == "pert" else Q, vac)
        a["drift"] = r["drift"]
        outs[i].setdefault("raw2", {})[Q] = (a, r)
    for o in outs:
        passed = o["passed"]
        if "raw1" not in o:
            passed.update(C2a=False, C2b=False, C3=False)
            o.update(n_pass=int(sum(passed.values())), F=0.0, rank=float(sum(passed.values())), S=None)
            continue
        a1, r1 = o["raw1"]
        raw2 = o.get("raw2", {})
        ap = raw2.get("pert", ({}, None))[0]
        passed["C2b"] = bool(ap.get("ok") and abs(ap["E"] - a1["E"]) / a1["E"] < 0.01 and ap["Qint"] == a1["Qint"])
        hard = all(passed.values())
        if a1.get("ok"):
            E = a1["E"]
            S = dict(a=float(np.tanh(max(a1["Epp"], 0) / E) * np.exp(-abs(a1["Ep"]) / (0.1 * E))),
                     b=1.0 if o["vac"]["rho"] == 0.0 else 0.0, d=float(min(1.0, 4 * max(a1["fEM"], 0))))
        else:
            S = dict(a=0.0, b=0.0, d=0.0)
        spec, nb = {}, 0
        for Q in (2, 3):
            if Q not in raw2:
                continue
            aQ, rQ = raw2[Q]
            if not aQ.get("ok"):
                spec[str(Q)] = dict(ok=False)
                continue
            okQ = bool(converged(aQ) and aQ["frac_in"] >= 0.9 and aQ["Qint"] == Q and aQ["rq"] > 2 * H
                       and aQ["E"] < Q * a1["E"] * 1.005)
            spec[str(Q)] = dict(E=aQ["E"], ratio=aQ["E"] / a1["E"], Q=aQ["Q"], frac_in=aQ["frac_in"],
                                fmax=aQ["fmax"], drift=aQ["drift"], rE=aQ["rE"], fEM=aQ["fEM"], ok=okQ)
            if keep_fields:
                spec[str(Q)]["dens"] = rQ["dens"].tolist()
            nb += okQ
        S["c"] = nb / 2 if hard else 0.0
        F = sum(WEIGHTS[k] * S[k] for k in "abcd") if hard else 0.0
        o["Q1"] = {k: v for k, v in a1.items()}
        o["Q1_pert"] = {k: ap.get(k) for k in ("E", "Q", "fmax", "drift")} if ap else None
        o["spectrum"] = spec
        if keep_fields:
            o["fields"] = dict(phi=r1["phi"].tolist(), A=r1["A"].tolist(), dens=r1["dens"].tolist(),
                               B=r1["B"].tolist())
        o.pop("raw1")
        o.pop("raw2", None)
        o.update(n_pass=int(sum(passed.values())), S=S, F=float(F), rank=float(sum(passed.values()) + F))
    print(f"   popolazione di {len(Ulist)} valutata in {time.time()-t0:.0f}s ({n1}+{len(jobs2)} rilassamenti)",
          flush=True)
    return outs


# ------------------------------------------------------------------------------- search
def ckpt_path(kind):
    return os.path.join(OUT, f"search_rules_{kind}.json")


def run_search(kind, budget_s=250, gens=10, lam=8):
    path = ckpt_path(kind)
    dim = len(NAMES)
    if os.path.exists(path):
        st = json.load(open(path))
    else:
        st = dict(kind=kind, gen=0, evals=[], history=[],
                  prereg_sha256=hashlib.sha256(open(os.path.join(OUT, "search_rules_prereg.json"), "rb").read()).hexdigest())
        if kind == "cma":
            st.update(mean=[0.5] * dim, sigma=0.3, C=np.eye(dim).tolist(), pc=[0.0] * dim, ps=[0.0] * dim)
    rng = np.random.default_rng(2026 if kind == "cma" else 2027)
    for _ in range(st["gen"]):  # riallinea il generatore (riproducibile alla ripresa)
        rng.standard_normal((lam, dim)) if kind == "cma" else rng.random((lam, dim))
    mu = lam // 2
    w = np.log(mu + 0.5) - np.log(np.arange(1, mu + 1))
    w /= w.sum()
    mueff = 1 / (w ** 2).sum()
    cc = (4 + mueff / dim) / (dim + 4 + 2 * mueff / dim)
    cs = (mueff + 2) / (dim + mueff + 5)
    c1 = 2 / ((dim + 1.3) ** 2 + mueff)
    cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((dim + 2) ** 2 + mueff))
    damps = 1 + 2 * max(0, np.sqrt((mueff - 1) / (dim + 1)) - 1) + cs
    chiN = np.sqrt(dim) * (1 - 1 / (4 * dim) + 1 / (21 * dim ** 2))
    t0 = time.time()
    while st["gen"] < gens and time.time() - t0 < budget_s:
        if kind == "cma":
            m, sig, C = np.array(st["mean"]), st["sigma"], np.array(st["C"])
            Dd, Bm = np.linalg.eigh(C)
            Dd = np.sqrt(np.maximum(Dd, 1e-20))
            Z = rng.standard_normal((lam, dim))
            U = np.clip(m + sig * (Z @ (Bm * Dd).T), 0, 1)
        else:
            U = rng.random((lam, dim))
        res = evaluate_population([list(u) for u in U], seed=st["gen"])
        for r in res:
            r["gen"] = st["gen"]
        st["evals"].extend(res)
        ranks = np.array([r["rank"] for r in res])
        Fs = np.array([r["F"] for r in res])
        if kind == "cma":
            order = np.argsort(-ranks, kind="stable")
            ysel = (U[order[:mu]] - m) / sig  # punti effettivamente valutati (clip => vincoli di scatola)
            yw = w @ ysel
            Cinvsqrt = Bm @ np.diag(1 / Dd) @ Bm.T
            ps = (1 - cs) * np.array(st["ps"]) + np.sqrt(cs * (2 - cs) * mueff) * Cinvsqrt @ yw
            hsig = float(np.linalg.norm(ps) / np.sqrt(1 - (1 - cs) ** (2 * (st["gen"] + 1))) / chiN
                         < 1.4 + 2 / (dim + 1))
            pc = (1 - cc) * np.array(st["pc"]) + hsig * np.sqrt(cc * (2 - cc) * mueff) * yw
            C = (1 - c1 - cmu) * C + c1 * (np.outer(pc, pc) + (1 - hsig) * cc * (2 - cc) * C) + cmu * (ysel.T * w) @ ysel
            new_sig = float(min(sig * np.exp((cs / damps) * (np.linalg.norm(ps) / chiN - 1)), 0.5))
            st.update(mean=list(np.clip(m + sig * yw, 0, 1)), sigma=new_sig, C=C.tolist(), pc=list(pc), ps=list(ps))
        allF = [r["F"] for r in st["evals"]]
        st["history"].append(dict(gen=st["gen"], best_F_gen=float(Fs.max()), best_F_sofar=float(max(allF)),
                                  mean_rank=float(ranks.mean()), max_rank=float(ranks.max()),
                                  n_feasible=int(sum(r["n_pass"] == 5 for r in res)), sigma=st.get("sigma"),
                                  t=time.time() - t0))
        st["gen"] += 1
        print(f"[{kind}] gen {st['gen']-1}: best F gen {Fs.max():.3f}, best so far {max(allF):.3f}, "
              f"feasible {sum(r['n_pass']==5 for r in res)}/{lam}, mean rank {ranks.mean():.2f} "
              f"[{time.time()-t0:.0f}s]", flush=True)
        with open(path, "w") as f:
            json.dump(st, f, default=float)
    print("fine blocco, generazioni fatte:", st["gen"])


def u_of(**kw):
    t = [kw["g"], kw["zeta"], kw["mu1"], kw["mu2"], kw["nu"], kw["kappa"], kw["a4"], kw["a6"]]
    return list((np.array(t, float) - BOX[:, 0]) / (BOX[:, 1] - BOX[:, 0]))


KNOWN = {  # teorie note dentro lo spazio delle regole (per il confronto, NON usate dalla ricerca)
    "sigma_O3_puro": dict(g=0.0, zeta=0.0, mu1=0.0, mu2=0.0, nu=-0.3, kappa=0.0, a4=0.0, a6=0.0),
    "baby_Skyrme": dict(g=0.0, zeta=0.0, mu1=0.5, mu2=0.0, nu=0.0, kappa=1.0, a4=0.0, a6=0.0),
    "sigma_gauge_Schroers": dict(g=1.0, zeta=0.0, mu1=0.0, mu2=0.7, nu=0.0, kappa=0.0, a4=0.0, a6=0.0),
    "baby_Skyrme_gauge_GPS": dict(g=1.0, zeta=0.0, mu1=0.5, mu2=0.0, nu=0.0, kappa=1.0, a4=0.0, a6=0.0),
    "Higgs_abeliano_NO": dict(g=1.0, zeta=0.0, mu1=0.0, mu2=0.0, nu=1.0, kappa=0.0, a4=0.0, a6=0.0),
}


# ------------------------------------------------------------------------------- dynamics
LAM_DYN = 20.0  # solo per tenere |phi| ~ 1 nella dinamica: l'energia dipende da n = phi/|phi|


def dynamics(u, tc_mult=100, perturb=0.0, seed=5, dt_frac=0.25, label="", budget_s=480, steps_relax=3000):
    """dinamica hamiltoniana vera nella gauge temporale A_0 = 0, RK4.
    H = sum h^2 [ 1/2 pi^T M^-1 pi + D^2/(2 eps(B)) ] + E_statica,
    M = 1 + kappa^2 sum_i (|D_i phi|^2 - D_i phi D_i phi^T)   (parte cinetica del termine di Skyrme)
    eps(B) = u'(B)/B = 1 + 4 a4 B^2 + 6 a6 B^4  (esatto all'ordine D^2; l'oggetto parte fermo, D = 0).
    Legge di Gauss: D = 0 e pi = 0 all'inizio => rho = 0, div D = 0; conservata dalla dinamica (H gauge-invariante).
    Spugna (smorzamento dei momenti) nelle ultime 2 unita' vicino al bordo: assorbe la radiazione uscente."""
    torch = _torch()
    P = decode(np.asarray(u, float))
    vac = vacuum(P)
    phi, A = init_config(1, vac)
    ck = os.path.join(OUT, f"search_rules_dynckpt_{label}_p{perturb:g}.npz")
    if os.path.exists(ck):  # ripresa: stato salvato alla fine del blocco precedente
        z = np.load(ck, allow_pickle=True)
        r1 = dict(phi=z["phi_static"], A=z["A_static"])
        a1 = json.loads(str(z["a1"]))
    else:
        r1 = flow([phi], [A], [P], [vac["Vmin"]], torch, steps_relax)[0]
        a1 = analyse(r1, 1, vac)
        a1["drift"] = r1["drift"]
    X, Y = grid()
    N, h = N_GRID, H
    bmask = np.zeros((N, N), bool)
    bmask[0, :] = bmask[-1, :] = bmask[:, 0] = bmask[:, -1] = True
    phi0, A0 = r1["phi"].copy(), r1["A"].copy()
    rng = np.random.default_rng(seed)
    if perturb > 0:
        phi0 = phi0 * (1 + perturb * rng.standard_normal(phi0.shape))
        phi0[:, bmask] = r1["phi"][:, bmask]
        A0 = A0 + perturb * (np.abs(A0).max() + 1e-3) * rng.standard_normal(A0.shape)
    dist = np.minimum.reduce([X - X.min(), X.max() - X, Y - Y.min(), Y.max() - Y])
    gam = torch.tensor(2.0 * np.clip(1 - dist / 2.0, 0, 1) ** 2)
    bm = torch.tensor(bmask)
    Pt = params_tensor([P], torch)
    vmin = torch.tensor([[[vac["Vmin"]]]])
    kap2 = P["kappa"] ** 2
    eye = torch.eye(3)[:, :, None, None]

    def Hfun(phi, A, pi, D):
        T, d = energy_terms(phi[None], A[None], Pt, vmin, torch)
        Es = sum(T.values())[0] + LAM_DYN * h * h * (((phi ** 2).sum(0) - 1) ** 2).sum()
        Dx = torch.nn.functional.pad(d["Dx"][0], (0, 0, 0, 1))
        Dy = torch.nn.functional.pad(d["Dy"][0], (0, 1, 0, 0))
        G = torch.zeros(3, 3, N, N)
        for Dd in (Dx, Dy):
            G = G + (Dd ** 2).sum(0)[None, None] * eye - Dd[:, None] * Dd[None, :]
        M = eye + kap2 * G
        v = torch.linalg.solve(M.permute(2, 3, 0, 1), pi.permute(1, 2, 0)[..., None])[..., 0]
        Tm = 0.5 * (pi.permute(1, 2, 0) * v).sum()
        B2 = d["B"][0] ** 2
        epsB = 1 + 4 * P["a4"] * B2 + 6 * P["a6"] * B2 ** 2
        ex = torch.nn.functional.pad(epsB, (1, 1, 0, 0), value=1.0)
        epx = 0.5 * (ex[:, :-1] + ex[:, 1:])
        ey = torch.nn.functional.pad(epsB, (0, 0, 1, 1), value=1.0)
        epy = 0.5 * (ey[:-1, :] + ey[1:, :])
        TE = 0.5 * (D[0, :-1, :] ** 2 / epx).sum() + 0.5 * (D[1, :, :-1] ** 2 / epy).sum()
        return Es + h * h * (Tm + TE)

    def rhs(st):
        phi, A, pi, D = [x.detach().requires_grad_(True) for x in st]
        Ht = Hfun(phi, A, pi, D)
        gphi, gA, gpi, gD = torch.autograd.grad(Ht, [phi, A, pi, D])
        s2 = h * h
        zero = torch.zeros(())
        dphi = torch.where(bm, zero, gpi / s2)
        dpi = torch.where(bm, zero, -gphi / s2 - gam * pi)
        return [dphi, gD / s2, dpi, -gA / s2 - gam * D], float(Ht)

    st = [torch.tensor(phi0), torch.tensor(A0), torch.zeros(3, N, N), torch.zeros(2, N, N)]
    rE = a1["rE"]
    tc = 2 * max(rE, 0.5)
    tmax = tc_mult * tc
    w = np.sqrt((8 / h ** 2) * (2 + 4 * kap2 + 2 * abs(P["zeta"])) + 16 * LAM_DYN)
    dt = min(dt_frac * h, 1.2 / w)
    nst = int(np.ceil(tmax / dt))
    xc, yc = a1["centroid"]
    disk = np.hypot(X - xc, Y - yc) < L_BOX / 4
    rec = []
    t = 0.0
    n_start = 0
    if os.path.exists(ck):
        z = np.load(ck, allow_pickle=True)
        st = [torch.tensor(z[k]) for k in ("s0", "s1", "s2", "s3")]
        rec = json.loads(str(z["rec"]))
        t = float(z["t"])
        n_start = int(z["n"])
    t0 = time.time()
    every = max(1, int(round(0.5 * tc / dt)))
    E_disk0 = rec[0]["E_disk"] if rec else None
    H0 = rec[0]["H"] if rec else None
    finished = False
    for n in range(n_start, nst + 1):
        if n % every == 0:
            with torch.no_grad():
                T, d = energy_terms(st[0][None], st[1][None], Pt, vmin, torch)
            dens = density_map(d, 0)
            pis, Dn = st[2].numpy(), st[3].numpy()
            dens = dens + 0.5 * (pis ** 2).sum(0)  # (mappa: cinetica con M ~ 1)
            dens[:-1, :] += 0.5 * Dn[0, :-1, :] ** 2
            dens[:, :-1] += 0.5 * Dn[1, :, :-1] ** 2
            Ed = float(dens[disk].sum() * h * h)
            E_disk0 = Ed if E_disk0 is None else E_disk0
            Qn = topo_charge(st[0].numpy())[0] if vac["rho"] == 0.0 else winding_boundary(st[0].numpy())
            _, Hn = rhs(st)
            H0 = Hn if H0 is None else H0
            rec.append(dict(t=t, E_disk_frac=Ed / E_disk0, E_disk=Ed, Q=float(Qn), H=Hn))
            if n % (every * 20) == 0:
                print(f"{label} t={t:.1f}/{tmax:.1f} Edisk/E0={Ed/E_disk0:.4f} Q={Qn:.3f} H={Hn:.5f} "
                      f"[{time.time()-t0:.0f}s]", flush=True)
            if time.time() - t0 > budget_s:
                print("budget di tempo esaurito: salvo e riprendo alla prossima esecuzione")
                np.savez(ck, s0=st[0].numpy(), s1=st[1].numpy(), s2=st[2].numpy(), s3=st[3].numpy(),
                         rec=json.dumps(rec), t=t, n=n, phi_static=r1["phi"], A_static=r1["A"],
                         a1=json.dumps(a1, default=float))
                return None
        k1, _ = rhs(st)
        k2, _ = rhs([s + 0.5 * dt * k for s, k in zip(st, k1)])
        k3, _ = rhs([s + 0.5 * dt * k for s, k in zip(st, k2)])
        k4, _ = rhs([s + dt * k for s, k in zip(st, k3)])
        st = [(s + dt / 6 * (a + 2 * b + 2 * c + e)).detach() for s, a, b, c, e in zip(st, k1, k2, k3, k4)]
        t += dt
    fr = np.array([r["E_disk_frac"] for r in rec])
    return dict(u=list(map(float, u)), P={k: float(v) for k, v in P.items()}, perturb=perturb, tc=tc, tmax=tmax,
                t_reached=rec[-1]["t"], dt=dt, min_frac=float(fr.min()), end_frac=float(fr[-1]),
                Q_start=rec[0]["Q"], Q_end=rec[-1]["Q"], static=a1,
                C2_dyn=bool(rec[-1]["t"] >= tmax - 0.5 * tc and fr.min() >= 0.9 and abs(rec[-1]["Q"] - rec[0]["Q"]) < 0.5),
                rec=rec)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "prereg":
        p = os.path.join(OUT, "search_rules_prereg.json")
        if os.path.exists(p):
            print("pre-registro gia' esistente: non lo riscrivo")
        else:
            with open(p, "w", encoding="utf8") as f:
                json.dump(PREREG, f, indent=1, ensure_ascii=False)
        print("sha256", hashlib.sha256(open(p, "rb").read()).hexdigest())
    elif cmd in ("cma", "rand"):
        run_search(cmd)
    elif cmd == "known":
        names = list(KNOWN)
        res = evaluate_population([u_of(**KNOWN[k]) for k in names], seed=99, keep_fields=False)
        out = {k: r for k, r in zip(names, res)}
        with open(os.path.join(OUT, "search_rules_known.json"), "w") as f:
            json.dump(out, f, default=float)
        for k, r in out.items():
            print(k, "pass", r["passed"], "F", round(r["F"], 3), "E1", r.get("Q1", {}).get("E"),
                  "Q", r.get("Q1", {}).get("Q"), "fEM", r.get("Q1", {}).get("fEM"),
                  "spec", {q: (s.get("ratio"), s.get("ok")) for q, s in (r.get("spectrum") or {}).items()})
    elif cmd == "dyn":
        # dyn <file.json con chiave u> <etichetta> <perturb>
        src = json.load(open(sys.argv[2]))
        lab = sys.argv[3]
        pert = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
        r = dynamics(src["u"], perturb=pert, label=lab)
        if r is None:
            sys.exit(3)
        with open(os.path.join(OUT, f"search_rules_dyn_{lab}_p{pert:g}.json"), "w") as f:
            json.dump(r, f, default=float)
        print(json.dumps({k: v for k, v in r.items() if k not in ("rec", "static")}, default=float, indent=1))
