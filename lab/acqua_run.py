"""acqua_run.py - un punto di calcolo (molecola o atomo), checkpoint JSON in lab/out/acqua_punti/.
Molecola H2O: O nell'origine, H1 = r1 (sin(t/2), 0, cos(t/2)), H2 = r2 (-sin(t/2), 0, cos(t/2)), asse C2 = z.
Uso:  python acqua_run.py mol X r1_A r2_A theta_deg [campo_au] [noaug]
      python acqua_run.py atom O|H X [noaug]
      python acqua_run.py h2 X R_bohr
Nuclei fermi (Born-Oppenheimer), niente relativita'. Tutti gli elettroni correlati (CCSD(T))."""
import os, sys, json, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import acqua_int as I, acqua_qc as Q, acqua_basis as AB

A2B = 1.0 / 0.529177210903          # bohr per angstrom (costante fisica)
PD = os.path.join(AB.OUT, "acqua_punti"); os.makedirs(PD, exist_ok=True)


def geom(r1, r2, th):
    t = math.radians(th) / 2
    return [[0, 0, 0], [r1 * A2B * math.sin(t), 0, r1 * A2B * math.cos(t)], [-r2 * A2B * math.sin(t), 0, r2 * A2B * math.cos(t)]]


def key(kind, *a):
    return kind + "_" + "_".join(f"{x:.6f}" if isinstance(x, float) else str(x) for x in a)


def compute(atoms, coords, X, aug, nd, ns, field=0.0, corr=True, extra=None):
    t0 = time.time()
    cen, Z, sh = AB.shells_for(atoms, coords, X, aug)
    B = I.Basis(cen, Z, sh); ints = I.integrals(B); t_int = time.time() - t0
    r = Q.scf(B, ints, nd, ns, field=field)
    out = dict(atoms=atoms, coords=np.asarray(coords).tolist(), X=X, aug=aug, field=field, M=B.M, E_hf=r["E"],
               scf_conv=bool(r["conv"]), t_int=t_int)
    # HF dipole (electronic + nuclear), for information
    C = r["C"]; Da = C[:, :nd + ns] @ C[:, :nd + ns].T; Db = C[:, :nd] @ C[:, :nd].T
    out["mu_hf_z"] = float(np.sum(B.Z * B.R[:, 2]) - np.sum((Da + Db) * ints["D"][2]))
    if corr and nd + ns > 0:
        h, G = Q.mo_ints(B, ints, r)
        if ns == 0:
            Ec, t1, t2, conv = Q.rccsd(G, r["eps"], nd)
            ET = Q.rccsd_t(G, r["eps"], nd, t1, t2)
            out.update(Ec_ccsd=Ec, Ec_T=ET, cc_conv=bool(conv), t1diag=float(np.linalg.norm(t1) / math.sqrt(2 * nd)))
        else:
            o = Q.so_ccsd_t(h, G, nd, ns)
            out.update(Ec_ccsd=o["Ec_ccsd"], Ec_T=o["Ec_T"], cc_conv=bool(o["conv"]),
                       dE_ref=o["E_ref"] + r["Enuc"] - r["E"])
        out["E_tot"] = out["E_hf"] + out["Ec_ccsd"] + out["Ec_T"]
    out["t"] = time.time() - t0
    if extra:
        out.update(extra)
    return out


def atom_calc(el, X, par=None, aug=True, corr=True, triples=True, ck=None, t_budget=None):
    """spherical ROHF (O: s,p block, O_h-averaged Fock; H: s block), functions of higher l appended as virtuals
    (Loewdin-orthonormalised: same centre, other l -> already orthogonal to the occupied); all-electron CC"""
    t0 = time.time()
    par = par or AB.load()[f"{el}{X}"]
    lmax_occ = 1 if el == "O" else 0
    sh = [(0, L, a) for L, a in AB.exps(el, X, par, aug)]
    lo = [s for s in sh if s[1] <= lmax_occ]; hi = [s for s in sh if s[1] > lmax_occ]
    B = I.Basis([[0, 0, 0]], [AB.ZNUC[el]], lo + hi); Blo = I.Basis([[0, 0, 0]], [AB.ZNUC[el]], lo)
    ints = I.integrals(B, eri=corr); m = Blo.M
    il = dict(S=ints["S"][:m, :m], T=ints["T"][:m, :m], V=ints["V"][:m, :m], D=ints["D"][:, :m, :m])
    il["G"] = ints["G"][:m, :m, :m, :m].copy() if corr else I.integrals(Blo)["G"]
    nd, ns = (3, 2) if el == "O" else (0, 1)
    r = Q.scf(Blo, il, nd, ns, sym=Q.sph_sym(Blo) if el == "O" else None)
    out = dict(atoms=[el], X=X, aug=aug, M=B.M, E_hf=r["E"], scf_conv=bool(r["conv"]))
    if corr and nd > 0:
        Sh = ints["S"][m:, m:]; w, U = np.linalg.eigh(Sh); k = w > 1e-8
        nlo = r["C"].shape[1]
        C = np.zeros((B.M, nlo + int(k.sum()))); C[:m, :nlo] = r["C"]; C[m:, nlo:] = U[:, k] / np.sqrt(w[k])
        hfull = ints["T"] + ints["V"]
        hm, G = Q.mo_ints(B, ints, dict(C=C, h=hfull))
        o = Q.so_ccsd_t(hm, G, nd, ns, triples=triples, ck=ck, t_budget=t_budget)
        if o.get("parziale"):
            return dict(parziale=True)
        out.update(Ec_ccsd=o["Ec_ccsd"], Ec_T=o["Ec_T"], cc_conv=bool(o["conv"]), dE_ref=o["E_ref"] - r["E"])
        out["E_tot"] = out["E_hf"] + out["Ec_ccsd"] + out["Ec_T"]
    elif el == "H":
        out.update(Ec_ccsd=0.0, Ec_T=0.0, E_tot=r["E"], cc_conv=True)
    out["t"] = time.time() - t0
    return out


def run(argv):
    kind = argv[0]
    if kind == "mol":
        X = int(argv[1]); r1, r2, th = map(float, argv[2:5])
        field = float(argv[5]) if len(argv) > 5 and argv[5] != "noaug" else 0.0
        aug = "noaug" not in argv
        k = key("mol", X, r1, r2, th, field, int(aug))
        f = os.path.join(PD, k + ".json")
        if os.path.exists(f):
            return json.load(open(f))
        out = compute(["O", "H", "H"], geom(r1, r2, th), X, aug, 5, 0, field, extra=dict(r1=r1, r2=r2, theta=th))
    elif kind == "atom":
        el, X = argv[1], int(argv[2]); aug = "noaug" not in argv
        k = key("atom", el, X, int(aug)); f = os.path.join(PD, k + ".json")
        if os.path.exists(f):
            return json.load(open(f))
        os.makedirs(os.path.join(AB.OUT, "acqua_amp"), exist_ok=True)
        out = atom_calc(el, X, aug=aug, ck=os.path.join(AB.OUT, "acqua_amp", k + ".npz"), t_budget=480)
        if out.get("parziale"):
            print("PARZIALE"); sys.exit(3)
    elif kind == "h2":
        X = int(argv[1]); R = float(argv[2]); aug = "noaug" not in argv
        k = key("h2", X, R, int(aug)); f = os.path.join(PD, k + ".json")
        if os.path.exists(f):
            return json.load(open(f))
        out = compute(["H", "H"], [[0, 0, 0], [0, 0, R]], X, aug, 1, 0, extra=dict(R=R))
    json.dump(out, open(f, "w"), indent=1)
    return out


if __name__ == "__main__":
    o = run(sys.argv[1:])
    print(json.dumps({k: v for k, v in o.items() if k != "coords"}))
