"""acqua_dimero.py - dimero dell'acqua (H2O)2 con integrali di Cholesky nostri + FNO + CCSD(T) (studio ACQUA).
Monomeri RIGIDI alla nostra geometria CBS (lab/out/acqua_stato.json). Piano di simmetria = xy.
  donatore A: O_A nell'origine, H legato a angolo td da +x, l'altro H a td + theta0 (nel piano xy)
  accettore B: O_B = (R, 0, 0); bisettrice (O -> centro degli H) a angolo ta da +x nel piano xy, H a +-z
Counterpoise (Boys-Bernardi): E_int = E_AB(base AB) - E_A(base AB) - E_B(base AB); D_e = -E_int (monomeri rigidi:
rilassamento dei monomeri trascurato, debito).
Uso: python acqua_dimero.py calc <kind AB|A|B> <X> <R_A> <td> <ta> [tau_no]     (a pezzi, exit 3 = non finito)
     python acqua_dimero.py mono <X> [tau_no]        controllo: monomero alla geometria CBS contro il calcolo esatto"""
import os, sys, json, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import acqua_int as I, acqua_qc as Q, acqua_basis as AB, acqua_run as AR, acqua_chol as CH

PD = os.path.join(AB.OUT, "acqua_dimero_punti"); os.makedirs(PD, exist_ok=True)
AMP = os.path.join(AB.OUT, "acqua_amp"); os.makedirs(AMP, exist_ok=True)
TAU_NO = 1e-6
A2B = AR.A2B


def mono_geom():
    s = json.load(open(os.path.join(AB.OUT, "acqua_stato.json")))
    return s["min_cbs"]


def dimer_xyz(R, td, ta):
    r, th = mono_geom(); r *= A2B; R *= A2B
    t0 = math.radians(th); td = math.radians(td); ta = math.radians(ta)
    OA = np.zeros(3)
    HA1 = r * np.array([math.cos(td), math.sin(td), 0]); HA2 = r * np.array([math.cos(td + t0), math.sin(td + t0), 0])
    OB = np.array([R, 0, 0]); u = np.array([math.cos(ta), math.sin(ta), 0]); z = np.array([0, 0, 1.0])
    HB1 = OB + r * (math.cos(t0 / 2) * u + math.sin(t0 / 2) * z); HB2 = OB + r * (math.cos(t0 / 2) * u - math.sin(t0 / 2) * z)
    return ["O", "H", "H", "O", "H", "H"], np.array([OA, HA1, HA2, OB, HB1, HB2])


def basis(atoms, xyz, X, real):
    bas = AB.load(); sh = []
    for i, el in enumerate(atoms):
        sh += [(i, L, a) for L, a in AB.exps(el, X, bas[f"{el}{X}"], True)]
    Z = [AB.ZNUC[el] if real[i] else 0 for i, el in enumerate(atoms)]
    return I.Basis(xyz, Z, sh)


def solve(B, nocc, ck, tau_no, t_budget=520):
    """Cholesky -> RHF -> FNO -> CCSD -> (T), resumable through ck (npz); returns dict or None if unfinished"""
    t_end = time.time() + t_budget; t0 = time.time()
    z = {}
    if os.path.exists(ck):
        try:
            z = dict(np.load(ck, allow_pickle=False))
        except Exception:                       # truncated checkpoint (e.g. disk full): start again
            os.remove(ck); z = {}
    if "Bn" not in z:
        ints = I.integrals(B, eri=False)
        L = CH.cholesky(B); nch = len(L)
        r = CH.scf_df(B, ints, L, nocc)
        Bm = CH.mo_B(L, r["C"]); del L
        if tau_no > 0:
            Bn, en, E2, E2t, k, occ = CH.mp2_fno(Bm, r["eps"], nocc, tau_no)
        else:
            Bn, en = Bm, r["eps"]; E2 = E2t = 0.0; k = Bm.shape[1] - nocc
        z = dict(Bn=Bn, en=en, E_hf=r["E"], E2=E2, E2t=E2t, k=k, nch=nch, scf_conv=r["conv"], t_setup=time.time() - t0)
        np.savez(ck, **z)
        if time.time() > t_end - 120:
            return None
    G = CH.DFG(z["Bn"]); en = z["en"]
    if not bool(z.get("cc_conv", False)):
        ti = (z["t1"], z["t2"]) if "t1" in z else None
        Ec, t1, t2, conv = Q.rccsd(G, en, nocc, t_init=ti, t_end=t_end)
        z.update(t1=t1, t2=t2, cc_conv=conv, Ec_ccsd=Ec); np.savez(ck, **z)
        if not conv:
            return None
    st = (float(z["TE"]), int(z["Ti"])) if "TE" in z else None
    ET = Q.rccsd_t(G, en, nocc, z["t1"], z["t2"], start=st, t_end=t_end)
    if isinstance(ET, tuple):
        z.update(TE=ET[0], Ti=ET[1]); np.savez(ck, **z); return None
    dmp2 = float(z["E2"] - z["E2t"])
    return dict(M=B.M, nch=int(z["nch"]), n_virt_fno=int(z["k"]), nocc=nocc, E_hf=float(z["E_hf"]),
                Ec_ccsd=float(z["Ec_ccsd"]), Ec_T=float(ET), dMP2_fno=dmp2, E_mp2_full=float(z["E2"]),
                E_tot=float(z["E_hf"] + z["Ec_ccsd"] + ET + dmp2), scf_conv=bool(z["scf_conv"]))


def key(*a):
    return "_".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in a)


def calc(kind, X, R, td, ta, tau_no=TAU_NO):
    k = key(kind, X, R, td, ta, f"{tau_no:.0e}"); f = os.path.join(PD, k + ".json")
    if os.path.exists(f):
        return json.load(open(f))
    atoms, xyz = dimer_xyz(R, td, ta)
    real = {"AB": [1] * 6, "A": [1, 1, 1, 0, 0, 0], "B": [0, 0, 0, 1, 1, 1]}[kind]
    B = basis(atoms, xyz, X, real)
    o = solve(B, 10 if kind == "AB" else 5, os.path.join(AMP, "dim_" + k + ".npz"), tau_no)
    if o is None:
        return None
    o.update(kind=kind, X=X, R=R, td=td, ta=ta, tau_no=tau_no)
    json.dump(o, open(f, "w"), indent=1)
    os.remove(os.path.join(AMP, "dim_" + k + ".npz"))          # checkpoint no longer needed (disk)
    return o


def mono(X, tau_no=TAU_NO):
    r, th = mono_geom()
    f = os.path.join(PD, key("mono", X, f"{tau_no:.0e}") + ".json")
    if os.path.exists(f):
        return json.load(open(f))
    cen, Z, sh = AB.shells_for(["O", "H", "H"], AR.geom(r, r, th), X, aug=True)
    B = I.Basis(cen, Z, sh)
    o = solve(B, 5, os.path.join(AMP, "dim_" + key("mono", X, f"{tau_no:.0e}") + ".npz"), tau_no)
    if o is None:
        return None
    ex = json.load(open(os.path.join(AR.PD, AR.key("mol", X, float(r), float(r), float(th), 0.0, 1) + ".json")))
    o.update(esatto=dict(E_hf=ex["E_hf"], Ec=ex["Ec_ccsd"] + ex["Ec_T"]),
             dE_hf=o["E_hf"] - ex["E_hf"], dEc=o["E_tot"] - o["E_hf"] - ex["Ec_ccsd"] - ex["Ec_T"])
    json.dump(o, open(f, "w"), indent=1)
    return o


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "calc":
        o = calc(a[1], int(a[2]), float(a[3]), float(a[4]), float(a[5]), float(a[6]) if len(a) > 6 else TAU_NO)
    else:
        o = mono(int(a[1]), float(a[2]) if len(a) > 2 else TAU_NO)
    if o is None:
        print("PARZIALE"); sys.exit(3)
    print(json.dumps(o))
