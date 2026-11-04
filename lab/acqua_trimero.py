"""acqua_trimero.py - trimero ciclico dell'acqua (H2O)3, stesso motore del dimero (Cholesky nostro + FNO + CCSD(T)).
Monomeri RIGIDI alla nostra geometria CBS. Ossigeni su un triangolo equilatero di lato R (piano xy); il monomero k
dona il suo H legato all'O del monomero k+1. Coordinate (3, simmetria C3 approssimata, schema "su su giu'"):
  R  = lato O-O;   td = rotazione nel piano dell'O-H legato rispetto all'asse O_k -> O_k+1 (+ = verso l'esterno);
  ps = angolo del piano del monomero attorno all'O-H legato (0 = H libero nel piano dell'anello, verso l'esterno;
       90 = H libero perpendicolare all'anello), con segni +,+,- per i tre H liberi.
Counterpoise a molti corpi (tutto nella base del trimero): tipi ABC, AB, AC, BC, A, B, C.
  E_int = E_ABC - E_A - E_B - E_C ;  E_3corpi = E_ABC - E_AB - E_AC - E_BC + E_A + E_B + E_C
Uso: python acqua_trimero.py calc <tipo> <X> <R> <td> <ps> [tau_no]   (a pezzi, exit 3 = non finito)"""
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import acqua_dimero as DM, acqua_basis as AB

PD = os.path.join(AB.OUT, "acqua_trimero_punti"); os.makedirs(PD, exist_ok=True)
SGN = (1, 1, -1)


def trimer_xyz(R, td, ps):
    r, th = DM.mono_geom(); r *= DM.A2B; R *= DM.A2B; t0 = math.radians(th)
    rho = R / math.sqrt(3); z = np.array([0, 0, 1.0])
    O = [rho * np.array([math.cos(math.radians(90 + 120 * k)), math.sin(math.radians(90 + 120 * k)), 0]) for k in range(3)]
    atoms, xyz = [], []
    for k in range(3):
        a = O[(k + 1) % 3] - O[k]; a /= np.linalg.norm(a)
        out = O[k] + 0.5 * (O[(k + 1) % 3] - O[k]); out /= np.linalg.norm(out)     # verso esterno, a meta' lato
        q = math.radians(td) * (1 if np.dot(np.cross(a, out), z) > 0 else -1)      # + ruota verso l'esterno
        e = math.cos(q) * a + math.sin(q) * np.cross(z, a)
        p = np.cross(z, e); p = p if np.dot(p, out) > 0 else -p                      # nel piano, verso l'esterno
        s = math.radians(ps)
        f = math.cos(t0) * e + math.sin(t0) * (math.cos(s) * p + SGN[k] * math.sin(s) * z)
        atoms += ["O", "H", "H"]; xyz += [O[k], O[k] + r * e, O[k] + r * f]
    return atoms, np.array(xyz)


def mp2(kind, X, R, td, ps):
    """solo RHF + MP2 completo (niente FNO/CC): per la geometria. Un pezzo unico, niente checkpoint."""
    import time, acqua_int as I, acqua_chol as CH
    k = DM.key(kind + "m", X, R, td, ps); f = os.path.join(PD, k + ".json")
    if os.path.exists(f):
        return json.load(open(f))
    t0 = time.time(); atoms, xyz = trimer_xyz(R, td, ps); real = [int(REAL[kind][i // 3]) for i in range(9)]
    B = DM.basis(atoms, xyz, X, real); nocc = 5 * REAL[kind].count("1")
    ints = I.integrals(B, eri=False); L = CH.cholesky(B); r = CH.scf_df(B, ints, L, nocc)
    Bm = CH.mo_B(L, r["C"]); nch = len(L); del L
    E2 = CH.mp2_fno(Bm, r["eps"], nocc, 1e-6)[2]
    o = dict(kind=kind, X=X, R=R, td=td, ps=ps, M=B.M, nch=nch, E_hf=r["E"], E_mp2=E2, E_tot=r["E"] + E2,
             scf_conv=r["conv"], t=time.time() - t0)
    json.dump(o, open(f, "w"), indent=1)
    return o


REAL = {"ABC": "111", "AB": "110", "AC": "101", "BC": "011", "A": "100", "B": "010", "C": "001"}


def calc(kind, X, R, td, ps, tau_no=DM.TAU_NO):
    k = DM.key(kind, X, R, td, ps, f"{tau_no:.0e}"); f = os.path.join(PD, k + ".json")
    if os.path.exists(f):
        return json.load(open(f))
    atoms, xyz = trimer_xyz(R, td, ps)
    real = [int(REAL[kind][i // 3]) for i in range(9)]
    B = DM.basis(atoms, xyz, X, real); ck = os.path.join(DM.AMP, "tri_" + k + ".npz")
    o = DM.solve(B, 5 * REAL[kind].count("1"), ck, tau_no)
    if o is None:
        return None
    o.update(kind=kind, X=X, R=R, td=td, ps=ps, tau_no=tau_no)
    json.dump(o, open(f, "w"), indent=1); os.remove(ck)
    return o


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "geo":
        at, x = trimer_xyz(*map(float, a[1:4])); x = x / DM.A2B
        for i in range(9):
            print(at[i], *[f"{v:8.4f}" for v in x[i]])
        d = lambda i, j: np.linalg.norm(x[i] - x[j])
        print("O-O", d(0, 3), d(3, 6), d(6, 0), " H-bond H...O", d(1, 3), d(4, 6), d(7, 0), " min H-H", min(d(i, j) for i in range(9) for j in range(i) if at[i] == at[j] == "H"))
        sys.exit()
    if a[0] == "mp2":
        print(json.dumps(mp2(a[1], int(a[2]), float(a[3]), float(a[4]), float(a[5])))); sys.exit()
    o = calc(a[1], int(a[2]), float(a[3]), float(a[4]), float(a[5]), float(a[6]) if len(a) > 6 else DM.TAU_NO)
    if o is None:
        print("PARZIALE"); sys.exit(3)
    print(json.dumps(o))
