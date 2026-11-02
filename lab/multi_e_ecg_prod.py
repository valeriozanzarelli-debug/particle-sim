"""multi_e_ecg_prod.py - produzione ECG dell'elio (pre-registrazione v2, sbloccata solo se G1' PASS).
  nr      corse non relativistiche a m_kin(theta_A): he1S Z = 1, 2, 6..10; he3S Z = 2; 1 elettrone Z = 1, 2, 10
  full    operatore pieno della banda (existing_A = livello 1 in uso, relearn_A = livello ri-imparato):
          base NR + 2 giri di raffinamento con T = S sum F(s_i); primo ordine sullo stato NR; G4 (parte oltre la zona)
Uscite: lab/out/multi_e_ecg.json (chiave 'prod'), basi in lab/out/multi_e_ecg_<job>.npz (non committate).
Run: python lab/multi_e_ecg_prod.py nr ; python lab/multi_e_ecg_prod.py full"""
import os, sys, time
import numpy as np
from multi_e_core import OUT, load_json, save_json
from multi_e_ecg import System, SVM, run_svm, JS

E0 = load_json(JS)
if not E0.get("G1_v2", {}).get("PASS"):
    raise SystemExit("G1' non PASS: produzione ECG ferma (pre-registrazione v2)")
KER = load_json("multi_e_kernel.json")
m_A = KER["R3"]["m_kin_theta_A_used"]
SEED = 2026

JOBS = {"he1S_Z1": dict(N=2, nup=1, Z=1.0, K=200, Lmax=30.0)}
for Zi in (2, 6, 7, 8, 9, 10):
    JOBS[f"he1S_Z{Zi}"] = dict(N=2, nup=1, Z=float(Zi), K=200, Lmax=20.0)
JOBS["he3S_Z2"] = dict(N=2, nup=2, Z=2.0, K=150, Lmax=20.0)
for Zi in (1, 2, 10):
    JOBS[f"one_Z{Zi}"] = dict(N=1, nup=1, Z=float(Zi), K=40, Lmax=20.0)

FT = np.load(os.path.join(OUT, "multi_e_Ftab.npz"))


def band(name):
    return (FT[f"{name}__s"], FT[f"{name}__R"]), float(FT[f"{name}__m_kin"]), FT[f"{name}__Fb"]


class FbSystem(System):
    """stesso operatore, ma la cinetica e' solo la parte da |k| oltre la zona (per G4)"""
    def __init__(self, *a, Fb=None, **k):
        super().__init__(*a, **k); self.Fb = Fb

    def F(self, s):
        sv = self.Ftab[0]
        return np.interp(np.log(np.clip(s, sv[0], sv[-1])), np.log(sv), self.Fb)


def expect(sv, sys_, c=None):
    """valori attesi delle parti con l'operatore sys_ sullo stato di coefficienti c (default: stato di sv)"""
    K = len(sv.A)
    if c is None:
        _, c = sv.state()
    I, J = np.triu_indices(K)
    S, H, P = sys_.me(sv.A[I], sv.A[J], parts=True)
    w = c[I] * c[J] * sv.nset[I] * sv.nset[J] * np.where(I == J, 1.0, 2.0)
    nrm = float(np.sum(w * S))
    out = {k: float(np.sum(w * v)) / nrm for k, v in P.items()}
    out["E"] = float(np.sum(w * H)) / nrm
    return out


def load_basis(name, j, sys_):
    z = np.load(os.path.join(OUT, f"multi_e_ecg_{name}.npz"))
    sv = SVM(sys_, np.random.default_rng(SEED + 1), j["Z"], j["Lmax"], 30)
    sv.A = z["A"]
    sv.rebuild()
    return sv


what = sys.argv[1] if len(sys.argv) > 1 else "nr"
names = sys.argv[2].split(",") if len(sys.argv) > 2 else None
res = load_json(JS)
prod = res.get("prod", {})

if what == "nr":
    for name, j in JOBS.items():
        if names and name not in names:
            continue
        if name in prod.get("nr", {}):
            continue
        t0 = time.time()
        s = System(j["N"], j["nup"], j["Z"], m_A, name=name)
        sv, EK, Es = run_svm(s, j["K"], j["Z"], SEED, Lmax=j["Lmax"], log=lambda *a: None)
        E = sv.rebuild()
        rec = dict(E=E, E_Kminus50=sv.energy_subset(j["K"] - 50) if j["K"] > 50 else sv.energy_subset(j["K"] - 10),
                   sweeps=Es, EK=EK[-5:], parts=expect(sv, s), m_kin=m_A, K=j["K"], cond_S=float(np.linalg.cond(sv.Smat)),
                   sec=time.time() - t0)
        rec["conv_err"] = abs(rec["E"] - rec["E_Kminus50"])
        np.savez(os.path.join(OUT, f"multi_e_ecg_{name}.npz"), A=sv.A)
        res = load_json(JS); prod = res.get("prod", {}); prod.setdefault("nr", {})[name] = rec; res["prod"] = prod; save_json(JS, res)
        print(name, f"E = {E:.10f}  conv {rec['conv_err']:.1e}  {rec['sec']:.0f} s", flush=True)

if what == "full":
    targets = [("existing_A", n) for n in ("one_Z10", "he1S_Z2", "he1S_Z6", "he1S_Z7", "he1S_Z8", "he1S_Z9", "he1S_Z10")] + \
              [("relearn_A", n) for n in ("one_Z10", "he1S_Z2", "he1S_Z6", "he1S_Z7", "he1S_Z8", "he1S_Z9", "he1S_Z10")]
    for bn, name in targets:
        key = f"{bn}:{name}"
        if names and key not in names:
            continue
        if key in prod.get("full", {}):
            continue
        t0 = time.time()
        j = JOBS[name]
        Ftab, m_b, Fb = band(bn)
        s_nr = System(j["N"], j["nup"], j["Z"], m_A, name=name)
        s_f = System(j["N"], j["nup"], j["Z"], m_A, kin="full", Ftab=Ftab, name=name)
        sv = load_basis(name, j, s_nr)
        E_nr_A = sv.E
        E_nr_b = E_nr_A * m_b / m_A                                   # NR alla massa della banda (scala esatta)
        p_nr = expect(sv, s_nr)
        pf1 = expect(sv, s_f)                                         # operatore pieno sullo stato NR (primo ordine)
        first = pf1["E"] - E_nr_b                                     # <H_full>_NR - E_NR(m_b)
        sfb = FbSystem(j["N"], j["nup"], j["Z"], m_A, kin="full", Ftab=Ftab, Fb=Fb, name=name)
        Tb = expect(sv, sfb)["T"]
        mv1 = -(7.2973525693e-3 ** 2) * p_nr["P4"] / (8 * m_b ** 3)
        # ricetta variazionale pre-registrata (base NR, poi operatore pieno): si prova e si registra se collassa
        E_f0 = sv.rebuild(s_f)
        collapsed = bool(E_f0 < E_nr_b - 10 * abs(first) - 1.0)
        E_f = E_f0
        if not collapsed:
            for sw in range(2):
                sv.refine_sweep()
            E_f = sv.rebuild()
            collapsed = bool(E_f < E_nr_b - 10 * abs(first) - 1.0)
        rec = dict(band=bn, m_band=m_b, E_NR_mA=E_nr_A, E_NR_mband=E_nr_b, delta_first_order=first,
                   P4_nr_state=p_nr["P4"], mass_velocity_nr_state=mv1, ratio_first_order=first / mv1,
                   G4_frac_beyond_zone_nr_state=Tb / pf1["T"], variational_E_full_same_basis=E_f0, variational_E_full=E_f,
                   variational_collapsed=collapsed, sec=time.time() - t0)
        if name == "one_Z10":
            rec["ratio_vs_5Z4_first_order"] = first / (-5 * (7.2973525693e-3 ** 2) * 10 ** 4 * m_b / 8)
        res = load_json(JS); prod = res.get("prod", {}); prod.setdefault("full", {})[key] = rec; res["prod"] = prod; save_json(JS, res)
        print(key, {k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in rec.items()}, flush=True)

if what == "dirac":
    # DIAGNOSI POST HOC (dichiarata): riferimento = cinetica di Dirac ESATTA sqrt(k^2 c^2 + m^2 c^4) - m c^2 (stessa m della
    # banda) sullo stesso stato NR. Domanda: la banda imparata segue Dirac? (il rapporto con -alpha^2 <p^4>/8 non e' 1
    # nemmeno per Dirac, perche' <p^4> prende peso da p ~ c nella coda della cuspide)
    a = 7.2973525693e-3; c = 1 / a
    u = np.linspace(0, 14, 7001)[1:]
    wu = np.sqrt(2 / np.pi) * u ** 2 * np.exp(-u ** 2 / 2) * (u[1] - u[0])
    out = prod.setdefault("dirac_ref_posthoc", {})
    for bn in ("existing_A", "relearn_A"):
        Ftab, m_b, _ = band(bn)
        sv_ = Ftab[0]
        Fd = np.array([np.sum(wu * (np.sqrt((u * np.sqrt(s)) ** 2 * c * c + m_b ** 2 * c ** 4) - m_b * c * c)) for s in sv_])
        Rd = Fd / (1.5 * sv_)
        for name in ("one_Z10", "he1S_Z2", "he1S_Z6", "he1S_Z7", "he1S_Z8", "he1S_Z9", "he1S_Z10"):
            j = JOBS[name]
            s_nr = System(j["N"], j["nup"], j["Z"], m_A, name=name)
            s_d = System(j["N"], j["nup"], j["Z"], m_A, kin="full", Ftab=(sv_, Rd), name=name)
            sv = load_basis(name, j, s_nr)
            E_nr_b = sv.E * m_b / m_A
            dD = expect(sv, s_d)["E"] - E_nr_b
            fb = prod["full"][f"{bn}:{name}"]["delta_first_order"]
            out[f"{bn}:{name}"] = dict(delta_dirac_first_order=dD, delta_band_first_order=fb, band_over_dirac=fb / dD)
            print(bn, name, f"Dirac {dD:.6g}  banda {fb:.6g}  banda/Dirac {fb / dD:.4f}", flush=True)
    res = load_json(JS); res.setdefault("prod", {})["dirac_ref_posthoc"] = out; save_json(JS, res)
