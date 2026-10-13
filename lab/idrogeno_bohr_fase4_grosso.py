"""Idrogeno di Bohr, fase 4 (docs/IDROGENO_BOHR.md §13): il GROSSO DI DEFAULT (emendamento dichiarato, prereg
lab/out/idrogeno_bohr_fase4_prereg.json, esame G).

Livello grosso di default = regola T imparata dal maestro scelto dall'agglomeratore (lab/out/idrogeno_bohr_fase4_regola.json,
livello a0/8) + bordi ritorti (media theta = 0 e pi su tutti gli assi) + coda di Born sulla cella del nucleo (0 parametri).
Motore: lab/agg2s_core.py IN SOLA LETTURA (campo di Galerkin senza alias, TwistAtom, LOBPCG), stessa configurazione di
AGGLOMERATO_2S: a0/8, N = 256 (L = 32 a0). Il dV del primo ordine dell'accoppiamento (spin-orbita + Darwin imparati) e' letto
da aggfull_atomK.json come in agg2s_analyze.py.

Uso: timeout 590 python3 lab/idrogeno_bohr_fase4_grosso.py JOB [JOB ...]   (JOB = 1s_000, 1s_111, 2s_000, 2s_111, 2p_000,
     2p_111; a pezzi: un job non converso salva il vettore e riprende alla chiamata dopo)
     python3 lab/idrogeno_bohr_fase4_grosso.py esame   (criteri G1-G3 -> lab/out/idrogeno_bohr_fase4_grosso_esame.json)
Vettori e campo in lab/out/idrogeno_bohr_fase4_grosso_*.npy (mai committati)."""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import agg2s_core as c
from zz_core import dirac_formula

OUT = c.OUT
N = 256
ASSET = os.path.join(OUT, "idrogeno_bohr_fase4_regola.json")
JS = os.path.join(OUT, "idrogeno_bohr_fase4_grosso.json")
TOL = 2e-8
res = json.load(open(JS)) if os.path.exists(JS) else {}
JOBS = {f"{s}_{t}": (s, t, (f"1s_{t}" if s == "2s" else None)) for s in ("1s", "2s", "2p") for t in ("000", "111")}


def livello():
    a = json.load(open(ASSET))
    lv = [l for l in a["chain_A_3d"]["levels"] if abs(l["cell_over_a0"] - 0.125) < 1e-9][0]
    G = lv["generator_G"]
    K = np.array(G["re"]) + 1j * np.array(G["im"])          # (2, 2, n_offsets)
    return {"K": K, "offs": np.array(lv["offsets"])}, lv["cell_lambdaC"], a["maestro_j"]


def taylor_iso(rec, h):
    K, offs = rec["K"], np.asarray(rec["offs"], float)
    s = (0.5 * (K[0, 0] + K[1, 1])).real
    j2 = np.sum(offs ** 2, 1)
    inv2m = -(h / 6) * np.sum(s * j2)
    return float(1 / (2 * inv2m))


def campo(h):
    fn = os.path.join(OUT, "idrogeno_bohr_fase4_grosso_campo.npy")
    if os.path.exists(fn):
        W = np.load(fn)
        cf = c.CompactField.__new__(c.CompactField)
        cf.N, cf.h, cf.x = N, h, (np.arange(N) - N // 2) * h
        cf.subs = [(sx, sy, sz) for sx in (0, 1) for sy in (0, 1) for sz in (0, 1)]
        cf.W, cf.parity_err, cf.idx = W, 0.0, {}
        return cf
    cf = c.galerkin_field(N, h, cache=False)
    np.save(fn, cf.W)
    print(f"  campo di Galerkin N={N}: errore di parita' {cf.parity_err:.1e}", flush=True)
    return cf


def vecfile(job):
    return os.path.join(OUT, f"idrogeno_bohr_fase4_grosso_vec_{job}.npy")


def run(job, rec, h, f):
    kind, tws, dep = JOBS[job]
    if res.get(job, {}).get("converged"):
        print(f"{job}: gia' fatto", flush=True)
        return
    t0 = time.time()
    tw = tuple(int(x) for x in tws)
    T, spin = c.kin_learned(rec, N, h, tw)
    m_kin = taylor_iso(rec, h)
    at = c.TwistAtom(T, f, tw, gamma0=1.0, sym="oddz" if kind == "2p" else "even")
    Y = None
    if dep is not None:
        if not res.get(dep, {}).get("converged"):
            print(f"{job}: serve prima {dep}", flush=True)
            return
        Y = np.load(vecfile(dep))[:, None]
    x0 = np.load(vecfile(job)) if os.path.exists(vecfile(job)) else c.guess(at, kind)
    prev = res.get(job, {})
    v, E, r, it, ok = c.lobpcg1(at, x0, Y=Y, tol_res=TOL, budget=max(520 - (time.time() - t0) - 60, 60),
                                log=lambda s: print(f"  {job}{s}", flush=True))
    np.save(vecfile(job), v)
    vv = v.reshape((N,) * 3)
    h_a0 = h / c.A0
    nrm = float(np.sum(vv ** 2)) * h_a0 ** 3
    rho = float(vv[N // 2, N // 2, N // 2] ** 2 / nrm) if kind != "2p" else 0.0
    res[job] = dict(state=kind, twist=tws, N=N, h_a0=h_a0, box_a0=N * h_a0, E_eV=E * c.MC2, residual=r,
                    iterations=prev.get("iterations", 0) + it, converged=bool(ok), rho_c_a0=rho, m_kin=m_kin,
                    seconds=prev.get("seconds", 0) + time.time() - t0)
    json.dump(res, open(JS, "w"), indent=1)
    print(f"{job}: E = {E * c.MC2:.9f} eV, residuo {r:.2e}, converso {ok}, rho_c {rho:.6g}, m_kin {m_kin:.9f} "
          f"({time.time() - t0:.0f} s)", flush=True)


def esame():
    P = json.load(open(os.path.join(OUT, "idrogeno_bohr_fase4_prereg.json")))["esami"]["G_grosso_di_default"]["criteri"]
    AK = json.load(open(os.path.join(OUT, "aggfull_atomK.json")))
    BC = c.born_constants()
    MC2, RY = c.MC2, c.RY
    dV = dict(s1=float(np.mean(AK["qK8_1s"]["dV_eigen_eV"])), s2=float(np.mean(AK["qK8_2s"]["dV_eigen_eV"])),
              p_lo=float(np.mean(AK["qK8_2p"]["dV_eigen_eV"][:2])), p_hi=float(np.mean(AK["qK8_2p"]["dV_eigen_eV"][2:])))
    dirac = lambda n, k: (dirac_formula(n, k) - 1) * MC2
    avg = lambda s, k="E_eV": 0.5 * (res[f"{s}_000"][k] + res[f"{s}_111"][k])
    rel = lambda E, R: E / R - 1
    m = res["1s_000"]["m_kin"]
    h_a0 = res["1s_000"]["h_a0"]
    dEB = lambda rc: -2 * m * (4 * np.pi) ** 2 * rc * BC["J_coef"] * h_a0 ** 3 * 2 * RY
    rB = lambda rc: rc / (1 - 8 * np.pi * BC["I_coef"] * m * h_a0) ** 2
    out = dict(criteri=P, dV_eV=dV, born=BC, m_kin=m, livelli={}, esami={})
    L = {"1s1/2": (avg("1s") + dV["s1"], dEB(avg("1s", "rho_c_a0")), dirac(1, -1), 1),
         "2s1/2": (avg("2s") + dV["s2"], dEB(avg("2s", "rho_c_a0")), dirac(2, -1), 2),
         "2p1/2": (avg("2p") + dV["p_lo"], 0.0, dirac(2, 1), 2), "2p3/2": (avg("2p") + dV["p_hi"], 0.0, dirac(2, -2), 2)}
    worst = 0.0
    for k, (E, dB, D, n) in L.items():
        out["livelli"][k] = dict(E_senza_born_eV=E, dE_born_eV=dB, E_eV=E + dB, dirac_eV=D, rel=rel(E + dB, D),
                                 rel_senza_born=rel(E, D), rel_mkin_dirac_senza_born=rel(E, m * D),
                                 periodico_000_rel=rel(res[f"{k[:2]}_000"]["E_eV"] + (E - avg(k[:2])) + dB, D))
        worst = max(worst, abs(rel(E + dB, D)))
        print(f"{k}: {E + dB:.7f} eV vs Dirac {D:.7f}: {100 * rel(E + dB, D):+.5f} % (senza Born {100 * rel(E, D):+.5f} %, "
              f"vs m_kin x Dirac senza Born {100 * rel(E, m * D):+.5f} %)")
    out["esami"]["G1_livelli"] = dict(passed=bool(worst <= 5e-5), peggiore_rel=worst)
    r1 = rB(avg("1s", "rho_c_a0")) / (m ** 3 / np.pi)
    r2 = rB(avg("2s", "rho_c_a0")) / (m ** 3 / (8 * np.pi))
    out["esami"]["G2_nucleo"] = dict(passed=bool(abs(r1 - 1) <= 5e-3), rho_c_su_esatto_1s=avg("1s", "rho_c_a0") / (m ** 3 / np.pi),
                                     rho_B_su_esatto_1s=r1, rho_c_su_esatto_2s=avg("2s", "rho_c_a0") / (m ** 3 / (8 * np.pi)),
                                     rho_B_su_esatto_2s=r2)
    out["esami"]["G3_massa"] = dict(descrittivo=True, s_senza_born_vs_mkin_dirac={k: out["livelli"][k]["rel_mkin_dirac_senza_born"]
                                                                                  for k in ("1s1/2", "2s1/2")})
    for k, v in out["esami"].items():
        print(k, v)
    json.dump(out, open(os.path.join(OUT, "idrogeno_bohr_fase4_grosso_esame.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    if sys.argv[1] == "esame":
        esame()
    else:
        rec, h, j = livello()
        print(f"regola a0/8 del maestro j = {j}: h = {h:.6f} lambda_C (a0/{c.A0 / h:.3f}), m_kin {taylor_iso(rec, h):.9f}", flush=True)
        f = campo(h)
        for job in sys.argv[1:]:
            run(job, rec, h, f)
