"""tav_run.py - corse di produzione (a pezzi, con checkpoint) dello studio TAVOLA (Li, Be).
Uso:  python lab/tav_run.py <job> [tmax_s]     riprende dal checkpoint lab/out/tav_<job>.npz, lavora fino a tmax_s
                                              (default 540 s), salva; ripetere finche' stampa DONE.
      python lab/tav_run.py list
Configurazioni: in lab/out/tav_prereg.json (lette da li', non da qui: la pre-registrazione comanda).
Risultati: lab/out/tav_runs.json (energie per K, giri, sottoinsiemi, valori d'attesa)."""
import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import tav_ecg as T

PRE = os.path.join(T.OUT, "tav_prereg.json")
RJ = os.path.join(T.OUT, "tav_runs.json")


def job_cfg(name):
    pre = json.load(open(PRE))
    j = dict(pre["svm_common"]); j.update(pre["jobs"][name])
    return j


def make_system(j):
    return T.System(j["Z"], j["spin"], L=j["L"], m_kin=T.M_T_THETA_A)


def ckpt_path(name):
    return os.path.join(T.OUT, f"tav_{name}.npz")


def save_ckpt(name, sv, st):
    tmp = ckpt_path(name) + ".tmp.npz"
    np.savez(tmp, A=sv.A, u=sv.u, st=json.dumps(st))
    os.replace(tmp, ckpt_path(name))


def run(name, tmax=540.0):
    j = job_cfg(name)
    t0 = time.time()
    s = make_system(j)
    st = dict(phase="grow", sweep=0, k=0, EK=[], sweeps=[], seed=j["seed"])
    rng = np.random.default_rng([j["seed"], 0])
    sv = T.SVM(s, rng, j["Z"], Lmax=j["Lmax"], trials=j["trials"], mut=j["mut"], p_pair=j["p_pair"],
               polish=j["polish"], npol=j["npol"], nm=j["nm"])
    if os.path.exists(ckpt_path(name)):
        z = np.load(ckpt_path(name))
        st = json.loads(str(z["st"]))
        sv.A = z["A"].copy(); sv.u = z["u"].copy()
        if len(sv.A):
            sv.rebuild()
        print(f"[{name}] ripreso: K={len(sv.A)} fase={st['phase']} giro={st['sweep']} k={st['k']} E={sv.E}", flush=True)
    if st["phase"] == "done":
        print("DONE", flush=True)
        return True
    # rng deterministico a ogni ripresa
    sv.rng = np.random.default_rng([j["seed"], len(sv.A), st["sweep"], st["k"]])
    last = time.time()
    while time.time() - t0 < tmax:
        if st["phase"] == "grow":
            if len(sv.A) >= j["K"]:
                st["phase"] = "sweep"; st["k"] = 0
                continue
            e = sv.add()
            if len(sv.A) % 25 == 0:
                st["EK"].append([len(sv.A), e])
                print(f"[{name}] K={len(sv.A)} E={e:.9f} ({time.time() - t0:.0f} s)", flush=True)
        elif st["phase"] == "sweep":
            if st["sweep"] >= j["sweeps"]:
                st["phase"] = "final"
                continue
            if st["k"] >= len(sv.A):
                sv.E = float(sv._eig(sv.Hmat, sv.Smat)[0][0])
                st["sweeps"].append(sv.E)
                print(f"[{name}] giro {st['sweep'] + 1}: E={sv.E:.9f} ({time.time() - t0:.0f} s)", flush=True)
                st["sweep"] += 1; st["k"] = 0
                continue
            sv.refine_one(st["k"]); st["k"] += 1
        elif st["phase"] == "final":
            E = sv.rebuild()
            K = len(sv.A)
            sub = {str(k): sv.energy_subset(k) for k in (int(0.8 * K), int(0.9 * K))}
            ex = sv.expect() if j["L"] == 0 else {}
            rec = dict(E_mT=E, m_T=T.M_T_THETA_A, m0=T.M0, E_theta0=E * T.M0 / T.M_T_THETA_A, K=K, EK=st["EK"],
                       sweeps=st["sweeps"], subsets=sub, expect=ex, cfg=j,
                       sec_last_chunk=time.time() - t0)
            # barra pre-registrata: max(|ultimo giro - penultimo|, 0,5 |E(K) - E(0,9 K)|)
            sw = st["sweeps"]
            b1 = abs(sw[-1] - sw[-2]) if len(sw) >= 2 else float("nan")
            rec["barra"] = float(max(b1, 0.5 * abs(E - sub[str(int(0.9 * K))])))
            res = T.load_json(RJ); res[name] = rec; T.save_json(RJ, res)
            st["phase"] = "done"
            save_ckpt(name, sv, st)
            print(f"[{name}] FINE E(mT)={E:.9f} E(theta->0)={rec['E_theta0']:.9f} barra={rec['barra']:.2e}", flush=True)
            print("DONE", flush=True)
            return True
        if time.time() - last > 60:
            save_ckpt(name, sv, st); last = time.time()
    save_ckpt(name, sv, st)
    print(f"[{name}] pausa: K={len(sv.A)} fase={st['phase']} giro={st['sweep']} k={st['k']}", flush=True)
    return False


if __name__ == "__main__":
    if sys.argv[1] == "list":
        print(list(json.load(open(PRE))["jobs"]))
    else:
        run(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 540.0)
