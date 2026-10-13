"""lvloop_verdetto.py - apply the pre-registered conditional table (lvloop_prereg.json) to the phase-3 numbers of
vuoto_ (origin/master 35afcc5, lab/out/vuoto_f3_exams.json, read with git show and passed as argument).
Usage: git show origin/master:lab/out/vuoto_f3_exams.json > /tmp/f3.json; python lvloop_verdetto.py /tmp/f3.json
-> lab/out/lvloop_verdetto.json"""
import os, sys, json
import lvloop_prereg as P

f3 = json.load(open(sys.argv[1]))["table"]
B = P.BOUNDS
best = {
    "aniso_e": B["aniso_electron"][1]["B"],             # Yb+ 2022, 1e-21
    "aniso_e_2019": B["aniso_electron"][0]["B"],        # Yb+ 2019, 1.6e-20
    "iso_e_faster_astro": B["iso_electron_faster"][0]["B"],   # 2.5e-20 (Crab, indicative)
    "iso_g_faster_astro": B["iso_photon_faster"][0]["B"],     # 6e-21 (photon decay, showers)
    "iso_two_sided_LEP": B["iso_two_sided"][0]["B"],          # 5e-15 (LEP synchrotron, robust)
}
out = {"source": "origin/master 35afcc5 lab/out/vuoto_f3_exams.json", "alpha": P.ALPHA, "bounds": best, "cases": {}}
for key, label in (("R|K12", "regola congelata, taglio K12"), ("R|K12/3", "regola congelata, taglio K12/3"),
                   ("C|K12/3", "ordine a giro (proposta), taglio K12/3")):
    t = f3[key]
    an0, an1 = t["aniso_fit"][0], t["aniso_fit"][1]
    is0, is1 = t["iso_fit"][0], t["iso_fit"][1]
    an_abs, is_abs = abs(an0) * P.ALPHA, abs(is0) * P.ALPHA
    c = {"label": label, "aniso_a0_over_alpha": an0, "aniso_slope": an1, "iso_a0_over_alpha": is0, "iso_slope": is1,
         "scaling_exponent": 0, "aniso_abs": an_abs, "iso_abs": is_abs,
         "iso_sign": "elettrone piu' veloce della luce nuda (Cherenkov nel vuoto)" if is0 > 0 else
                     "elettrone piu' lento della luce nuda (decadimento del fotone)"}
    if an_abs > 1e-10:   # cyclic order: anisotropy is zero by symmetry (fit 5.6e-13 alpha = numerical zero)
        c["aniso_over_bound_Yb2022"] = an_abs / best["aniso_e"]
        c["aniso_over_bound_Yb2019"] = an_abs / best["aniso_e_2019"]
        c["aniso_verdict"] = "p = 0: ESCLUSA per qualunque a (pre-registrato R2)"
    else:
        c["aniso_verdict"] = "zero per simmetria (T_h): nessun limite"
    astro = best["iso_e_faster_astro"] if is0 > 0 else best["iso_g_faster_astro"]
    c["iso_over_bound_astro"] = is_abs / astro
    c["iso_over_bound_LEP"] = is_abs / best["iso_two_sided_LEP"]
    c["iso_verdict"] = ("p = 0: ESCLUSA per qualunque a CON RISERVA (misurato contro la luce NUDA): per salvarla il "
                        "vestito della luce deve compensare quello dell'elettrone entro 1 parte su %.0e (astro) / %.0e (LEP)"
                        % (is_abs / astro, is_abs / best["iso_two_sided_LEP"]))
    out["cases"][key] = c
    print(label, json.dumps({k: v for k, v in c.items() if k != "label"}, indent=0, default=float))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "lvloop_verdetto.json"), "w"),
          indent=1, ensure_ascii=False, default=float)
