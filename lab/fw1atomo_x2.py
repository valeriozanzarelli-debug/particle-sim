"""fw1atomo_x2.py - exploration AFTER the pre-registration (declared): size of the box effect in P3.  P3 uses a periodic
box of ~16 lambda_C; the 1s tail beyond 8 lambda_C is ~3 % and sees the periodic images of the nucleus (more binding).
Same Ritz as P3 (FW1, start 'dirac' on sublattice 0) with a box of ~20 lambda_C: a = 0,3 (n 66) and a = 0,24 (n 84).
Run: OMP_NUM_THREADS=1 timeout 590 python lab/fw1atomo_x2.py <a>  (exit 3 = checkpoint); sum -> lab/out/fw1atomo_x2.json"""
import os, sys, json
sys.path.insert(0, os.path.dirname(__file__))
import fw1atomo_p3 as p3

BOX = {"0.3": 66, "0.24": 84, "0.3b24": 80, "0.24b24": 100, "0.2b20": 100, "0.4b24": 60}
p3.AS = {0.3: 66, 0.24: 84}
p3.RES = os.path.join(p3.fc.SCR, "x2")
os.makedirs(p3.RES, exist_ok=True)
if __name__ == "__main__":
    if sys.argv[1] == "sum":
        out = {}
        for key in BOX:
            a = float(key.split("b")[0])
            d_ = p3.RES if "b" not in key else os.path.join(p3.fc.SCR, "x2_" + key)
            f = os.path.join(d_, f"ritz_FW1_dirac_{a}_c0.json")
            if os.path.exists(f):
                lv = [l for l in json.load(open(f))["levels"] if l["resid"] < 1e-5 and l["w"] > 0.01]   # 1e-5: at a = 0,2 the partner has 1,6e-6
                lv = sorted(lv, key=lambda l: -l["w"])[:2]
                om = sum(l["omega"] for l in lv) / len(lv)
                out[key] = {"n": BOX[key], "box": BOX[key] * a, "pair": [l["omega"] for l in lv], "B": 1 - om, "err": (1 - om) / p3.B_DIRAC - 1}
        json.dump(out, open(os.path.join(os.path.dirname(__file__), "out", "fw1atomo_x2.json"), "w"), indent=1)
        print(out)
    elif "b" in sys.argv[1]:
        # second set (box 24 at a = 0,3 / 0,24, box 20 at a = 0,2), Ritz grid dw = 0,02 for memory
        a = float(sys.argv[1].split("b")[0])
        p3.AS = {a: BOX[sys.argv[1]]}
        p3.RES = os.path.join(p3.fc.SCR, "x2_" + sys.argv[1])
        os.makedirs(p3.RES, exist_ok=True)
        orig = p3.fc.ritz
        p3.fc.ritz = lambda R, S, a_, lo, hi, dw, nper, **k: orig(R, S, a_, lo, hi, 0.02, nper, **k)
        p3.run(f"ritz_FW1_dirac_{a}_c0")
    else:
        p3.run(f"ritz_FW1_dirac_{sys.argv[1]}_c0")
