"""tav_basi.py - salva le basi ECG finali (lab/out/tav_<job>.npz, non tracciati) in JSON tracciabili lab/out/tav_basi/tav_<job>.json
(formato: ogni campo {dtype, shape, data piatto}; 'st' = stato del driver). Uso: python lab/tav_basi.py be1S be1P ..."""
import os, sys, json
import numpy as np
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(os.path.join(OUT, "tav_basi"), exist_ok=True)
for job in sys.argv[1:]:
    z = np.load(os.path.join(OUT, f"tav_{job}.npz"))
    d = {k: dict(dtype=str(z[k].dtype), shape=list(z[k].shape), data=z[k].ravel().tolist()) for k in ("A", "u")}
    d["st"] = json.loads(str(z["st"]))
    json.dump(d, open(os.path.join(OUT, "tav_basi", f"tav_{job}.json"), "w"))
    print("salvato", job, d["A"]["shape"], d["st"]["phase"])
