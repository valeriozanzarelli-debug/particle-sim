"""fp_amend.py - DECLARED amendment to fp_prereg.json (written before any exam was judged).
What happened: in the first production chunk (circ+, t = 0) the spin read DIRECTLY FROM THE FIELD was 0.8945 instead
of ~1. Diagnosis (scratch test, same pulse builder): the transverse box (6 lambda) is too narrow for w0 = 2 lambda; the
spin flow lives at the edge of the beam and the periodic images cancel part of it. Field-only values of J_z*omega/U:
box 6 lambda / w0 2 lambda: 0.894; box 12 lambda / w0 2 lambda: 1.0063; box 8 lambda / w0 1.5 lambda: 1.0012
(= k0 <1/k> of the pulse spectrum, 1.0012: the exact expectation). The balls were never judged. Exams and thresholds
unchanged; only the box and the waist change. The partial run was deleted."""
import json, os
AM = {
    "emendamento": 1,
    "data": "2026-10-12",
    "motivo": __doc__,
    "configurazioni_nuove": {"scatola_celle": [192, 192, 192], "w0_lambda": 1.5},
    "atteso_F2_dal_campo": "sigma * k0 <1/k> = +-1.0012 (spettro del pacchetto)",
}
if __name__ == "__main__":
    json.dump(AM, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "fp_amend.json"), "w"), indent=1, ensure_ascii=False)
