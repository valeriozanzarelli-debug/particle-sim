"""mg_amend.py - DECLARED amendment 1 to mg_prereg.json, written before any R4 measurement.
The presets reuse the quality steps already in fotscena.rs (the pre-registered photon tests of docs/ATOMO.md are tied
to them): bassa = 8 cells/lambda 96^3, media = 12 / 128^3, alta = 24 / 192^3, ultra = 32 / 384^3 (new step).
Balls: 20k / 60k / 250k / 1M (not 4M: the balls still move on the CPU; 4M needs the balls on the GPU, not done).
The rule runs on the CPU in every preset for now: the GPU rule passed R1 but is not wired into the scene yet (debt).
R4 thresholds unchanged (bassa on the Intel UHD: >= 20 frames/s, <= 50 ms per tick)."""
import json, os
AM = {"emendamento": 1, "data": "2026-10-13", "motivo": __doc__,
      "preset": {"bassa": {"celle_per_lambda": 8, "n": 96, "palline": 20000, "regola": "cpu"},
                 "media": {"celle_per_lambda": 12, "n": 128, "palline": 60000, "regola": "cpu"},
                 "alta": {"celle_per_lambda": 24, "n": 192, "palline": 250000, "regola": "cpu"},
                 "ultra": {"celle_per_lambda": 32, "n": 384, "palline": 1000000, "regola": "cpu"}}}
if __name__ == "__main__":
    json.dump(AM, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "mg_amend.json"), "w"), indent=1, ensure_ascii=False)
