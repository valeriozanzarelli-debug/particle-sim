"""emiss6_prereg.py - pre-registration v6 (phase 4 of docs/EMISSIONE.md, §13): absorption and spontaneous emission
from the time evolution of the zig-zag rule, true charge K = 1.  Written BEFORE the production runs.
-> lab/out/emiss6_prereg.json"""
import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
G = (1 - (60 * 7.2973525693e-3)**2)**0.5
P = {
 "study": "emiss6_ (EMISSIONE.md fase 4)",
 "date": "2026-11-02",
 "rule": "electron: radial form of the zig-zag rule (lab/emiss_radial.py, no valleys); light: 3-D spin-1 zig-zag light walk "
         "128^3, sponge 12 cells (e^-7.5/pass), Gaussian average sigma = 1 cell, 30 light cells per wavelength; light entered "
         "through its exact linear response at the atom (lab/emiss6_kernel.py, 700 steps, two parities); coupling "
         "centered + transverse as production (lab/emiss_run.py) WITH AMENDMENT A6-1 (transverse projector of the deposit's "
         "own tick parity, transverse_odd); K = 1 (e = sqrt(4 pi alpha)); nucleus point charge Z; no fitted parameter.",
 "amendment_A6_1": "found in development (declared): Light3D.make_transverse uses the static axis of the two-tick light walk "
         "starting on an EVEN tick for every deposit; the walk alternates the axis order, so at odd light-tick counts the static "
         "mode is S0 n_even. Production kernel: odd-parity deposits leave a permanent field at the atom (|K(700)|/|K(1)| = 0.40, "
         "anisotropic x/y/z first response 9.1/7.2/4.6e-3 vs 7.0e-3); parity-correct kernel: 1.2e-3, isotropic. A6-1 = use "
         "T_odd = S0 P_even S0^T (w e_c) for deposits at odd counts. Same rule, no new parameter.",
 "development_runs_declared": [
   "dev_q60: Z=60 a=0.1 QE with production (blind) kernel to t=6000: windowed rate/A_rule W2 = 1.008, 1.002, 1.014, 0.96, then 0.21 "
   "after step ~700 (kernel cut of the permanent odd-parity field): the blind coupling cannot be cut, the parity one can.",
   "dev_r60: Z=60 a=0.1 cw Rabi: Omega_R/w = 0.02/0.01/0.005 -> W sqrt(P)/Omega_rule = 0.9912/0.9979/0.9995 (error ~ Omega^2): choose 0.005; "
   "P_max at Delta = +-Omega: 0.489/0.511.",
   "dev_pulse/packet: Z=60 a=0.1: analytic weak pulse R_dyn/R = 0.99996; light-walk packet helicity +1: 0.99990, helicity -1: 0.889 "
   "(wrong helicity: resonant component 1/0.19 of the counter-rotating one) -> helicity +1."],
 "configs": {
   "Z60": {"rmax": 40.0, "a": [0.1, 0.05, 0.025], "QE_T": "6/(A_rule W2); a=0.025: 10/(A_rule W2)",
           "rabi": "cw circular, Omega_R = w/200, 1.3 generalized cycles", "packet": "helicity +1, k = 2pi/30 cells^-1, "
           "Gaussian width 15 cells, start -35 cells, 240 light steps, amplitude for P_e = 1e-4"},
   "Z1": {"rmax": 2500.0, "a": [0.5, 0.25, 0.125], "QE_T": "6 optical periods, rate window [2, 6] periods",
          "packet": "as Z60"},
   "live3D": "Rust test emissione6_rabi_vivo: Esperimento::new(600) (Z=60, 48^3, a=0.5, learned nucleus, 1s in valley k=0), "
             "drive e E0 cos(w t) z at the strongest exact 2p level of livelli_di_dipolo(1s, z, 0.96-0.99, 1200 ticks), "
             "Omega_R = 2e-4, 1.3 cycles"},
 "estimators": {
   "A_fit": "least-squares slope of ln(norm) (norm = excited-world norm of the one-excitation state)",
   "omega_emit": "minus the slope of the unwrapped phase of D(t) = <1s(t)|r|psi(t)> . VEC* over the rate window",
   "rabi": "P_e = P sin^2(W t/2) least squares; Omega_meas = W sqrt(P); R_dyn = Omega_meas sqrt6 / (2 e E0)",
   "packet": "R_dyn = sqrt(P_e) / (e |E_c(w)|), E_c(w) = sum_s b E(s).VEC* e^{i w t_s} sinc(w b/2pi), E(s) the light walk's field at the atom",
   "f": "f(1s->2p3/2) = (2/3) w (4/6) R^2; Z=1 total f = (2/3) w R^2 (2p1/2 and 2p3/2 degenerate to 3e-5 at Z=1: declared)",
   "line_shape_natural": "power spectrum of J(t) over the whole a=0.025 Z=60 run, zero padding x8, Lorentzian fit within +-5 FWHM",
   "extrapolation_Z1": "c0 + c1 a^2 + c2 a^4 through the three a",
   "extrapolation_Z60_info": f"c0 + c1 a^(2 gamma) + c2 a^2, 2 gamma = {2*G:.4f}"},
 "exams": [
   {"id": "G0", "what": "kernel = full 128^3 box (production blind coupling), Z=60 a=0.1, T=1500: max|norm_k - norm_box|/(1-norm_box)", "threshold": "< 1e-3", "predicted": "PASS"},
   {"id": "GK", "what": "parity kernel decays: max |K(700)|/|K(1)|", "threshold": "< 2e-3", "predicted": "PASS (dev 1.2e-3)"},
   {"id": "X1", "what": "Z=60 a=0.025: exponential: max_t |norm - exp(-A_fit t)|, A_fit over [0.5/A_rW, T]", "threshold": "<= 0.005", "predicted": "PASS"},
   {"id": "X2", "what": "Z=60, each a: |A_fit/(A_rule W2) - 1|", "threshold": "<= 0.01", "predicted": "PASS"},
   {"id": "X3", "what": "Z=60 a=0.025: A_fit/W2 vs Dirac A_LW(2p3/2) = 1.08865e-5 (8.4517e15 s^-1)", "threshold": "<= 1.5 %", "predicted": "PASS (~+0.9 %)"},
   {"id": "X4", "what": "Z=60 a=0.025: max |p_2p - norm|", "threshold": "<= 2e-3", "predicted": "PASS"},
   {"id": "X5", "what": "Z=60 a=0.025: |omega_emit/omega_rule - 1| over [0.5, 3]/A", "threshold": "<= 2e-4", "predicted": "PASS"},
   {"id": "X6", "what": "Z=60 a=0.025: natural line: |FWHM/A_fit - 1| and max Lorentzian residual/peak", "threshold": "<= 0.05 and <= 0.05", "predicted": "PASS"},
   {"id": "B1", "what": "Z=60 each a: |Omega_meas/Omega_rule - 1| (cw Rabi)", "threshold": "<= 2e-3", "predicted": "PASS"},
   {"id": "B3", "what": "Z=60 a=0.025: f(1s->2p3/2) from the light-walk packet vs Dirac 0.25365", "threshold": "<= 1.5 %", "predicted": "PASS"},
   {"id": "B4", "what": "Z=60 a=0.05: power-broadened line: |P_max(D) - 1/(1+(D/Omega)^2)| for D/Omega = +-0.5, +-1, +-2", "threshold": "<= 0.02", "predicted": "PASS"},
   {"id": "B7", "what": "packet vs rule, every Z and a: |R_dyn/R_rule - 1|", "threshold": "<= 1e-3", "predicted": "PASS"},
   {"id": "H1", "what": "Z=1 each a: |A_fit/(A_rule W2) - 1|", "threshold": "<= 0.01", "predicted": "PASS"},
   {"id": "H2", "what": "Z=1: extrapolated A/W2 vs Dirac A_LW(2p3/2) = 6.26838e8 s^-1", "threshold": "<= 0.5 %", "predicted": "PASS (risk: a=0.5 not in a^2 regime)"},
   {"id": "H3", "what": "Z=1: |omega_emit/omega_rule - 1| each a; extrapolated lambda vs Dirac 121.50025 nm (2p3/2, infinite mass)", "threshold": "<= 1e-4 and <= 0.02 %", "predicted": "PASS"},
   {"id": "B2", "what": "Z=1: f(1s->2p) total from the packet, extrapolated, vs 0.4162", "threshold": "<= 0.5 %", "predicted": "PASS"},
   {"id": "B6", "what": "Z=1: line-centre cross-section 4 pi alpha f / A (from the two dynamic numbers) vs 3 lambda^2/(2 pi)", "threshold": "<= 1 %", "predicted": "PASS"},
   {"id": "H4", "what": "Z=1 info: with the proton mass (mu/m = 0.99945568): A mu/m and lambda m/mu vs NIST 6.2649e8 s^-1, 121.567 nm", "threshold": "info", "predicted": "-"},
   {"id": "V1", "what": "live 3-D atom: f(1s->2p3/2) = 2 w |d_z|^2 from the Rabi frequency vs Dirac 0.25365", "threshold": "<= 5 %", "predicted": "FAIL (cell 0.5 and box 24 lambda_C: dipoles 1.4-1.5x Dirac in phase 2)"},
   {"id": "V2", "what": "live 3-D atom: Omega_meas / (e E0 |d_z static|)", "threshold": "within 3 %", "predicted": "PASS"},
   {"id": "C5", "what": "cascade uses the dynamic A(2p->1s) when assets/dinamica/lyman_alpha.json exists (unit test)", "threshold": "exact", "predicted": "PASS"}],
 "notes": "Semiclassical absorption (coherent packet, first order, back-reaction off) and one-excitation quantum emission (debts of EMISSIONE §4.2). "
          "Radial channel cut: sigma+ from 1s m=+1/2 reaches only the p3/2 m=3/2 channel; d channels cut (declared). W2 = Gaussian-average form "
          "factor (computed, not fitted)."}
json.dump(P, open(os.path.join(OUT, "emiss6_prereg.json"), "w"), indent=1, ensure_ascii=False)
print("written")
