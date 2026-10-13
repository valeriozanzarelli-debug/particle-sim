"""surf_f4_prereg.py -> lab/out/surf_f4_prereg.json : PHASE 4 of SURFISTI_ESATTI, written BEFORE the runs of surf_f4.py."""
import json
import os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = dict(
    study="SURFISTI_ESATTI fase 4", date="2026-10-13",
    question="phase 3 found the Dirac rotation of the valley-0 1s at 0.85-0.87 Za c sin(theta) at a = 0.5: does the "
             "missing ~13 % go down with a, and with what law?",
    setup=dict(
        lattices=[dict(a=0.5, n=48), dict(a=0.4, n=60), dict(a=0.3, n=80), dict(a=0.25, n=96)],
        box="24 lambda_C (same as the Laboratorio at a = 0.5), Z = 60, nucleus on the corner of 8 cells",
        nucleus="NUCLEO_1S sec. 14 transferred nucleus: s(a) = s(0.5) + (s(0.3) - s(0.5)) (a - 0.5)/(0.3 - 0.5), "
                "s(0.5) = (0.554593, 1.003549), s(0.3) = (0.57482, 1.002538), factors on the 8 cells at r0 = sqrt3/2 a and the "
                "24 cells at r1 = sqrt11/2 a. a = 0.4 and 0.25 were never used to learn.",
        preparation="as phase 3, rescaled to the same physical times: seed = radial Dirac 1s (spin up); pass 0: "
                    "autocorrelation spectrum over 1024/a ticks (= 2048 at a = 0.5), Hann, peak in 0.86..0.94; passes 1 "
                    "and 2: filtering at the peak (1024/a ticks each, pass 2 from the pass-1 state); projection on valley "
                    "k = 0 (same Fourier cells as phase 2/3), normalised. Change from phase 3, declared: pass 1 of phase 3 "
                    "accumulated 33 candidate frequencies at once; here the peak comes from the spectrum alone (memory at "
                    "96^3), then 2 filtering passes.",
        run="physical time 2000 lambda_C/c (4000 ticks at a = 0.5, 8000 at a = 0.25), 2e5 exact surfers (two-tick flow), Born "
            "start; 5 surfer windows of ~5 lambda_C/c (10 a/c at a = 0.5) starting at t = 250, 500, 1000, 1500 lambda_C/c and "
            "one closing at the end of the run (the 3990 window of phase 3, now measured); records every 50 lambda_C/c.",
        v_of_a="v(a) = mean of the surfer v_phi/(Za sin) of the two early windows (250, 500), before the slow valley beat; "
               "deficit d(a) = 1 - v(a)"),
    exams=dict(
        G1=dict(what="every a: prepared state valley-0 weight >= 0.99 and <Sigma_z> >= 0.9", thr=[0.99, 0.9]),
        G2=dict(what="every a: fidelity >= 0.9 and valley-0 weight >= 0.9 at every record and at the end", thr=0.9),
        G3=dict(what="every a: surfer v_phi/(Za sin) in [0.75, 1.25] in all 5 windows (including the end window)",
                thr=[0.75, 1.25]),
        G4=dict(what="every a: max |z| <= 4", thr=4.0),
        G5=dict(what="the deficit goes down with a: d(0.5) > d(0.4) > d(0.3) > d(0.25)"),
        G6=dict(what="continuum limit: least-squares d(a) = d0 + c1 a on the 4 lattices, |d0| <= 0.05 (Dirac recovered "
                     "within 5 %); also reported: d = C a^p (p, R^2) and d0 + c2 a^2", thr=0.05)),
)
json.dump(P, open(os.path.join(OUT, "surf_f4_prereg.json"), "w"), indent=1)
print("written lab/out/surf_f4_prereg.json")
